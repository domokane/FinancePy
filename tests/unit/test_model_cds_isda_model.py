"""ISDA vs OT kernel regression. Run with pytest from FinancePy root.
No official C library, market grids, curve bootstrap or SciPy required.
Times ACT/365F, premiums ACT/360, protection PV per unit notional.
Requires OT fixes: inclusive-end Q and Z, honour default-accrual flag.
"""
import numpy as np
import pytest
from financepy.models import cds_isda_model as isda
from financepy.models.cds_okane_turnbull_model import (
    rpv01_ot_numba, prot_leg_pv_ot_numba,
)

from financepy.utils.global_vars import CLEAN, DIRTY
isda_rpv = getattr(isda, 'rpv01_isda_numba', None) or isda.rpv01_numba
isda_prot = getattr(isda, 'prot_leg_pv_isda_numba', None) or isda.prot_leg_pv_numba
DAY = 1.0 / 365.0
ATOL = 3e-12


def close(x, y, atol=ATOL):
    np.testing.assert_allclose(x, y, rtol=0.0, atol=atol)


def schedule(step=1, maturity=1825, accrued=43, delay=2):
    """Synthetic quarterly schedule with final maturity day included."""
    starts = [step-accrued]
    ends = []
    while starts[-1]+91 < maturity+1:
        ends.append(starts[-1]+91)
        starts.append(ends[-1])
    ends.append(maturity+1)
    a = np.array(starts, dtype=float)
    b = np.array(ends, dtype=float)
    p = b.copy()
    p[-1] = maturity+delay
    return a*DAY, b*DAY, p*DAY, (b-a)/360.


def flat(r, h):
    t = np.array([0., .1, 1., 3., 10., 40.])
    return t, np.exp(-r*t), t.copy(), np.exp(-h*t)


def shaped():
    dt = np.array([0., .13, .72, 1.4, 3.3, 7., 30.])
    qt = np.array([0., .21, .91, 2.2, 4.1, 8., 30.])
    z = np.exp(-np.r_[0., np.cumsum(np.diff(dt)*[.02, .045, -.005, .035, .06, .04])])
    q = np.exp(-np.r_[0., np.cumsum(np.diff(qt)*[.01, .035, .08, .02, .12, .05])])
    return dt, z, qt, q


def factor(t, nodes, values):
    i = np.clip(np.searchsorted(nodes, t, side='right')-1, 0, len(nodes)-2)
    return np.exp(np.log(values[i])+(t-nodes[i])*np.log(values[i+1]/values[i])/(nodes[i+1]-nodes[i]))


def premiums(curves, sched, step, accrued, df, aod):
    a, b, p, f = sched
    dt, z, qt, q = curves
    ot = rpv01_ot_numba((step-1)*DAY, accrued, p, f, dt, z, qt, q, a, b, df, aod)
    ix = isda_rpv(step*DAY, accrued, p, f, dt, z, qt, q, a, b, df, aod)
    return ot, ix


def protections(curves, step, maturity, recovery):
    dt, z, qt, q = curves
    return (prot_leg_pv_ot_numba((step-1)*DAY, maturity*DAY, dt, z, qt, q, recovery),
            isda_prot(step*DAY, maturity*DAY, dt, z, qt, q, recovery))


@pytest.mark.parametrize('r,h', [(0., 0.), (.04, 0.), (0., .03), (.04, .02), (-.03, .03), (-.05, .01), (.15, .5), (.03, 2.), (1e-10, 1e-10)])
@pytest.mark.parametrize('step,maturity', [(1, 1), (1, 30), (1, 365), (1, 1825), (91, 365), (365, 3650)])
@pytest.mark.parametrize('recovery', [0., .4, 1.])
def test_flat_protection(r, h, step, maturity, recovery):
    ot, ix = protections(flat(r, h), step, maturity, recovery)
    start = (step-1)*DAY
    length = maturity*DAY-start
    alpha = r+h
    moment = length if alpha == 0. else -np.expm1(-alpha*length)/alpha
    expected = (1-recovery)*h*np.exp(-alpha*start)*moment
    close(ot, expected)
    close(ix, expected)


@pytest.mark.parametrize('step,maturity,accrued,delay', [(1, 1, 0, 2), (1, 2, 89, 3), (1, 90, 0, 2), (1, 91, 1, 1), (1, 1825, 43, 2), (1, 3650, 90, 4), (91, 365, 0, 2), (365, 1825, 42, 3)])
@pytest.mark.parametrize('curve', ['flat', 'shaped', 'zero', 'negative'])
def test_scheduled_exact(step, maturity, accrued, delay, curve):
    curves = {'flat': flat(.04, .02), 'shaped': shaped(), 'zero': flat(0, 0), 'negative': flat(-.03, .03)}[curve]
    s = schedule(step, maturity, accrued, delay)
    a, b, p, f = s
    dt, z, qt, q = curves
    df = factor(3*DAY, dt, z)
    expected = np.array([factor(t-DAY, qt, q)*factor(pay, dt, z)*frac for t, pay, frac in zip(b, p, f)])
    (ov, ot), (iv, it) = premiums(curves, s, step, accrued/360, df, False)
    close(ov, iv)
    close(ov[DIRTY], sum(expected))
    close(ot, it)
    close(ot[:, 4], np.cumsum(expected))
    close(ov[DIRTY]-ov[CLEAN], df*accrued/360)


@pytest.mark.parametrize('r,h', [(0., 0.), (.04, .02), (-.03, .03), (.08, .15), (.03, 1.), (.15, 2.)])
@pytest.mark.parametrize('step,maturity,accrued', [(1, 1, 0), (1, 30, 43), (1, 1825, 43), (91, 365, 0), (365, 3650, 90)])
def test_flat_aod_difference(r, h, step, maturity, accrued):
    """Account explicitly for midpoint and half-day differences."""
    curves = flat(r, h)
    s = schedule(step, maturity, accrued)
    a, b, p, f = s
    (ov, _), (iv, _) = premiums(curves, s, step, accrued/360, .997, True)
    (on, _), (inn, _) = premiums(curves, s, step, accrued/360, .997, False)
    nodes, weights = np.polynomial.legendre.leggauss(32)
    expected_ot = 0.
    expected_ix = 0.
    for i in range(len(p)):
        left = max((step-1)*DAY, a[i]-DAY, 0.)
        right = b[i]-DAY
        if right <= left:
            continue
        x = left+(nodes+1)*(right-left)/2
        w = weights*(right-left)/2
        density = h*np.exp(-(r+h)*x)
        kappa = f[i]/(b[i]-a[i])
        expected_ix += np.sum(w*density*kappa*(x-a[i]+1.5*DAY))
        if i == 0:
            dq = np.exp(-h*left)-np.exp(-h*right)
            expected_ot += np.exp(-r*p[i])*dq*(accrued/360+.5*(f[i]-accrued/360))
        else:
            expected_ot += np.sum(w*density*kappa*(x-left))
    close(ov[DIRTY]-on[DIRTY], expected_ot)
    close(iv[DIRTY]-inn[DIRTY], expected_ix)
    close(ov-iv, expected_ot-expected_ix)
    assert ov[DIRTY] >= on[DIRTY]-ATOL
    assert iv[DIRTY] >= inn[DIRTY]-ATOL


@pytest.mark.parametrize('step,maturity,accrued', [(1, 30, 43), (1, 365, 0), (1, 1825, 90), (91, 1825, 0), (365, 3650, 42)])
def test_shaped_curves(step, maturity, accrued):
    curves = shaped()
    s = schedule(step, maturity, accrued)
    (ov, _), (iv, _) = premiums(curves, s, step, accrued/360, .999, True)
    (on, _), (inn, _) = premiums(curves, s, step, accrued/360, .999, False)
    close(on, inn)
    # Practical limits for these specific moderate piecewise curves.
    close(ov, iv, atol=.0015)
    op, ip = protections(curves, step, maturity, .4)
    close(op, ip, atol=2e-5)


@pytest.mark.parametrize('aod', [False, True])
@pytest.mark.parametrize('direction', [-1., 1.])
@pytest.mark.parametrize('coupon', [0., .01, .05])
def test_settlement_and_value_identities(aod, direction, coupon):
    curves = shaped()
    s = schedule()
    accrued = 43/360
    n = 1e6
    op, ip = protections(curves, 1, 1825, .4)
    old = premiums(curves, s, 1, accrued, 1., aod)
    new = premiums(curves, s, 1, accrued, .995, aod)
    for (v1, _), (v2, _), prot in zip(old, new, (op, ip)):
        close(v1[DIRTY], v2[DIRTY])
        close(v2[CLEAN]-v1[CLEAN], .005*accrued)
        dirty = direction*n*(prot-coupon*v2[DIRTY])
        clean = direction*n*(prot-coupon*v2[CLEAN])
        upfront = clean/(.995*n)
        csa = dirty/.995
        close(clean-dirty, direction*n*coupon*accrued*.995, 1e-7)
        close(csa, n*upfront-direction*n*coupon*accrued, 1e-7)
        if v2[CLEAN] > 0:
            par = prot/v2[CLEAN]
            close(n*(prot-par*v2[CLEAN]), 0., 1e-7)


@pytest.mark.parametrize('step', [0., -DAY, .5*DAY])
def test_isda_rejects_bad_stepin(step):
    dt, z, qt, q = flat(.04, .02)
    a, b, p, f = schedule()
    with pytest.raises(ValueError):
        isda_prot(step, 1., dt, z, qt, q, .4)
    with pytest.raises(ValueError):
        isda_rpv(step, 0., p, f, dt, z, qt, q, a, b, 1.)


@pytest.mark.parametrize('df', [0., -1., np.nan, np.inf])
def test_isda_rejects_bad_settlement(df):
    dt, z, qt, q = flat(.04, .02)
    a, b, p, f = schedule()
    with pytest.raises(ValueError):
        isda_rpv(DAY, 0., p, f, dt, z, qt, q, a, b, df)


def test_isda_empty_and_expired():
    dt, z, qt, q = flat(.04, .02)
    empty = np.array([])
    v, table = isda_rpv(DAY, 0., empty, empty, dt, z, qt, q, empty, empty, 1.)
    close(v, [0., 0.])
    assert table.shape == (0, 5)
    a, b, p, f = schedule(1, 30, 0)
    v, table = isda_rpv(31*DAY, 0., p, f, dt, z, qt, q, a, b, 1.)
    close(v, [0., 0.])
    assert table.shape == (0, 5)
    close(isda_prot(31*DAY, 30*DAY, dt, z, qt, q, .4), 0.)
