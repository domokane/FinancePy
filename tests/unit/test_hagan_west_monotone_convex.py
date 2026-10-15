import numpy as np
import pytest
from scipy.integrate import quad
from financepy.market.curves.hagan_west_monotone_convex import HaganWestMonotoneConvex
from financepy.market.curves.interpolator import Interpolator, InterpTypes
from financepy.utils.error import FinError


def test_hagan_west_matches_published_forward_example():
    times=np.arange(0.,6.)
    zero=np.array([0.,.0275,.04,.0475,.05,.0475])
    hw=HaganWestMonotoneConvex(times,np.exp(-times*zero))
    # Example from Hagan-West 2008, reproduced in TF Quant Finance.
    np.testing.assert_allclose(
        hw.forward_rate(np.array([.3,1.3,2.1,4.5])),
        [.0229375,.05010625,.0609,.03625],atol=1e-10
    )


@pytest.mark.parametrize("dfs",[
    [1.,.98,.95,.92,.87],
    [1.,1.01,1.03,1.05,1.06],
    [1.,.98,.94,.94,.91],
])
def test_preserves_knots_and_all_interval_integrals(dfs):
    times=np.array([0.,.25,.75,2.,5.])
    hw=HaganWestMonotoneConvex(times,dfs)
    np.testing.assert_allclose(hw.discount_factor(times),dfs,atol=1e-13,rtol=0)
    for i in range(1,len(times)):
        integrated=quad(hw.forward_rate,times[i-1],times[i],epsabs=1e-11)[0]
        assert integrated == pytest.approx(-np.log(dfs[i]/dfs[i-1]),abs=2e-8)


def test_flat_forward_curve_and_extrapolation():
    times=np.array([0.,.5,1.,2.,3.])
    dfs=np.exp(-.04*times)
    hw=HaganWestMonotoneConvex(times,dfs)
    t=np.linspace(0,4,101)
    np.testing.assert_allclose(hw.forward_rate(t),.04,atol=1e-12)
    np.testing.assert_allclose(hw.discount_factor(t),np.exp(-.04*t),atol=1e-12)


def test_interpolator_and_discount_curve_type():
    from financepy.market.curves.discount_curve import DiscountCurve
    from financepy.utils.date import Date
    t=np.array([0.,.4,1.,1.8,3.])
    dfs=np.array([1.,.99,.965,.93,.89])
    interp=Interpolator(InterpTypes.MONOTONE_CONVEX)
    interp.fit(t,dfs)
    hw=HaganWestMonotoneConvex(t,dfs)
    grid=np.linspace(0,4,50)
    np.testing.assert_allclose(interp.interpolate(grid),hw.discount_factor(grid))
    assert interp.interpolate(1.) == pytest.approx(dfs[2])
    anchor=Date(5,1,2026)
    curve=DiscountCurve(anchor,[anchor.add_years(1),anchor.add_years(2)],
                        np.exp(-.04*np.array([1.,2.])),
                        interp_type=InterpTypes.MONOTONE_CONVEX)
    assert curve.df(anchor)==pytest.approx(1.)
    assert curve.df(anchor.add_years(2))==pytest.approx(np.exp(-.08))


@pytest.mark.parametrize("times,dfs",[
    ([0,1,1],[1,.9,.8]),([0,2,1],[1,.9,.8]),
    ([0,1],[1,-.8]),([1,2],[.9,.8]),([0,1],[.99,.95]),
])
def test_invalid_curves_rejected(times,dfs):
    with pytest.raises(FinError):
        HaganWestMonotoneConvex(times,dfs)

