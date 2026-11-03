#############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
#############################################################################

from copy import deepcopy

import numpy as np
import pytest

from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_fixed_recovery import CDSFixedRecovery
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_vars import DIRTY, ONE_MILLION

CURVE_RECOVERY = 0.40
FIXED_RECOVERY = 0.20
ONE_BP = 1.0e-4
ONE_PCT = 0.01


@pytest.fixture
def market():
    """Bootstrapped issuer curve from ordinary CDS market quotes."""
    value_dt = Date(20, 6, 2018)
    step_in_dt = value_dt.add_days(1)

    times = np.arange(11, dtype=float)
    dates = value_dt.add_years(times)
    discount_factors = np.exp(-0.03 * times)
    discount_curve = DiscountCurve(
        value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES
    )

    maturities = [value_dt.next_cds_date(m) for m in (12, 36, 60)]
    quotes = [0.0080, 0.0100, 0.0120]
    calibration_contracts = [
        CDS(step_in_dt, maturity, spread)
        for maturity, spread in zip(maturities, quotes)
    ]
    curve = CDSCurve(
        value_dt, calibration_contracts, discount_curve, CURVE_RECOVERY
    )
    return value_dt, step_in_dt, maturities[-1], curve


def make_contracts(market, fixed_recovery=FIXED_RECOVERY,
                   coupon=0.0100, long_protect=True):
    _, step_in_dt, maturity_dt, _ = market
    fixed = CDSFixedRecovery(
        step_in_dt, maturity_dt, coupon,
        fixed_recovery=fixed_recovery,
        notional=ONE_MILLION,
        long_protect=long_protect,
    )
    standard = CDS(
        step_in_dt, maturity_dt, coupon,
        notional=ONE_MILLION,
        long_protect=long_protect,
    )
    return fixed, standard


def assert_close(actual, expected, rtol=1e-10, atol=1e-7):
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)


@pytest.mark.parametrize("fixed_recovery", [0.0, 0.20, 0.40, 1.0])
def test_valuation_matches_standard_cds(market, fixed_recovery):
    """Only the recovery applied to the protection payoff differs."""
    value_dt, _, _, curve = market
    fixed, standard = make_contracts(market, fixed_recovery=fixed_recovery)
    assert_close(
        fixed.value(value_dt, curve),
        standard.value(value_dt, curve, recovery_rate=fixed_recovery),
    )
    assert_close(
        fixed.prot_leg_pv(value_dt, curve),
        standard.prot_leg_pv(value_dt, curve, recovery_rate=fixed_recovery),
    )


def test_matching_recovery_argument_is_accepted(market):
    value_dt, _, _, curve = market
    fixed, _ = make_contracts(market)
    assert_close(
        fixed.value(value_dt, curve),
        fixed.value(value_dt, curve, recovery_rate=FIXED_RECOVERY),
    )
    assert_close(
        fixed.prot_leg_pv(value_dt, curve),
        fixed.prot_leg_pv(value_dt, curve, recovery_rate=FIXED_RECOVERY),
    )


@pytest.mark.parametrize("method", ["value", "prot_leg_pv"])
def test_conflicting_recovery_argument_is_rejected(market, method):
    value_dt, _, _, curve = market
    fixed, _ = make_contracts(market)
    with pytest.raises(FinError, match="[Rr]ecovery"):
        getattr(fixed, method)(value_dt, curve, recovery_rate=CURVE_RECOVERY)


@pytest.mark.parametrize("recovery", [-0.01, 1.01])
def test_invalid_fixed_recovery_is_rejected(market, recovery):
    _, step_in_dt, maturity_dt, _ = market
    with pytest.raises(FinError):
        CDSFixedRecovery(step_in_dt, maturity_dt, 0.01, recovery)


def test_premium_leg_does_not_depend_on_fixed_recovery(market):
    value_dt, _, _, curve = market
    low, _ = make_contracts(market, fixed_recovery=0.0)
    high, _ = make_contracts(market, fixed_recovery=0.80)
    assert_close(low.rpv01(value_dt, curve), high.rpv01(value_dt, curve))
    assert_close(
        low.premium_leg_pv(value_dt, curve),
        high.premium_leg_pv(value_dt, curve),
    )


def test_protection_pv_scales_with_loss_given_default(market):
    value_dt, _, _, curve = market
    cds_20, _ = make_contracts(market, fixed_recovery=0.20)
    cds_60, _ = make_contracts(market, fixed_recovery=0.60)
    pv_20 = cds_20.prot_leg_pv(value_dt, curve)
    pv_60 = cds_60.prot_leg_pv(value_dt, curve)
    assert pv_20 > 0.0
    assert_close(pv_60 / pv_20, (1.0 - 0.60) / (1.0 - 0.20))


def test_full_recovery_has_zero_protection_leg(market):
    value_dt, _, _, curve = market
    fixed, _ = make_contracts(market, fixed_recovery=1.0)
    assert_close(fixed.prot_leg_pv(value_dt, curve), 0.0)
    assert_close(fixed.par_spread(value_dt, curve, recovery_rate=1.0), 0.0)


def test_long_short_pv_symmetry(market):
    value_dt, _, _, curve = market
    long_cds, _ = make_contracts(market, long_protect=True)
    short_cds, _ = make_contracts(market, long_protect=False)
    assert_close(long_cds.value(value_dt, curve),
                 -np.asarray(short_cds.value(value_dt, curve)))


def test_clean_dirty_difference_matches_accrued_premium(market):
    value_dt, _, _, curve = market
    fixed, _ = make_contracts(market)
    dirty, clean = fixed.value(value_dt, curve)
    dirty_rpv01, clean_rpv01 = fixed.rpv01(value_dt, curve)
    assert dirty_rpv01 >= clean_rpv01
    assert_close(
        clean - dirty,
        fixed.running_cpn * fixed.notional * (dirty_rpv01 - clean_rpv01),
    )


def test_par_spread_matches_standard_cds(market):
    value_dt, _, _, curve = market
    fixed, standard = make_contracts(market)
    assert_close(
        fixed.par_spread(value_dt, curve, recovery_rate=FIXED_RECOVERY),
        standard.par_spread(value_dt, curve, recovery_rate=FIXED_RECOVERY),
    )


def test_inherited_spread_and_interest_dv01(market):
    """The inherited methods take the recovery argument explicitly."""
    value_dt, _, _, curve = market
    fixed, standard = make_contracts(market)
    assert_close(
        fixed.spread_dv01(value_dt, curve, recovery_rate=FIXED_RECOVERY),
        standard.spread_dv01(value_dt, curve, recovery_rate=FIXED_RECOVERY),
    )
    assert_close(
        fixed.ir_dv01(value_dt, curve, recovery_rate=FIXED_RECOVERY),
        standard.ir_dv01(value_dt, curve, recovery_rate=FIXED_RECOVERY),
    )


def test_curve_recovery_dv01_holds_fixed_recovery_constant(market):
    value_dt, _, _, curve = market
    fixed, _ = make_contracts(market)
    bumped_curve = deepcopy(curve)
    bumped_curve.recovery_rate += ONE_PCT
    bumped_curve.build_curve()

    expected = (
        fixed.value(value_dt, bumped_curve)[DIRTY]
        - fixed.value(value_dt, curve)[DIRTY]
    )
    actual = fixed.recovery_dv01(value_dt, curve)
    assert_close(actual, expected)
    assert abs(actual) > 1e-6  # Should not silently be zero.


def test_fast_approximation_is_consistent_with_fixed_recovery(market):
    value_dt, _, _, _ = market
    fixed, standard = make_contracts(market)
    rate, spread = 0.03, 0.012

    fixed_results = fixed.value_fast_approx(
        value_dt, rate, spread, curve_recovery=CURVE_RECOVERY
    )
    standard_results = standard.value_fast_approx(
        value_dt, rate, spread,
        curve_recovery=CURVE_RECOVERY,
        recovery_rate=FIXED_RECOVERY,
    )
    assert len(fixed_results) == 7
    # PV, RPV01, credit DV01 and IR DV01 agree; recovery DV01 must differ.
    assert_close(fixed_results[:6], standard_results[:6], rtol=1e-9)

    bumped = fixed.value_fast_approx(
        value_dt, rate, spread,
        curve_recovery=CURVE_RECOVERY + ONE_PCT,
    )
    assert_close(fixed_results[6], bumped[0] - fixed_results[0])


def test_fast_approximation_supports_zero_coupon(market):
    value_dt, _, _, _ = market
    fixed, _ = make_contracts(market, coupon=0.0)
    results = fixed.value_fast_approx(
        value_dt, 0.03, 0.012, curve_recovery=CURVE_RECOVERY
    )
    assert len(results) == 7
    assert np.all(np.isfinite(results))
    assert_close(results[0], results[1])  # No accrued premium.


def test_fixed_recovery_is_shown_in_repr(market):
    fixed, _ = make_contracts(market)
    assert "FIXED RECOVERY" in repr(fixed)
