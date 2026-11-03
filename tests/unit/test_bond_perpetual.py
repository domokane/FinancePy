# Copyright (C) 2018-2026 Dominic O'Kane

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.zero_rates_discount_curve import ZeroRatesDiscountCurve
from financepy.products.bonds.bond import YTMCalcType
from financepy.products.bonds.bond_perpetual import BondPerpetual
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.error import FinError
from financepy.utils.frequency import FrequencyTypes


def test_perpetual_bond_pricing_on_coupon_date():
    issue_dt = Date(15, 1, 1923)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA

    bond = BondPerpetual(issue_dt, coupon, freq_type, dc_type)
    settle_dt = Date(15, 7, 2026)

    # Par pricing
    y_par = 0.05
    dp = bond.dirty_price_from_ytm(settle_dt, y_par)
    cp = bond.clean_price_from_ytm(settle_dt, y_par)
    ai = bond.accrued_interest(settle_dt)

    assert dp == pytest.approx(100.0, rel=1e-10)
    assert cp == pytest.approx(100.0, rel=1e-10)
    assert ai == pytest.approx(0.0, abs=1e-10)

    # Inversion: clean price -> YTM
    ytm_calc = bond.yield_to_maturity(settle_dt, cp)
    assert ytm_calc == pytest.approx(y_par, abs=1e-7)

    # Premium pricing (y < c)
    y_prem = 0.04
    dp_prem = bond.dirty_price_from_ytm(settle_dt, y_prem)
    assert dp_prem == pytest.approx(100.0 * (0.05 / 0.04), rel=1e-10)
    ytm_prem = bond.yield_to_maturity(settle_dt, dp_prem)
    assert ytm_prem == pytest.approx(y_prem, abs=1e-7)

    # Discount pricing (y > c)
    y_disc = 0.08
    dp_disc = bond.dirty_price_from_ytm(settle_dt, y_disc)
    assert dp_disc == pytest.approx(100.0 * (0.05 / 0.08), rel=1e-10)
    ytm_disc = bond.yield_to_maturity(settle_dt, dp_disc)
    assert ytm_disc == pytest.approx(y_disc, abs=1e-7)


def test_perpetual_bond_pricing_between_coupon_dates():
    issue_dt = Date(15, 1, 1923)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA

    bond = BondPerpetual(issue_dt, coupon, freq_type, dc_type)
    # Settle midway through the semi-annual period
    settle_dt = Date(15, 10, 2026)
    y = 0.05

    dp = bond.dirty_price_from_ytm(settle_dt, y, convention=YTMCalcType.UK_DMO)
    cp = bond.clean_price_from_ytm(settle_dt, y, convention=YTMCalcType.UK_DMO)
    ai = bond.accrued_interest(settle_dt)

    assert dp == pytest.approx(cp + ai, rel=1e-10)
    assert ai > 0.0

    # Invert YTM from clean price
    ytm_calc = bond.yield_to_maturity(settle_dt, cp, convention=YTMCalcType.UK_DMO)
    assert ytm_calc == pytest.approx(y, abs=1e-7)

    # Test US_TREASURY and US_STREET conventions
    dp_ust = bond.dirty_price_from_ytm(settle_dt, y, convention=YTMCalcType.US_TREASURY)
    dp_str = bond.dirty_price_from_ytm(settle_dt, y, convention=YTMCalcType.US_STREET)
    assert dp_ust > 0.0
    assert dp_str > 0.0


def test_perpetual_bond_duration_and_convexity():
    issue_dt = Date(15, 1, 1923)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA

    bond = BondPerpetual(issue_dt, coupon, freq_type, dc_type)
    settle_dt = Date(15, 7, 2026)  # on coupon date
    y = 0.05
    f = 2.0

    dv01 = bond.dv01(settle_dt, y)
    assert dv01 > 0.0

    dollar_dur = bond.dollar_duration(settle_dt, y)
    assert dollar_dur == pytest.approx(10000.0 * dv01, rel=1e-12)

    mod_dur = bond.modified_duration(settle_dt, y)
    dp = bond.dirty_price_from_ytm(settle_dt, y)
    assert mod_dur == pytest.approx(dollar_dur / dp, rel=1e-12)

    # On coupon date, analytical modified duration is approximately 1/y
    assert mod_dur == pytest.approx(1.0 / y, rel=0.01)

    mac_dur = bond.macaulay_duration(settle_dt, y)
    assert mac_dur == pytest.approx(mod_dur * (1.0 + y / f), rel=1e-12)
    assert mac_dur == pytest.approx((1.0 + y / f) / y, rel=0.01)

    conv = bond.convexity_from_ytm(settle_dt, y)
    assert conv > 0.0
    # Analytical convexity order of magnitude: 2 / (y^2 * (1 + y/f))
    analytical_conv = 2.0 / (y * y * (1.0 + y / f))
    assert conv == pytest.approx(analytical_conv, rel=0.05)


def test_perpetual_bond_current_yield():
    issue_dt = Date(1, 1, 2020)
    coupon = 0.04
    bond = BondPerpetual(issue_dt, coupon, FrequencyTypes.ANNUAL)
    settle_dt = Date(1, 1, 2026)

    assert bond.current_yield(settle_dt, 100.0) == pytest.approx(0.04, rel=1e-10)
    assert bond.current_yield(settle_dt, 80.0) == pytest.approx(0.05, rel=1e-10)
    assert bond.current_yield(settle_dt, 125.0) == pytest.approx(0.032, rel=1e-10)


def test_perpetual_bond_discount_curve():
    settle_dt = Date(1, 1, 2026)
    issue_dt = Date(1, 1, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = BondPerpetual(issue_dt, coupon, freq_type)

    # Semi-annually compounded curve matches ~100 par for 5% coupon
    curve_semi = FlatDiscountCurve(settle_dt, 0.05, freq_type=FrequencyTypes.SEMI_ANNUAL)
    dp_semi = bond.dirty_price_from_discount_curve(settle_dt, curve_semi)
    cp_semi = bond.clean_price_from_discount_curve(settle_dt, curve_semi)
    val_semi = bond.value(settle_dt, curve_semi)

    assert dp_semi == pytest.approx(100.0, rel=0.005)
    assert cp_semi == pytest.approx(100.0, rel=0.005)
    assert val_semi == pytest.approx(cp_semi, rel=1e-12)

    # Continuously compounded curve: r_continuous=0.05 corresponds to r_semi ~ 5.063%
    curve_cont = FlatDiscountCurve(settle_dt, 0.05)
    dp_cont = bond.dirty_price_from_discount_curve(settle_dt, curve_cont)
    assert dp_cont == pytest.approx(98.70, rel=0.005)


def test_perpetual_bond_vectorization():
    issue_dt = Date(15, 1, 1923)
    bond = BondPerpetual(issue_dt, 0.05, FrequencyTypes.SEMI_ANNUAL)
    settle_dt = Date(15, 7, 2026)

    ytms = np.array([0.04, 0.05, 0.06])
    dps = bond.dirty_price_from_ytm(settle_dt, ytms)
    assert len(dps) == 3
    assert dps[1] == pytest.approx(100.0, rel=1e-10)

    cps = [125.0, 100.0, 83.33333333]
    calc_ytms = bond.yield_to_maturity(settle_dt, cps)
    assert len(calc_ytms) == 3
    assert calc_ytms[0] == pytest.approx(0.04, abs=1e-6)
    assert calc_ytms[1] == pytest.approx(0.05, abs=1e-6)
    assert calc_ytms[2] == pytest.approx(0.06, abs=1e-6)


def test_perpetual_bond_input_validation():
    issue_dt = Date(15, 1, 2020)

    # Coupon must be positive
    with pytest.raises(FinError, match="Coupon must be strictly positive"):
        BondPerpetual(issue_dt, 0.0, FrequencyTypes.ANNUAL)
    with pytest.raises(FinError, match="Coupon must be strictly positive"):
        BondPerpetual(issue_dt, -0.05, FrequencyTypes.ANNUAL)

    # Zero frequency invalid
    with pytest.raises(FinError, match="Zero coupon perpetual bond is not valid"):
        BondPerpetual(issue_dt, 0.05, FrequencyTypes.ZERO)

    # Settlement before issue date
    bond = BondPerpetual(issue_dt, 0.05, FrequencyTypes.ANNUAL)
    with pytest.raises(FinError, match="before issue date"):
        bond.dirty_price_from_ytm(Date(15, 1, 2019), 0.05)

    # Non-positive YTM
    with pytest.raises(FinError, match="strictly positive"):
        bond.dirty_price_from_ytm(Date(15, 1, 2025), 0.0)
    with pytest.raises(FinError, match="strictly positive"):
        bond.dirty_price_from_ytm(Date(15, 1, 2025), -0.01)

    # Ex-div days > 90
    with pytest.raises(FinError, match="Ex-dividend days cannot exceed 90"):
        BondPerpetual(issue_dt, 0.05, FrequencyTypes.ANNUAL, ex_div_days=100)

    # First coupon before issue
    with pytest.raises(FinError, match="after issue date"):
        BondPerpetual(issue_dt, 0.05, FrequencyTypes.ANNUAL, first_coupon_dt=Date(1, 1, 2020))


def test_perpetual_bond_repr_and_print():
    issue_dt = Date(15, 1, 1923)
    bond = BondPerpetual(issue_dt, 0.025, FrequencyTypes.QUARTERLY)
    rep = repr(bond)
    assert "BondPerpetual" in rep
    assert "15-JAN-1923" in rep
    assert "0.025" in rep
    bond._print()

def _reference_dirty_price(bond, settle_dt, curve, last_year=2190):
    # Discount every coupon off the curve to the end of the supported date
    # range, then add the geometric remainder using the forward rate the curve
    # extrapolates over the final year, averaged over a year so that the day
    # counts of single periods do not bias it.
    bond._calc_pcd_ncd(settle_dt)
    freq = int(bond.freq)
    months_step = 12 // freq
    pv = 0.0
    dfs = []
    i = 0
    while True:
        flow_dt = bond._ncd.add_months(i * months_step)
        if flow_dt.y > last_year:
            break
        dfs.append(curve.df(flow_dt))
        if flow_dt > settle_dt:
            pv += bond.cpn / freq * dfs[-1]
        i += 1
    q = (dfs[-1] / dfs[-1 - freq]) ** (1.0 / freq)
    pv += bond.cpn / freq * dfs[-1] * q / (1.0 - q)
    return pv / curve.df(settle_dt) * bond.par


def test_perpetual_bond_sloped_discount_curve_all_frequencies():
    # An upward-sloping curve: the perpetuity tail must continue the curve's
    # forward rates rather than the spot zero rate at the last summed coupon.
    settle_dt = Date(15, 1, 2026)
    issue_dt = Date(15, 7, 2025)
    years = [1, 2, 5, 10, 20, 30, 50]
    zero_rates = [0.02, 0.023, 0.03, 0.035, 0.04, 0.043, 0.045]
    curve = ZeroRatesDiscountCurve(
        settle_dt,
        [settle_dt.add_years(y) for y in years],
        zero_rates,
        FrequencyTypes.CONTINUOUS,
        DayCountTypes.ACT_365F,
    )

    for freq_type in [
        FrequencyTypes.ANNUAL,
        FrequencyTypes.SEMI_ANNUAL,
        FrequencyTypes.QUARTERLY,
        FrequencyTypes.MONTHLY,
    ]:
        bond = BondPerpetual(issue_dt, 0.05, freq_type)
        dp = bond.dirty_price_from_discount_curve(settle_dt, curve)
        expected = _reference_dirty_price(bond, settle_dt, curve)
        assert dp == pytest.approx(expected, abs=5e-3), freq_type

        # The horizon only moves coupons between the sum and the tail
        dp_60 = bond.dirty_price_from_discount_curve(settle_dt, curve, horizon_years=60.0)
        assert dp_60 == pytest.approx(dp, abs=5e-3), freq_type


def test_perpetual_bond_flat_curve_matches_ytm_formula():
    settle_dt = Date(15, 1, 2026)
    issue_dt = Date(15, 1, 2020)
    for freq_type in [FrequencyTypes.ANNUAL, FrequencyTypes.QUARTERLY, FrequencyTypes.MONTHLY]:
        bond = BondPerpetual(issue_dt, 0.05, freq_type)
        curve = FlatDiscountCurve(settle_dt, 0.04, freq_type=freq_type)
        dp_curve = bond.dirty_price_from_discount_curve(settle_dt, curve)
        # On a coupon date the dirty price is the coupon over the yield
        assert dp_curve == pytest.approx(100.0 * 0.05 / 0.04, rel=2e-3), freq_type
