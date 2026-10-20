# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.models.hw_tree import HWTree
from financepy.utils.date import Date

########################################################################################


def test_hull_white_example_two():

    # HULL BOOK ZERO COUPON BOND EXAMPLE 28.1 SEE TABLE 28.3
    # Replication may not be exact as I am using dates rather than times

    zero_days = [
        0,
        3,
        31,
        62,
        94,
        185,
        367,
        731,
        1096,
        1461,
        1826,
        2194,
        2558,
        2922,
        3287,
        3653,
    ]

    zero_rates = [
        5.0,
        5.01772,
        4.98282,
        4.97234,
        4.96157,
        4.99058,
        5.09389,
        5.79733,
        6.30595,
        6.73464,
        6.94816,
        7.08807,
        7.27527,
        7.30852,
        7.39790,
        7.49015,
    ]

    times = np.array(zero_days) / 365.0
    zeros = np.array(zero_rates) / 100.0
    dfs = np.exp(-zeros * times)

    start_dt = Date(1, 12, 2019)
    sigma = 0.01
    a = 0.1
    strike = 63.0
    face = 100.0

    expiry_dt = start_dt.add_tenor("3Y")
    maturity_dt = start_dt.add_tenor("9Y")

    t_exp = (expiry_dt - start_dt) / G_DAYS_IN_YEAR
    t_mat = (maturity_dt - start_dt) / G_DAYS_IN_YEAR

    num_time_steps = None
    model = HWTree(sigma, a, num_time_steps)
    v_anal = model.option_on_zcb(t_exp, t_mat, strike, face, times, dfs)

    num_time_steps = 200

    model = HWTree(sigma, a, num_time_steps)
    model.build_tree(t_exp, times, dfs)
    v_tree1 = model.option_on_zero_cpn_bond_tree(t_exp, t_mat, strike, face)

    model = HWTree(sigma, a, num_time_steps + 1)
    model.build_tree(t_exp, times, dfs)
    v_tree2 = model.option_on_zero_cpn_bond_tree(t_exp, t_mat, strike, face)

    v_tree_call = (v_tree1[0] + v_tree2[0]) / 2.0
    v_tree_put = (v_tree1[1] + v_tree2[1]) / 2.0

    assert round(v_tree_call, 4) == 1.0450
    assert round(v_anal[0], 4) == 1.0448
    assert round(v_tree_put, 4) == 1.8237
    assert round(v_anal[1], 4) == 1.8239


########################################################################################


def test_hw_df_tree_matches_input_curve():
    """The sum of the Arrow-Debreu prices at the tree step for time t is the
    discount factor to t, not to t + dt."""

    r = 0.04
    times = np.linspace(0.0, 10.0, 41)
    dfs = np.exp(-r * times)

    hw = HWTree(0.01, 0.1, 200)
    hw.build_tree(5.0, times, dfs)

    for k in (1, 10, 50, 100):
        t = k * hw.dt
        df, zero_rate = hw.df_tree(t)
        assert abs(df - np.exp(-r * t)) < 1e-10
        assert abs(zero_rate - r) < 1e-10


def test_hw_expiry_tree_bond_option_flow_just_after_expiry():
    """A coupon that is adjusted to fall between the option expiry and the
    nearest tree time must not be counted in the dirty price without being
    counted in the accrued. The three European calculation types should agree."""

    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.day_count import DayCountTypes
    from financepy.utils.global_types import OptionTypes
    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
    from financepy.models.hw_tree import HWEuropeanCalcTypes
    from financepy.products.bonds.bond import Bond
    from financepy.products.bonds.bond_option import BondOption

    value_dt = Date(20, 10, 2026)
    curve = FlatDiscountCurve(
        value_dt, 0.04, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    bond = Bond(
        Date(20, 10, 2024),
        Date(20, 10, 2034),
        0.05,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
    )

    # Saturday coupon date: the payment is adjusted to 22 October, two days
    # after the option expiry.
    expiry_dt = Date(20, 10, 2029)

    for option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        option = BondOption(bond, expiry_dt, 100.0, option_type)
        values = {}
        for calc_type in HWEuropeanCalcTypes:
            model = HWTree(0.01, 0.1, 500, calc_type)
            values[calc_type] = option.value(value_dt, curve, model)

        jamshidian = values[HWEuropeanCalcTypes.JAMSHIDIAN]
        expiry_tree = values[HWEuropeanCalcTypes.EXPIRY_TREE]
        expiry_only = values[HWEuropeanCalcTypes.EXPIRY_ONLY]

        assert abs(expiry_tree - expiry_only) < 1e-3
        assert abs(expiry_tree - jamshidian) / jamshidian < 5e-3


def test_hw_jamshidian_swaption_put_call_parity():
    """Payer minus receiver equals the forward swap value when the root search
    and the final pricing use the same discount factor bump."""

    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.day_count import DayCountTypes
    from financepy.utils.global_types import SwapTypes
    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
    from financepy.products.rates.ibor_swaption import IborSwaption

    value_dt = Date(20, 10, 2026)
    curve = FlatDiscountCurve(
        value_dt, 0.04, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    exercise_dt = Date(20, 10, 2028)
    maturity_dt = Date(20, 10, 2033)
    notional = 1_000_000.0

    for strike in (0.03, 0.04, 0.05):
        payer = IborSwaption(
            value_dt, exercise_dt, maturity_dt, SwapTypes.PAY, strike,
            FrequencyTypes.ANNUAL, DayCountTypes.THIRTY_E_360, notional,
        )
        receiver = IborSwaption(
            value_dt, exercise_dt, maturity_dt, SwapTypes.RECEIVE, strike,
            FrequencyTypes.ANNUAL, DayCountTypes.THIRTY_E_360, notional,
        )
        v_pay = payer.value(value_dt, curve, HWTree(0.01, 0.1, 100))
        v_rec = receiver.value(value_dt, curve, HWTree(0.01, 0.1, 100))
        forward = payer.pv01 * (payer.fwd_swap_rate - strike) * notional
        assert abs((v_pay - v_rec) - forward) < 1e-6 * notional


def test_hw_expiry_only_bond_option_on_coupon_date():
    """When the option expires exactly on a coupon date the dirty price at
    expiry includes that coupon, so the accrued must count it as fully
    accrued. The expiry-only calculation used to leave it out and overstate
    the clean price by a full coupon."""

    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.day_count import DayCountTypes
    from financepy.utils.global_types import OptionTypes
    from financepy.market.curves.discount_curve import DiscountCurve
    from financepy.models.hw_tree import HWEuropeanCalcTypes
    from financepy.products.bonds.bond import Bond
    from financepy.products.bonds.bond_option import BondOption

    settle_dt = Date(1, 12, 2019)
    bond = Bond(
        Date(1, 12, 2018),
        settle_dt.add_tenor("10Y"),
        0.05,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
    )
    t_mat = (bond.maturity_dt - settle_dt) / G_DAYS_IN_YEAR
    times = np.linspace(0, t_mat, 20)
    curve = DiscountCurve(settle_dt, settle_dt.add_years(times), np.exp(-0.05 * times))

    expiry_dt = settle_dt.add_tenor("18m")
    assert expiry_dt in bond.cpn_dts

    for option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        option = BondOption(bond, expiry_dt, 100.0, option_type)
        values = {
            calc_type: option.value(settle_dt, curve, HWTree(0.01, 0.1, 100, calc_type))
            for calc_type in HWEuropeanCalcTypes
        }
        jamshidian = values[HWEuropeanCalcTypes.JAMSHIDIAN]
        for calc_type in (HWEuropeanCalcTypes.EXPIRY_ONLY, HWEuropeanCalcTypes.EXPIRY_TREE):
            assert abs(values[calc_type] - jamshidian) / jamshidian < 0.02
