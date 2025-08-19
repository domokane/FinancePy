##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.hw_tree import HWTree, FinHWEuropeanCalcType
from financepy.utils.date import Date
import numpy as np


def test_HullWhiteExampleTwo():
    # HULL BOOK ZERO COUPON BOND EXAMPLE 28.1 SEE TABLE 28.3
    # Replication may not be exact as I am using dates rather than times

    zeroDays = [
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

    times = np.array(zeroDays) / 365.0
    zeros = np.array(zero_rates) / 100.0
    dfs = np.exp(-zeros * times)

    start_dt = Date(1, 12, 2019)
    sigma = 0.01
    a = 0.1
    strike = 63.0
    face = 100.0

    expiry_dt = start_dt.add_tenor("3Y")
    maturity_dt = start_dt.add_tenor("9Y")

    t_exp = (expiry_dt - start_dt) / G_DAYS_IN_YEARS
    t_mat = (maturity_dt - start_dt) / G_DAYS_IN_YEARS

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

    v_tree_call = (v_tree1["call"] + v_tree2["call"]) / 2.0
    v_tree_put = (v_tree1["put"] + v_tree2["put"]) / 2.0

    assert round(v_tree_call, 4) == 1.0450
    assert round(v_anal["call"], 4) == 1.0448
    assert round(v_tree_put, 4) == 1.8237
    assert round(v_anal["put"], 4) == 1.8239
