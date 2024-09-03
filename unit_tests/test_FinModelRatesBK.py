##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.global_types import FinExerciseTypes
from financepy.models.bk_tree import BKTree
from financepy.utils.global_vars import g_days_in_year
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date
import numpy as np


def test_BKExampleOne():
    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    start_dt = Date(1, 12, 2019)
    end_dt = Date(1, 6, 2021)
    sigma = 0.25
    a = 0.22
    num_time_steps = 3
    t_mat = (end_dt - start_dt)/g_days_in_year
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)

    assert [round(x, 4) for x in model.Q[2]] == \
        [0.0190, 0.2126, 0.5009, 0.2112, 0.0187]

    assert [round(x, 4) for x in model.rt[2]] == \
        [0.0259, 0.0351, 0.0477, 0.0648, 0.0881]

    assert [round(x, 4) for x in model.pu] == \
        [0.0808, 0.2278, 0.1667, 0.1177, 0.8606]

    assert [round(x, 4) for x in model.pm] == \
        [0.0586, 0.6545, 0.6667, 0.6545, 0.0586]

    assert [round(x, 4) for x in model.pd] == \
        [0.8606, 0.1177, 0.1667, 0.2278, 0.0808]


def test_BKExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2018)
    expiry_dt = settle_dt.add_tenor("18m")
    maturity_dt = settle_dt.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []

    cpn = bond.cpn / bond.freq

    num_flows = len(bond.cpn_dts)

    for i in range(1, num_flows):
        pcd = bond.cpn_dts[i-1]
        ncd = bond.cpn_dts[i]
        if pcd < settle_dt and ncd > settle_dt:
            flow_time = (pcd - settle_dt) / g_days_in_year
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    for flow_dt in bond.cpn_dts:
        if flow_dt > settle_dt:
            flow_time = (flow_dt - settle_dt) / g_days_in_year
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 105.0
    face = 100.0

    t_mat = (maturity_dt - settle_dt) / g_days_in_year
    t_exp = (expiry_dt - settle_dt) / g_days_in_year
    times = np.linspace(0, t_mat, 11)
    dates = settle_dt.add_years(times)
    dfs = np.exp(-0.05*times)
    curve = DiscountCurve(settle_dt, dates, dfs)

    price = bond.clean_price_from_discount_curve(settle_dt, curve)
    assert round(price, 4) == 99.5420

    sigma = 0.20
    a = 0.05
    num_time_steps = 26

    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    exercise_type = FinExerciseTypes.AMERICAN
    v = model.bond_option(t_exp, strike_price, face, cpn_times,
                          cpn_flows, exercise_type)

    # Test convergence
    num_time_steps = 200
    exercise_type = FinExerciseTypes.AMERICAN

    treeVector = []
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(t_mat, times, dfs)
    v = model.bond_option(t_exp, strike_price,
                          face, cpn_times, cpn_flows, exercise_type)
    treeVector.append(v)

    assert round(v['call'], 4) == 0.6998
    assert round(v['put'], 4) == 7.9605

test_BKExampleOne()
