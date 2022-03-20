##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.global_types import FinExerciseTypes
from financepy.models.bk_tree import BKTree
from financepy.utils.global_vars import gDaysInYear
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

    start_date = Date(1, 12, 2019)
    end_date = Date(1, 6, 2021)
    sigma = 0.25
    a = 0.22
    num_time_steps = 3
    tmat = (end_date - start_date)/gDaysInYear
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(tmat, times, dfs)

    assert [round(x, 4) for x in model._Q[2]] == \
        [0.0190, 0.2126, 0.5009, 0.2112, 0.0187]

    assert [round(x, 4) for x in model._rt[2]] == \
        [0.0259, 0.0351, 0.0477, 0.0648, 0.0881]

    assert [round(x, 4) for x in model._pu] == \
        [0.0808, 0.2278, 0.1667, 0.1177, 0.8606]

    assert [round(x, 4) for x in model._pm] == \
        [0.0586, 0.6545, 0.6667, 0.6545, 0.0586]

    assert [round(x, 4) for x in model._pd] == \
        [0.8606, 0.1177, 0.1667, 0.2278, 0.0808]


def test_BKExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    expiry_date = settlement_date.add_tenor("18m")
    maturity_date = settlement_date.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    coupon_times = []
    coupon_flows = []
    cpn = bond._coupon/bond._frequency
    num_flows = len(bond._coupon_dates)

    for i in range(1, num_flows):
        pcd = bond._coupon_dates[i-1]
        ncd = bond._coupon_dates[i]
        if pcd < settlement_date and ncd > settlement_date:
            flow_time = (pcd - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    for flow_date in bond._coupon_dates:
        if flow_date > settlement_date:
            flow_time = (flow_date - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    strike_price = 105.0
    face = 100.0

    tmat = (maturity_date - settlement_date) / gDaysInYear
    texp = (expiry_date - settlement_date) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settlement_date.add_years(times)
    dfs = np.exp(-0.05*times)
    curve = DiscountCurve(settlement_date, dates, dfs)

    price = bond.clean_price_from_discount_curve(settlement_date, curve)
    assert round(price, 4) == 99.5420

    sigma = 0.20
    a = 0.05
    num_time_steps = 26

    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(tmat, times, dfs)
    exercise_type = FinExerciseTypes.AMERICAN
    v = model.bond_option(texp, strike_price, face, coupon_times,
                          coupon_flows, exercise_type)

    # Test convergence
    num_time_steps = 200
    exercise_type = FinExerciseTypes.AMERICAN

    treeVector = []
    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(tmat, times, dfs)
    v = model.bond_option(texp, strike_price,
                          face, coupon_times, coupon_flows, exercise_type)
    treeVector.append(v)

    assert round(v['call'], 4) == 0.6998
    assert round(v['put'], 4) == 7.9605
