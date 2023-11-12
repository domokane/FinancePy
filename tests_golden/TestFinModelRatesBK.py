##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
sys.path.append("..")

import time
import numpy as np
from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.helpers import print_tree
from financepy.models.bk_tree import BKTree
from financepy.utils.global_types import FinExerciseTypes
from FinTestCases import FinTestCases, globalTestCaseMode


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_BKExampleOne():

    testCases.banner("=== HULL INITIAL EXAMPLE SECTION 28.7 ED 6 PG 668 ====")

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

    # Agrees with Figure 28.10 - Not exact as we have dt not exactly 0.50
    if num_time_steps < 5:
        testCases.header("LABEL", "VALUE")
        testCases.print("QTREE", model._Q)
        testCases.print("RTREE", model._rt)
#        print_tree(model._rt)
        testCases.print("PU AT LAST TIME", model._pu)
        testCases.print("PDM AT LAST TIME", model._pm)
        testCases.print("PD AT LAST TIME", model._pd)

###############################################################################


def test_BKExampleTwo():
    # Valuation of a European option on a coupon bearing bond
    # This follows example in Fig 28.11 of John Hull's book but does not
    # have the exact same dt so there are some differences

    settle_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    expiry_date = settle_date.add_tenor("18m")
    maturity_date = settle_date.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []
    cpn = bond._cpn / bond._freq
    num_flows = len(bond._cpn_dates)

    for i in range(1, num_flows):
        pcd = bond._cpn_dates[i-1]
        ncd = bond._cpn_dates[i]
        if pcd < settle_date and ncd > settle_date:
            flow_time = (pcd - settle_date) / gDaysInYear
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    for flow_date in bond._cpn_dates:
        if flow_date > settle_date:
            flow_time = (flow_date - settle_date) / gDaysInYear
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 105.0
    face = 100.0

    tmat = (maturity_date - settle_date) / gDaysInYear
    t_exp = (expiry_date - settle_date) / gDaysInYear
    times = np.linspace(0, tmat, 11)
    dates = settle_date.add_years(times)
    dfs = np.exp(-0.05*times)
    curve = DiscountCurve(settle_date, dates, dfs)

    price = bond.clean_price_from_discount_curve(settle_date, curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    sigma = 0.20
    a = 0.05
    num_time_steps = 26

    model = BKTree(sigma, a, num_time_steps)
    model.build_tree(tmat, times, dfs)
    exercise_type = FinExerciseTypes.AMERICAN
    v = model.bond_option(t_exp, strike_price, face, cpn_times,
                          cpn_flows, exercise_type)

    # Test convergence
    num_steps_list = [100, 200, 300, 500, 1000]
    exercise_type = FinExerciseTypes.AMERICAN

    testCases.header("TIMESTEPS", "TIME", "VALUE")
    treeVector = []
    for num_time_steps in num_steps_list:
        start = time.time()
        model = BKTree(sigma, a, num_time_steps)
        model.build_tree(tmat, times, dfs)
        v = model.bond_option(t_exp, strike_price,
                              face, cpn_times, cpn_flows, exercise_type)
        end = time.time()
        period = end-start
        treeVector.append(v)
        testCases.print(num_time_steps, period, v)

#    plt.plot(num_steps_list, treeVector)

    # Value in Hill converges to 0.699 with 100 time steps while I get 0.700

    if 1 == 0:
        print("RT")
        print_tree(model._rt, 5)
        print("Q")
        print_tree(model._Q, 5)

###############################################################################


test_BKExampleOne()
test_BKExampleTwo()
testCases.compareTestCases()
