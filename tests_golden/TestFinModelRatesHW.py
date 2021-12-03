##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import time
import numpy as np
from financepy.utils.date import Date
from financepy.models.hw_tree import HWTree, FinHWEuropeanCalcType
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.helpers import print_tree
from financepy.utils.global_types import FinExerciseTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_HullWhiteExampleOne():
    # HULL BOOK INITIAL EXAMPLE SECTION 28.7 HW EDITION 6

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros*times)

    start_date = Date(1, 12, 2019)
    end_date = Date(1, 12, 2022)
    sigma = 0.01
    a = 0.1
    num_time_steps = 3
    model = HWTree(sigma, a, num_time_steps)
    treeMat = (end_date - start_date)/gDaysInYear
    model.build_tree(treeMat, times, dfs)
#   print_tree(model._Q)
#   print("")
#   print_tree(model._rt)
#   print("")

###############################################################################


def test_HullWhiteExampleTwo():
    # HULL BOOK ZERO COUPON BOND EXAMPLE 28.1 SEE TABLE 28.3
    # Replication may not be exact as I am using dates rather than times

    zeroDays = [0, 3, 31, 62, 94, 185, 367, 731, 1096, 1461,
                1826, 2194, 2558, 2922, 3287, 3653]

    zero_rates = [5.0, 5.01772, 4.98282, 4.97234, 4.96157, 4.99058, 5.09389,
                  5.79733, 6.30595, 6.73464, 6.94816, 7.08807, 7.27527,
                  7.30852, 7.39790, 7.49015]

    times = np.array(zeroDays) / 365.0
    zeros = np.array(zero_rates) / 100.0
    dfs = np.exp(-zeros*times)

    start_date = Date(1, 12, 2019)
    sigma = 0.01
    a = 0.1
    strike = 63.0
    face = 100.0

    expiry_date = start_date.add_tenor("3Y")
    maturity_date = start_date.add_tenor("9Y")

    texp = (expiry_date - start_date)/gDaysInYear
    tmat = (maturity_date - start_date)/gDaysInYear

    num_time_steps = None
    model = HWTree(sigma, a, num_time_steps)
    vAnal = model.option_on_zcb(texp, tmat, strike, face, times, dfs)

    # Test convergence
    num_steps_list = range(100, 500, 100)
    analVector = []
    treeVector = []

    testCases.banner("Comparing option on zero coupon bond analytical vs Tree")

    testCases.header("NUMTIMESTEP", "TIME", "VTREE_CALL", "VTREE_PUT",
                     "VANAL CALL", "VANAL_PUT", "CALLDIFF", "PUTDIFF")

    for num_time_steps in num_steps_list:

        start = time.time()

        model = HWTree(sigma, a, num_time_steps)
        model.build_tree(texp, times, dfs)
        vTree1 = model.option_on_zero_coupon_bond_tree(
            texp, tmat, strike, face)

        model = HWTree(sigma, a, num_time_steps+1)
        model.build_tree(texp, times, dfs)
        vTree2 = model.option_on_zero_coupon_bond_tree(
            texp, tmat, strike, face)

        end = time.time()
        period = end-start
        treeVector.append(vTree1['put'])
        analVector.append(vAnal['put'])
        vTreeCall = (vTree1['call'] + vTree2['call']) / 2.0
        vTreePut = (vTree1['put'] + vTree2['put']) / 2.0
        diffC = vTreeCall - vAnal['call']
        diffP = vTreePut - vAnal['put']

        testCases.print(num_time_steps, period, vTreeCall, vAnal['call'],
                        vTreePut, vAnal['put'], diffC, diffP)

 #   plt.plot(num_steps_list, treeVector)
 #   plt.plot(num_steps_list, analVector)

###############################################################################


def test_HullWhiteBondOption():
    # Valuation of a European option on a coupon bearing bond

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

    num_flows = len(bond._flow_dates)
    for i in range(1, num_flows):

        pcd = bond._flow_dates[i-1]
        ncd = bond._flow_dates[i]

        if ncd > settlement_date:

            if len(coupon_times) == 0:
                flow_time = (pcd - settlement_date) / gDaysInYear
                coupon_times.append(flow_time)
                coupon_flows.append(cpn)

            flow_time = (ncd - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    strike_price = 100.0
    face = 100.0
    y = 0.05
    times = np.linspace(0, 10, 21)
    dfs = np.power(1+y/2, -times*2)

    sigma = 0.0000001
    a = 0.1
    model = HWTree(sigma, a, None)

    #  Test convergence
    num_steps_list = range(50, 500, 50)
    texp = (expiry_date - settlement_date)/gDaysInYear

    vJam = model.european_bond_option_jamshidian(texp, strike_price, face,
                                                 coupon_times, coupon_flows,
                                                 times, dfs)

    testCases.banner(
        "Pricing bond option on tree that goes to bond maturity and one using european bond option tree that goes to expiry.")

    testCases.header("NUMSTEPS", "TIME", "EXPIRY_ONLY",
                     "EXPIRY_TREE", "JAMSHIDIAN")

    for num_time_steps in num_steps_list:

        start = time.time()
        model = HWTree(sigma, a, num_time_steps,
                       FinHWEuropeanCalcType.EXPIRY_ONLY)
        model.build_tree(texp, times, dfs)

        exercise_type = FinExerciseTypes.EUROPEAN

        v1 = model.bond_option(texp, strike_price, face,
                               coupon_times, coupon_flows, exercise_type)

        model = HWTree(sigma, a, num_time_steps,
                       FinHWEuropeanCalcType.EXPIRY_TREE)
        model.build_tree(texp, times, dfs)

        v2 = model.bond_option(texp, strike_price, face,
                               coupon_times, coupon_flows, exercise_type)

        end = time.time()
        period = end-start

        testCases.print(num_time_steps, period, v1, v2, vJam)

#    plt.plot(num_steps_list, treeVector)

    if 1 == 0:
        print("RT")
        print_tree(model._rt, 5)
        print("BOND")
        print_tree(model._bond_values, 5)
        print("OPTION")
        print_tree(model._option_values, 5)


###############################################################################


def test_HullWhiteCallableBond():
    # Valuation of a European option on a coupon bearing bond

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    maturity_date = settlement_date.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    coupon_times = []
    coupon_flows = []
    cpn = bond._coupon/bond._frequency

    for flow_date in bond._flow_dates[1:]:

        if flow_date > settlement_date:
            flow_time = (flow_date - settlement_date) / gDaysInYear
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    call_dates = []
    call_prices = []
    callPx = 120.0
    call_dates.append(settlement_date.add_tenor("2Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("3Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("4Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("5Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("6Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("7Y"))
    call_prices.append(callPx)
    call_dates.append(settlement_date.add_tenor("8Y"))
    call_prices.append(callPx)

    call_times = []
    for dt in call_dates:
        t = (dt - settlement_date) / gDaysInYear
        call_times.append(t)

    put_dates = []
    put_prices = []
    putPx = 98.0
    put_dates.append(settlement_date.add_tenor("2Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("3Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("4Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("5Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("6Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("7Y"))
    put_prices.append(putPx)
    put_dates.append(settlement_date.add_tenor("8Y"))
    put_prices.append(putPx)

    put_times = []
    for dt in put_dates:
        t = (dt - settlement_date) / gDaysInYear
        put_times.append(t)

    ###########################################################################

    tmat = (maturity_date - settlement_date) / gDaysInYear
    curve = DiscountCurveFlat(settlement_date, 0.05, FrequencyTypes.CONTINUOUS)

    dfs = []
    times = []

    for dt in bond._flow_dates:
        if dt > settlement_date:
            t = (dt - settlement_date) / gDaysInYear
            df = curve.df(dt)
            times.append(t)
            dfs.append(df)

    dfs = np.array(dfs)
    times = np.array(times)

    ###########################################################################

    v1 = bond.clean_price_from_discount_curve(settlement_date, curve)

    sigma = 0.02  # basis point volatility
    a = 0.01

    # Test convergence
    num_steps_list = [100, 200, 500, 1000]
    tmat = (maturity_date - settlement_date)/gDaysInYear

    testCases.header("NUMSTEPS", "TIME", "BOND_ONLY", "CALLABLE_BOND")

    for num_time_steps in num_steps_list:

        start = time.time()
        model = HWTree(sigma, a, num_time_steps)
        model.build_tree(tmat, times, dfs)

        v2 = model.callable_puttable_bond_tree(coupon_times, coupon_flows,
                                               call_times, call_prices,
                                               put_times, put_prices, 100.0)

        end = time.time()
        period = end-start
        testCases.print(num_time_steps, period, v1, v2)

###############################################################################


test_HullWhiteExampleOne()
test_HullWhiteExampleTwo()
test_HullWhiteBondOption()
test_HullWhiteCallableBond()
testCases.compareTestCases()
