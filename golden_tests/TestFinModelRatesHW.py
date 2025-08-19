##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys

sys.path.append("..")

import time
import numpy as np
from financepy.utils.date import Date
from financepy.models.hw_tree import HWTree, FinHWEuropeanCalcType
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.utils.helpers import print_tree
from financepy.utils.global_types import FinExerciseTypes
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_HullWhiteExampleOne():
    # HULL BOOK INITIAL EXAMPLE SECTION 28.7 HW EDITION 6

    times = [0.0, 0.5000, 1.00000, 1.50000, 2.00000, 2.500000, 3.00000]
    zeros = [0.03, 0.0343, 0.03824, 0.04183, 0.04512, 0.048512, 0.05086]
    times = np.array(times)
    zeros = np.array(zeros)
    dfs = np.exp(-zeros * times)

    start_dt = Date(1, 12, 2019)
    end_dt = Date(1, 12, 2022)
    sigma = 0.01
    a = 0.1
    num_time_steps = 3
    model = HWTree(sigma, a, num_time_steps)
    tree_mat = (end_dt - start_dt) / G_DAYS_IN_YEARS
    model.build_tree(tree_mat, times, dfs)


#   print_tree(model._Q)
#   print("")
#   print_tree(model._rt)
#   print("")

###############################################################################


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

    # Test convergence
    num_steps_list = range(100, 500, 100)
    anal_vector = []
    tree_vector = []

    test_cases.banner("Comparing option on zero coupon bond analytical vs Tree")

    test_cases.header(
        "NUMTIMESTEP",
        "TIME",
        "VTREE_CALL",
        "VTREE_PUT",
        "VANAL CALL",
        "VANAL_PUT",
        "CALLDIFF",
        "PUTDIFF",
    )

    for num_time_steps in num_steps_list:

        start = time.time()

        model = HWTree(sigma, a, num_time_steps)
        model.build_tree(t_exp, times, dfs)
        v_tree1 = model.option_on_zero_cpn_bond_tree(t_exp, t_mat, strike, face)

        model = HWTree(sigma, a, num_time_steps + 1)
        model.build_tree(t_exp, times, dfs)
        v_tree2 = model.option_on_zero_cpn_bond_tree(t_exp, t_mat, strike, face)

        end = time.time()
        period = end - start
        tree_vector.append(v_tree1["put"])
        anal_vector.append(v_anal["put"])
        v_tree_call = (v_tree1["call"] + v_tree2["call"]) / 2.0
        v_tree_put = (v_tree1["put"] + v_tree2["put"]) / 2.0
        diff_c = v_tree_call - v_anal["call"]
        diff_p = v_tree_put - v_anal["put"]

        test_cases.print(
            num_time_steps,
            period,
            v_tree_call,
            v_anal["call"],
            v_tree_put,
            v_anal["put"],
            diff_c,
            diff_p,
        )


#   plt.plot(num_steps_list, tree_vector)
#   plt.plot(num_steps_list, anal_vector)

###############################################################################


def test_HullWhiteBondOption():
    # Valuation of a European option on a cpn bearing bond

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2018)
    expiry_dt = settle_dt.add_tenor("18m")
    maturity_dt = settle_dt.add_tenor("10Y")
    cpn = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, cpn, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []
    cpn = bond.cpn / bond.freq

    num_flows = len(bond.cpn_dts)

    for i in range(1, num_flows):

        pcd = bond.cpn_dts[i - 1]
        ncd = bond.cpn_dts[i]

        if ncd > settle_dt:

            if len(cpn_times) == 0:
                flow_time = (pcd - settle_dt) / G_DAYS_IN_YEARS
                cpn_times.append(flow_time)
                cpn_flows.append(cpn)

            flow_time = (ncd - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    strike_price = 100.0
    face = 100.0
    y = 0.05
    times = np.linspace(0, 10, 21)
    dfs = np.power(1 + y / 2, -times * 2)

    sigma = 0.0000001
    a = 0.1
    model = HWTree(sigma, a, None)

    #  Test convergence
    num_steps_list = range(50, 500, 50)
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEARS

    v_jam = model.european_bond_option_jamshidian(
        t_exp, strike_price, face, cpn_times, cpn_flows, times, dfs
    )

    test_cases.banner(
        "Pricing bond option on tree that goes to bond maturity and one using european bond option tree that goes to expiry."
    )

    test_cases.header("NUMSTEPS", "TIME", "EXPIRY_ONLY", "EXPIRY_TREE", "JAMSHIDIAN")

    for num_time_steps in num_steps_list:

        start = time.time()
        model = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        model.build_tree(t_exp, times, dfs)

        exercise_type = FinExerciseTypes.EUROPEAN

        v1 = model.bond_option(
            t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
        )

        model = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_TREE)
        model.build_tree(t_exp, times, dfs)

        v2 = model.bond_option(
            t_exp, strike_price, face, cpn_times, cpn_flows, exercise_type
        )

        end = time.time()
        period = end - start

        test_cases.print(num_time_steps, period, v1, v2, v_jam)

    #    plt.plot(num_steps_list, tree_vector)

    if 1 == 0:
        print("RT")
        print_tree(model._rt, 5)
        print("BOND")
        print_tree(model._bond_values, 5)
        print("OPTION")
        print_tree(model._option_values, 5)


###############################################################################


def test_HullWhiteCallableBond():
    # Valuation of a European option on a cpn bearing bond

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2018)
    maturity_dt = settle_dt.add_tenor("10Y")
    cpn = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, cpn, freq_type, dc_type)

    cpn_times = []
    cpn_flows = []
    cpn = bond.cpn / bond.freq

    for flow_dt in bond.cpn_dts[1:]:

        if flow_dt > settle_dt:
            flow_time = (flow_dt - settle_dt) / G_DAYS_IN_YEARS
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    call_dts = []
    call_prices = []
    call_price = 120.0
    call_dts.append(settle_dt.add_tenor("2Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("3Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("4Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("5Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("6Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("7Y"))
    call_prices.append(call_price)
    call_dts.append(settle_dt.add_tenor("8Y"))
    call_prices.append(call_price)

    call_times = []
    for dt in call_dts:
        t = (dt - settle_dt) / G_DAYS_IN_YEARS
        call_times.append(t)

    put_dts = []
    put_prices = []
    put_price = 98.0
    put_dts.append(settle_dt.add_tenor("2Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("3Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("4Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("5Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("6Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("7Y"))
    put_prices.append(put_price)
    put_dts.append(settle_dt.add_tenor("8Y"))
    put_prices.append(put_price)

    put_times = []
    for dt in put_dts:
        t = (dt - settle_dt) / G_DAYS_IN_YEARS
        put_times.append(t)

    ###########################################################################

    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS
    curve = DiscountCurveFlat(settle_dt, 0.05, FrequencyTypes.CONTINUOUS)

    dfs = []
    times = []

    for dt in bond.cpn_dts:
        if dt > settle_dt:
            t = (dt - settle_dt) / G_DAYS_IN_YEARS
            df = curve.df(dt)
            times.append(t)
            dfs.append(df)

    dfs = np.array(dfs)
    times = np.array(times)

    ###########################################################################

    v1 = bond.clean_price_from_discount_curve(settle_dt, curve)

    sigma = 0.02  # basis point volatility
    a = 0.01

    # Test convergence
    num_steps_list = [100, 200, 500, 1000]
    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS

    test_cases.header("NUMSTEPS", "TIME", "BOND_ONLY", "CALLABLE_BOND")

    for num_time_steps in num_steps_list:

        start = time.time()
        model = HWTree(sigma, a, num_time_steps)
        model.build_tree(t_mat, times, dfs)

        v2 = model.callable_puttable_bond_tree(
            cpn_times, cpn_flows, call_times, call_prices, put_times, put_prices, 100.0
        )

        end = time.time()
        period = end - start
        test_cases.print(num_time_steps, period, v1, v2)


###############################################################################


test_HullWhiteExampleOne()
test_HullWhiteExampleTwo()
test_HullWhiteBondOption()
test_HullWhiteCallableBond()
test_cases.compareTestCases()
