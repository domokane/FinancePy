# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time

import matplotlib.pyplot as plt

import numpy as np

import add_fp_to_path

from financepy.utils.global_types import ExerciseTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_types import OptionTypes
from financepy.models.hw_tree import HWTree, FinHWEuropeanCalcType
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

########################################################################################


def test_bond_option():

    settle_dt = Date(1, 12, 2019)
    issue_dt = Date(1, 12, 2018)
    maturity_dt = settle_dt.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    times = np.linspace(0, 10.0, 21)
    dfs = np.exp(-0.05 * times)
    dates = settle_dt.add_years(times)
    discount_curve = DiscountCurve(settle_dt, dates, dfs)

    expiry_dt = settle_dt.add_tenor("18m")
    strike_price = 105.0

    strikes = [80, 90, 100, 110, 120]

    opt_type = OptionTypes.EUROPEAN_CALL

    test_cases.header("LABEL", "VALUE")

    price = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
    test_cases.print("Fixed Income Price:", price)

    num_time_steps = 50

    test_cases.banner("HW EUROPEAN CALL")
    test_cases.header("STRIKE", "VALUE")

    for strike_price in strikes:

        sigma = 0.01
        a = 0.1

        bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)

        model = HWTree(sigma, a, num_time_steps)
        v = bond_option.value(settle_dt, discount_curve, model)
        test_cases.print(strike_price, v)

    opt_type = OptionTypes.AMERICAN_CALL

    price = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
    test_cases.header("LABEL", "VALUE")
    test_cases.print("Fixed Income Price:", price)

    test_cases.banner("HW AMERICAN CALL")
    test_cases.header("STRIKE", "VALUE")

    for strike_price in strikes:

        sigma = 0.01
        a = 0.1

        bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)

        model = HWTree(sigma, a)
        v = bond_option.value(settle_dt, discount_curve, model)
        test_cases.print(strike_price, v)

    opt_type = OptionTypes.EUROPEAN_PUT
    test_cases.banner("HW EUROPEAN PUT")
    test_cases.header("STRIKE", "VALUE")

    price = bond.clean_price_from_discount_curve(settle_dt, discount_curve)

    for strike_price in strikes:

        sigma = 0.01
        a = 0.1

        bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)

        model = HWTree(sigma, a)
        v = bond_option.value(settle_dt, discount_curve, model)
        test_cases.print(strike_price, v)

    opt_type = OptionTypes.AMERICAN_PUT
    test_cases.banner("HW AMERICAN PUT")
    test_cases.header("STRIKE", "VALUE")

    price = bond.clean_price_from_discount_curve(settle_dt, discount_curve)

    for strike_price in strikes:

        sigma = 0.02
        a = 0.1

        bond_option = BondOption(bond, expiry_dt, strike_price, opt_type)

        model = HWTree(sigma, a)
        v = bond_option.value(settle_dt, discount_curve, model)
        test_cases.print(strike_price, v)


########################################################################################


def test_bond_option_european_convergence():

    # CONVERGENCE TESTS
    # COMPARE AMERICAN TREE VERSUS JAMSHIDIAN IN EUROPEAN LIMIT TO CHECK THAT
    # TREE HAS BEEN CORRECTLY CONSTRUCTED. FIND VERY GOOD AGREEMENT.

    # Build discount curve
    settle_dt = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(settle_dt, 0.05, FrequencyTypes.CONTINUOUS)

    # Bond details
    issue_dt = Date(1, 12, 2015)
    maturity_dt = Date(1, 12, 2020)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    # Option Details - put expiry in the middle of a coupon period
    expiry_dt = Date(1, 3, 2020)
    strike_price = 100.0
    face = 100.0

    time_steps = range(100, 400, 100)
    strike_price = 100.0

    test_cases.header("TIME", "N", "PUT_JAM", "PUT_TREE", "CALL_JAM", "CALL_TREE")

    for num_time_steps in time_steps:

        sigma = 0.05
        a = 0.1

        start = time.time()
        opt_type = OptionTypes.EUROPEAN_PUT

        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model1 = HWTree(sigma, a, num_time_steps)
        v1put = bond_option1.value(settle_dt, discount_curve, model1)

        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bond_option2.value(settle_dt, discount_curve, model2)

        opt_type = OptionTypes.EUROPEAN_CALL

        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model1 = HWTree(sigma, a, num_time_steps)
        v1call = bond_option1.value(settle_dt, discount_curve, model1)

        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bond_option2.value(settle_dt, discount_curve, model2)

        end = time.time()
        period = end - start
        test_cases.print(period, num_time_steps, v1put, v2put, v1call, v2call)


########################################################################################


def test_bond_option_american_convergence_one():

    # Build discount curve
    settle_dt = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(settle_dt, 0.05)

    # Bond details
    issue_dt = Date(1, 9, 2014)
    maturity_dt = Date(1, 9, 2025)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    # Option Details
    expiry_dt = Date(1, 12, 2020)
    strike_price = 100.0
    face = 100.0

    test_cases.header("TIME", "N", "PUT_AMER", "PUT_EUR", "CALL_AME", "CALL_EUR")

    time_steps = range(100, 400, 100)

    for num_time_steps in time_steps:

        sigma = 0.05
        a = 0.1

        start = time.time()

        opt_type = OptionTypes.AMERICAN_PUT
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model1 = HWTree(sigma, a, num_time_steps)
        v1put = bond_option1.value(settle_dt, discount_curve, model1)

        opt_type = OptionTypes.EUROPEAN_PUT
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bond_option2.value(settle_dt, discount_curve, model2)

        opt_type = OptionTypes.AMERICAN_CALL
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model1 = HWTree(sigma, a, num_time_steps)
        v1call = bond_option1.value(settle_dt, discount_curve, model1)

        opt_type = OptionTypes.EUROPEAN_CALL
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bond_option2.value(settle_dt, discount_curve, model2)

        end = time.time()

        period = end - start

        test_cases.print(period, num_time_steps, v1put, v2put, v1call, v2call)


########################################################################################


def test_bond_option_american_convergence_two():

    # Build discount curve
    settle_dt = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(settle_dt, 0.05)

    # Bond details
    issue_dt = Date(1, 12, 2015)
    maturity_dt = settle_dt.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
    expiry_dt = settle_dt.add_tenor("18m")
    face = 100.0

    spot_value = bond.clean_price_from_discount_curve(settle_dt, discount_curve)
    test_cases.header("LABEL", "VALUE")
    test_cases.print("BOND PRICE", spot_value)

    test_cases.header("TIME", "N", "EUR_CALL", "AMER_CALL", "EUR_PUT", "AMER_PUT")

    sigma = 0.01
    a = 0.1
    hw_model = HWTree(sigma, a)
    k = 102.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    num_steps_vector = range(100, 400, 100)

    for num_steps in num_steps_vector:
        hw_model = HWTree(sigma, a, num_steps)

        start = time.time()

        euro_call_bond_option = BondOption(
            bond, expiry_dt, k, OptionTypes.EUROPEAN_CALL
        )

        v_ec = euro_call_bond_option.value(settle_dt, discount_curve, hw_model)

        amer_call_bond_option = BondOption(
            bond, expiry_dt, k, OptionTypes.AMERICAN_CALL
        )

        v_ac = amer_call_bond_option.value(settle_dt, discount_curve, hw_model)

        euro_put_bond_option = BondOption(bond, expiry_dt, k, OptionTypes.EUROPEAN_PUT)

        v_ep = euro_put_bond_option.value(settle_dt, discount_curve, hw_model)

        amer_put_bond_option = BondOption(bond, expiry_dt, k, OptionTypes.AMERICAN_PUT)

        v_ap = amer_put_bond_option.value(settle_dt, discount_curve, hw_model)

        end = time.time()
        period = end - start

        test_cases.print(period, num_steps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(num_steps_vector, vec_ac, label="American Call")
        plt.legend()

        plt.figure()
        plt.plot(num_steps_vector, vec_ap, label="American Put")
        plt.legend()


########################################################################################


def test_bond_option_zerovol_convergence():

    # Build discount curve
    settle_dt = Date(1, 9, 2019)
    rate = 0.05
    discount_curve = DiscountCurveFlat(settle_dt, rate, FrequencyTypes.ANNUAL)

    # Bond details
    issue_dt = Date(1, 9, 2014)
    maturity_dt = Date(1, 9, 2025)
    coupon = 0.06
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

    # Option Details
    expiry_dt = Date(1, 12, 2021)

    df_expiry = discount_curve.df(expiry_dt)
    fwd_clean_value = bond.clean_price_from_discount_curve(expiry_dt, discount_curve)
    #    fwd_full_value = bond.dirty_price_from_discount_curve(expiry_dt, discount_curve)
    #    print("BOND FwdCleanBondPx", fwd_clean_value)
    #    print("BOND FwdFullBondPx", fwd_full_value)
    #    print("BOND Accrued:", bond.accrued_int)

    spot_clean_value = bond.clean_price_from_discount_curve(settle_dt, discount_curve)

    test_cases.header(
        "STRIKE",
        "STEPS",
        "CALL_INT",
        "CALL_INT_PV",
        "CALL_EUR",
        "CALL_AMER",
        "PUT_INT",
        "PUT_INT_PV",
        "PUT_EUR",
        "PUT_AMER",
    )

    num_time_steps = range(100, 400, 100)
    strike_prices = [90, 120]

    for strike_price in strike_prices:

        call_intrinsic = max(spot_clean_value - strike_price, 0)
        put_intrinsic = max(strike_price - spot_clean_value, 0)
        call_intrinsic_pv = max(fwd_clean_value - strike_price, 0) * df_expiry
        put_intrinsic_pv = max(strike_price - fwd_clean_value, 0) * df_expiry

        for num_steps in num_time_steps:

            sigma = 0.0000001
            a = 0.1
            model = HWTree(sigma, a, num_steps)

            opt_type = OptionTypes.EUROPEAN_CALL
            bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
            v1 = bond_option1.value(settle_dt, discount_curve, model)

            opt_type = OptionTypes.AMERICAN_CALL
            bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)
            v2 = bond_option2.value(settle_dt, discount_curve, model)

            opt_type = OptionTypes.EUROPEAN_PUT
            bond_option3 = BondOption(bond, expiry_dt, strike_price, opt_type)
            v3 = bond_option3.value(settle_dt, discount_curve, model)

            opt_type = OptionTypes.AMERICAN_PUT
            bond_option4 = BondOption(bond, expiry_dt, strike_price, opt_type)
            v4 = bond_option4.value(settle_dt, discount_curve, model)

            test_cases.print(
                strike_price,
                num_steps,
                call_intrinsic,
                call_intrinsic_pv,
                v1,
                v2,
                put_intrinsic,
                put_intrinsic_pv,
                v3,
                v4,
            )


########################################################################################


def test_bond_option_deriva_gem():

    # See https://github.com/domokane/FinancePy/issues/98

    settle_dt = Date(1, 12, 2019)

    rate = 0.05
    dc_type = DayCountTypes.THIRTY_360_BOND
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    discount_curve = DiscountCurveFlat(settle_dt, rate, fixed_freq, dc_type)

    issue_dt = Date(1, 12, 2018)
    expiry_dt = settle_dt.add_tenor("18m")
    maturity_dt = settle_dt.add_tenor("10Y")

    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_360_BOND
    bond = Bond(issue_dt, maturity_dt, coupon, freq_type, accrual_type)
    strike_price = 100.0
    face = 100.0

    euro_call_bond_option = BondOption(
        bond, expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
    )

    # cp = bond.clean_price_from_discount_curve(expiry_dt, discount_curve)
    # fp = bond.dirty_price_from_discount_curve(expiry_dt, discount_curve)
    #    print("Fixed Income Clean Price: %9.3f"% cp)
    #    print("Fixed Income Full  Price: %9.3f"% fp)

    num_steps = 500
    sigma = 0.0125
    a = 0.1
    model_hw = HWTree(sigma, a, num_steps)

    # ec = euro_call_bond_option.value(settle_dt, discount_curve, model_hw)

    coupon_times = []
    coupon_flows = []
    cpn = bond.cpn / bond.freq

    num_flows = len(bond.cpn_dts)
    for i in range(0, num_flows):

        pcd = bond.cpn_dts[i - 1]
        ncd = bond.cpn_dts[i]

        if ncd > settle_dt:

            if len(coupon_times) == 0:
                flow_time = (pcd - settle_dt) / G_DAYS_IN_YEAR
                coupon_times.append(flow_time)
                coupon_flows.append(cpn)

            flow_time = (ncd - settle_dt) / G_DAYS_IN_YEAR
            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    y = 0.05
    times = np.linspace(0, 10, 21)
    dfs = np.power(1 + y / 2, -times * 2)

    sigma = 0.0125
    a = 0.1
    model = HWTree(sigma, a, None)

    #  Test convergence
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR
    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR

    # Jamshidian approach
    v_jam = model.european_bond_option_jamshidian(
        t_exp, strike_price, face, coupon_times, coupon_flows, times, dfs
    )
    # print("Jamshidian:", v_jam)

    model.num_time_steps = 100
    model.build_tree(t_mat, times, dfs)
    exercise_type = ExerciseTypes.EUROPEAN

    v_hw = model.bond_option(
        t_exp, strike_price, face, coupon_times, coupon_flows, exercise_type
    )

    # print("Full Tree:", v_hw)


########################################################################################

test_bond_option_deriva_gem()

test_bond_option_zerovol_convergence()
test_bond_option()
test_bond_option_european_convergence()
test_bond_option_american_convergence_one()
test_bond_option_american_convergence_two()
test_cases.compare_test_cases()
