###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import time
import numpy as np

import sys

sys.path.append("..")

from financepy.utils.global_types import FinExerciseTypes
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.utils.date import Date
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_types import OptionTypes
from financepy.models.hw_tree import HWTree, FinHWEuropeanCalcType
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

plot_graphs = False

###############################################################################


def test_BondOption():

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
    face = 100.0

    ###########################################################################

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

    ###########################################################################

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

    ###########################################################################

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

    ###########################################################################

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


###############################################################################


def test_BondOptionEuropeanConvergence():

    # CONVERGENCE TESTS
    # COMPARE AMERICAN TREE VERSUS JAMSHIDIAN IN EUROPEAN LIMIT TO CHECK THAT
    # TREE HAS BEEN CORRECTLY CONSTRUCTED. FIND VERY GOOD AGREEMENT.

    # Build discount curve
    settle_dt = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(
        settle_dt, 0.05, FrequencyTypes.CONTINUOUS
    )

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

    test_cases.header(
        "TIME", "N", "PUT_JAM", "PUT_TREE", "CALL_JAM", "CALL_TREE"
    )

    for num_time_steps in time_steps:

        sigma = 0.05
        a = 0.1

        start = time.time()
        opt_type = OptionTypes.EUROPEAN_PUT

        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)
        model1 = HWTree(sigma, a, num_time_steps)
        v1put = bond_option1.value(settle_dt, discount_curve, model1)

        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(
            sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_ONLY
        )
        v2put = bond_option2.value(settle_dt, discount_curve, model2)

        opt_type = OptionTypes.EUROPEAN_CALL

        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model1 = HWTree(sigma, a, num_time_steps)
        v1call = bond_option1.value(settle_dt, discount_curve, model1)

        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(
            sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_TREE
        )
        v2call = bond_option2.value(settle_dt, discount_curve, model2)

        end = time.time()
        period = end - start
        test_cases.print(period, num_time_steps, v1put, v2put, v1call, v2call)


###############################################################################


def test_BondOptionAmericanConvergenceONE():

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

    test_cases.header(
        "TIME", "N", "PUT_AMER", "PUT_EUR", "CALL_AME", "CALL_EUR"
    )

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

        model2 = HWTree(
            sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_ONLY
        )
        v2put = bond_option2.value(settle_dt, discount_curve, model2)

        opt_type = OptionTypes.AMERICAN_CALL
        bond_option1 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model1 = HWTree(sigma, a, num_time_steps)
        v1call = bond_option1.value(settle_dt, discount_curve, model1)

        opt_type = OptionTypes.EUROPEAN_CALL
        bond_option2 = BondOption(bond, expiry_dt, strike_price, opt_type)

        model2 = HWTree(
            sigma, a, num_time_steps, FinHWEuropeanCalcType.EXPIRY_TREE
        )
        v2call = bond_option2.value(settle_dt, discount_curve, model2)

        end = time.time()

        period = end - start

        test_cases.print(period, num_time_steps, v1put, v2put, v1call, v2call)


###############################################################################


def test_BondOptionAmericanConvergenceTWO():

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

    spot_value = bond.clean_price_from_discount_curve(
        settle_dt, discount_curve
    )
    test_cases.header("LABEL", "VALUE")
    test_cases.print("BOND PRICE", spot_value)

    test_cases.header(
        "TIME", "N", "EUR_CALL", "AMER_CALL", "EUR_PUT", "AMER_PUT"
    )

    sigma = 0.01
    a = 0.1
    hwModel = HWTree(sigma, a)
    K = 102.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    num_stepsVector = range(100, 400, 100)

    for num_steps in num_stepsVector:
        hwModel = HWTree(sigma, a, num_steps)

        start = time.time()

        europeanCallBondOption = BondOption(
            bond, expiry_dt, K, OptionTypes.EUROPEAN_CALL
        )

        v_ec = europeanCallBondOption.value(settle_dt, discount_curve, hwModel)

        americanCallBondOption = BondOption(
            bond, expiry_dt, K, OptionTypes.AMERICAN_CALL
        )

        v_ac = americanCallBondOption.value(settle_dt, discount_curve, hwModel)

        europeanPutBondOption = BondOption(
            bond, expiry_dt, K, OptionTypes.EUROPEAN_PUT
        )

        v_ep = europeanPutBondOption.value(settle_dt, discount_curve, hwModel)

        americanPutBondOption = BondOption(
            bond, expiry_dt, K, OptionTypes.AMERICAN_PUT
        )

        v_ap = americanPutBondOption.value(settle_dt, discount_curve, hwModel)

        end = time.time()
        period = end - start

        test_cases.print(period, num_steps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    if plot_graphs:
        plt.figure()
        plt.plot(num_stepsVector, vec_ac, label="American Call")
        plt.legend()

        plt.figure()
        plt.plot(num_stepsVector, vec_ap, label="American Put")
        plt.legend()


###############################################################################


def test_BondOptionZEROVOLConvergence():

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
    fwd_clean_value = bond.clean_price_from_discount_curve(
        expiry_dt, discount_curve
    )
    #    fwd_full_value = bond.dirty_price_from_discount_curve(expiry_dt, discount_curve)
    #    print("BOND FwdCleanBondPx", fwd_clean_value)
    #    print("BOND FwdFullBondPx", fwd_full_value)
    #    print("BOND Accrued:", bond.accrued_int)

    spot_clean_value = bond.clean_price_from_discount_curve(
        settle_dt, discount_curve
    )

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


###############################################################################


def test_BondOptionDerivaGem():

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

    europeanCallBondOption = BondOption(
        bond, expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
    )

    # cp = bond.clean_price_from_discount_curve(expiry_dt, discount_curve)
    # fp = bond.dirty_price_from_discount_curve(expiry_dt, discount_curve)
    #    print("Fixed Income Clean Price: %9.3f"% cp)
    #    print("Fixed Income Full  Price: %9.3f"% fp)

    num_steps = 500
    sigma = 0.0125
    a = 0.1
    modelHW = HWTree(sigma, a, num_steps)

    # ec = europeanCallBondOption.value(settle_dt, discount_curve, modelHW)

    ###########################################################################

    couponTimes = []
    couponFlows = []
    cpn = bond.cpn / bond.freq

    numFlows = len(bond.cpn_dts)
    for i in range(0, numFlows):

        pcd = bond.cpn_dts[i - 1]
        ncd = bond.cpn_dts[i]

        if ncd > settle_dt:

            if len(couponTimes) == 0:
                flowTime = (pcd - settle_dt) / G_DAYS_IN_YEARS
                couponTimes.append(flowTime)
                couponFlows.append(cpn)

            flowTime = (ncd - settle_dt) / G_DAYS_IN_YEARS
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    y = 0.05
    times = np.linspace(0, 10, 21)
    dfs = np.power(1 + y / 2, -times * 2)

    sigma = 0.0125
    a = 0.1
    model = HWTree(sigma, a, None)

    #  Test convergence
    t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEARS
    t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEARS

    # Jamshidian approach
    vjam = model.european_bond_option_jamshidian(
        t_exp, strike_price, face, couponTimes, couponFlows, times, dfs
    )
    # print("Jamshidian:", vjam)

    model.num_time_steps = 100
    model.build_tree(t_mat, times, dfs)
    exerciseType = FinExerciseTypes.EUROPEAN

    vHW = model.bond_option(
        t_exp, strike_price, face, couponTimes, couponFlows, exerciseType
    )

    # print("Full Tree:", vHW)


###############################################################################

test_BondOptionDerivaGem()

test_BondOptionZEROVOLConvergence()
test_BondOption()
test_BondOptionEuropeanConvergence()
test_BondOptionAmericanConvergenceONE()
test_BondOptionAmericanConvergenceTWO()
test_cases.compareTestCases()
