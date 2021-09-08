###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.bk_tree import BKTree
from financepy.utils.global_types import OptionTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date
import numpy as np
import time
import matplotlib.pyplot as plt

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_BondOption():

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    maturity_date = settlement_date.add_tenor("10Y")
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    tmat = (maturity_date - settlement_date) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dates = settlement_date.add_years(times)
    dfs = np.exp(-0.05*times)
    discount_curve = DiscountCurve(settlement_date, dates, dfs)

    expiry_date = settlement_date.add_tenor("18m")
    strike_price = 105.0
    face = 100.0

    ###########################################################################

    strikes = [80, 85, 90, 95, 100, 105, 110, 115, 120]

    option_type = OptionTypes.EUROPEAN_CALL

    testCases.header("LABEL", "VALUE")

    price = bond.full_price_from_discount_curve(
        settlement_date, discount_curve)
    testCases.print("Fixed Income Price:", price)

    num_time_steps = 20

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strike_price in strikes:

        sigma = 0.20
        a = 0.1

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a, num_time_steps)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("EUROPEAN CALL - BK", strike_price, v)

    for strike_price in strikes:

        sigma = 0.20
        a = 0.05

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a, num_time_steps)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("EUROPEAN CALL - BK", strike_price, v)

    ###########################################################################

    option_type = OptionTypes.AMERICAN_CALL

    price = bond.full_price_from_discount_curve(
        settlement_date, discount_curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strike_price in strikes:

        sigma = 0.01
        a = 0.1

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("AMERICAN CALL - BK", strike_price, v)

    for strike_price in strikes:

        sigma = 0.20
        a = 0.05

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("AMERICAN CALL - BK", strike_price, v)

    ###########################################################################

    option_type = OptionTypes.EUROPEAN_PUT

    price = bond.full_price_from_discount_curve(
        settlement_date, discount_curve)

    for strike_price in strikes:

        sigma = 0.01
        a = 0.1

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("EUROPEAN PUT - BK", strike_price, v)

    for strike_price in strikes:

        sigma = 0.20
        a = 0.05

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("EUROPEAN PUT - BK", strike_price, v)

    ###########################################################################

    option_type = OptionTypes.AMERICAN_PUT

    price = bond.full_price_from_discount_curve(
        settlement_date, discount_curve)

    for strike_price in strikes:

        sigma = 0.02
        a = 0.1

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("AMERICAN PUT - BK", strike_price, v)

    for strike_price in strikes:

        sigma = 0.20
        a = 0.05

        bond_option = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model = BKTree(sigma, a)
        v = bond_option.value(settlement_date, discount_curve, model)
        testCases.print("AMERICAN PUT - BK", strike_price, v)

###############################################################################


def test_BondOptionAmericanConvergenceONE():

    # Build discount curve
    settlement_date = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(settlement_date, 0.05)

    # Bond details
    maturity_date = Date(1, 9, 2025)
    issue_date = Date(1, 9, 2016)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    # Option Details
    expiry_date = Date(1, 12, 2020)
    strike_price = 100.0
    face = 100.0

    testCases.header("TIME", "N", "PUT_AMER", "PUT_EUR",
                     "CALL_AME", "CALL_EUR")

    timeSteps = range(30, 100, 10)

    for num_time_steps in timeSteps:

        sigma = 0.20
        a = 0.1

        start = time.time()

        option_type = OptionTypes.AMERICAN_PUT
        bond_option1 = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model1 = BKTree(sigma, a, num_time_steps)
        v1put = bond_option1.value(settlement_date, discount_curve, model1)

        option_type = OptionTypes.EUROPEAN_PUT
        bond_option2 = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model2 = BKTree(sigma, a, num_time_steps)
        v2put = bond_option2.value(settlement_date, discount_curve, model2)

        option_type = OptionTypes.AMERICAN_CALL
        bond_option1 = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model1 = BKTree(sigma, a, num_time_steps)
        v1call = bond_option1.value(settlement_date, discount_curve, model1)

        option_type = OptionTypes.EUROPEAN_CALL
        bond_option2 = BondOption(
            bond, expiry_date, strike_price, face, option_type)
        model2 = BKTree(sigma, a, num_time_steps)
        v2call = bond_option2.value(settlement_date, discount_curve, model2)

        end = time.time()

        period = end - start

        testCases.print(period, num_time_steps, v1put, v2put, v1call, v2call)

###############################################################################


def test_BondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlement_date = Date(1, 12, 2019)
    discount_curve = DiscountCurveFlat(settlement_date,
                                       0.05,
                                       FrequencyTypes.CONTINUOUS)

    # Bond details
    issue_date = Date(1, 9, 2014)
    maturity_date = Date(1, 9, 2025)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)
    expiry_date = settlement_date.add_tenor("18m")
    face = 100.0

    spotValue = bond.full_price_from_discount_curve(
        settlement_date, discount_curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("BOND PRICE", spotValue)

    testCases.header("TIME", "N", "EUR_CALL", "AMER_CALL",
                     "EUR_PUT", "AMER_PUT")

    sigma = 0.2
    a = 0.1
    bkModel = BKTree(sigma, a)
    K = 101.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    if 1 == 1:
        K = 100.0
        bkModel = BKTree(sigma, a, 100)
        europeanCallBondOption = BondOption(bond, expiry_date, K, face,
                                            OptionTypes.EUROPEAN_CALL)

        v_ec = europeanCallBondOption.value(settlement_date, discount_curve,
                                            bkModel)
        testCases.header("LABEL", "VALUE")
        testCases.print("OPTION", v_ec)

    num_stepsVector = range(100, 100, 1)  # should be 100-400

    for num_steps in num_stepsVector:

        bkModel = BKTree(sigma, a, num_steps)

        start = time.time()

        europeanCallBondOption = BondOption(bond, expiry_date, K, face,
                                            OptionTypes.EUROPEAN_CALL)
        v_ec = europeanCallBondOption.value(settlement_date, discount_curve,
                                            bkModel)

        americanCallBondOption = BondOption(bond, expiry_date, K, face,
                                            OptionTypes.AMERICAN_CALL)
        v_ac = americanCallBondOption.value(settlement_date, discount_curve,
                                            bkModel)

        europeanPutBondOption = BondOption(bond, expiry_date, K, face,
                                           OptionTypes.EUROPEAN_PUT)
        v_ep = europeanPutBondOption.value(settlement_date, discount_curve,
                                           bkModel)

        americanPutBondOption = BondOption(bond, expiry_date, K, face,
                                           OptionTypes.AMERICAN_PUT)
        v_ap = americanPutBondOption.value(settlement_date, discount_curve,
                                           bkModel)

        end = time.time()
        period = end - start

        testCases.print(period, num_steps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    if plotGraphs:

        plt.figure()
        plt.plot(num_stepsVector, vec_ec, label="European Call")
        plt.legend()

        plt.figure()
        plt.plot(num_stepsVector, vec_ac, label="American Call")
        plt.legend()

        plt.figure()
        plt.plot(num_stepsVector, vec_ep, label="European Put")
        plt.legend()

        plt.figure()
        plt.plot(num_stepsVector, vec_ap, label="American Put")
        plt.legend()

###############################################################################


def test_BondOptionZEROVOLConvergence():

    # Build discount curve
    settlement_date = Date(1, 9, 2019)
    rate = 0.05
    discount_curve = DiscountCurveFlat(
        settlement_date, rate, FrequencyTypes.ANNUAL)

    # Bond details
    issue_date = Date(1, 9, 2014)
    maturity_date = Date(1, 9, 2025)
    coupon = 0.06
    freq_type = FrequencyTypes.ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    # Option Details
    expiry_date = Date(1, 12, 2021)
    face = 100.0

    dfExpiry = discount_curve.df(expiry_date)
    fwdCleanValue = bond.clean_price_from_discount_curve(
        expiry_date, discount_curve)
    fwdFullValue = bond.full_price_from_discount_curve(
        expiry_date, discount_curve)
#    print("BOND FwdCleanBondPx", fwdCleanValue)
#    print("BOND FwdFullBondPx", fwdFullValue)
#    print("BOND Accrued:", bond._accrued_interest)

    spotCleanValue = bond.clean_price_from_discount_curve(
        settlement_date, discount_curve)

    testCases.header("STRIKE", "STEPS",
                     "CALL_INT", "CALL_INT_PV", "CALL_EUR", "CALL_AMER",
                     "PUT_INT", "PUT_INT_PV", "PUT_EUR", "PUT_AMER")

    num_time_steps = range(100, 1000, 100)
    strike_prices = [90, 100, 110, 120]

    for strike_price in strike_prices:

        callIntrinsic = max(spotCleanValue - strike_price, 0)
        putIntrinsic = max(strike_price - spotCleanValue, 0)
        callIntrinsicPV = max(fwdCleanValue - strike_price, 0) * dfExpiry
        putIntrinsicPV = max(strike_price - fwdCleanValue, 0) * dfExpiry

        for num_steps in num_time_steps:

            sigma = 0.0000001
            a = 0.1
            model = BKTree(sigma, a, num_steps)

            option_type = OptionTypes.EUROPEAN_CALL
            bond_option1 = BondOption(
                bond, expiry_date, strike_price, face, option_type)
            v1 = bond_option1.value(settlement_date, discount_curve, model)

            option_type = OptionTypes.AMERICAN_CALL
            bond_option2 = BondOption(
                bond, expiry_date, strike_price, face, option_type)
            v2 = bond_option2.value(settlement_date, discount_curve, model)

            option_type = OptionTypes.EUROPEAN_PUT
            bond_option3 = BondOption(
                bond, expiry_date, strike_price, face, option_type)
            v3 = bond_option3.value(settlement_date, discount_curve, model)

            option_type = OptionTypes.AMERICAN_PUT
            bond_option4 = BondOption(
                bond, expiry_date, strike_price, face, option_type)
            v4 = bond_option4.value(settlement_date, discount_curve, model)

            testCases.print(strike_price, num_steps,
                            callIntrinsic, callIntrinsicPV, v1, v2,
                            putIntrinsic, putIntrinsicPV, v3, v4)

###############################################################################


test_BondOptionZEROVOLConvergence()
test_BondOption()
test_BondOptionAmericanConvergenceONE()
test_BondOptionAmericanConvergenceTWO()
testCases.compareTestCases()
