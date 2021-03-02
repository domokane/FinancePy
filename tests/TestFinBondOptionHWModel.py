###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.utils.Date import Date
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from financepy.products.bonds.FinBond import FinBond
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.FinGlobalVariables import gDaysInYear
from financepy.products.bonds.FinBondOption import FinBondOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.models.FinModelRatesHW import FinModelRatesHW, FinHWEuropeanCalcType

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_FinBondOption():

    settlement_date = Date(1, 12, 2019)
    issue_date = Date(1, 12, 2018)
    maturity_date = settlement_date.addTenor("10Y")
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    times = np.linspace(0, 10.0, 21)
    dfs = np.exp(-0.05*times)
    dates = settlement_date.addYears(times)
    discount_curve = FinDiscountCurve(settlement_date, dates, dfs)

    expiry_date = settlement_date.addTenor("18m")
    strikePrice = 105.0
    face = 100.0

    ###########################################################################

    strikes = [80, 85, 90, 95, 100, 105, 110, 115, 120]

    optionType = FinOptionTypes.EUROPEAN_CALL

    testCases.header("LABEL", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)
    testCases.print("Fixed Income Price:", price)

    numTimeSteps = 100

    testCases.banner("HW EUROPEAN CALL")
    testCases.header("STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiry_date, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a, numTimeSteps)
        v = bondOption.value(settlement_date, discount_curve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_CALL

    price = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    testCases.banner("HW AMERICAN CALL")
    testCases.header("STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiry_date, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlement_date, discount_curve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.EUROPEAN_PUT
    testCases.banner("HW EUROPEAN PUT")
    testCases.header("STRIKE", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiry_date, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlement_date, discount_curve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_PUT
    testCases.banner("HW AMERICAN PUT")
    testCases.header("STRIKE", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)

    for strikePrice in strikes:

        sigma = 0.02
        a = 0.1

        bondOption = FinBondOption(bond, expiry_date, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlement_date, discount_curve, model)
        testCases.print(strikePrice, v)

###############################################################################


def test_FinBondOptionEuropeanConvergence():

    # CONVERGENCE TESTS
    # COMPARE AMERICAN TREE VERSUS JAMSHIDIAN IN EUROPEAN LIMIT TO CHECK THAT
    # TREE HAS BEEN CORRECTLY CONSTRUCTED. FIND VERY GOOD AGREEMENT.

    # Build discount curve
    settlement_date = Date(1, 12, 2019)
    discount_curve = FinDiscountCurveFlat(settlement_date, 0.05,
                                         FinFrequencyTypes.CONTINUOUS)

    # Bond details
    issue_date = Date(1, 12, 2015)
    maturity_date = Date(1, 12, 2020)
    coupon = 0.05
    freq_type = FinFrequencyTypes.ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    # Option Details - put expiry in the middle of a coupon period
    expiry_date = Date(1, 3, 2020)
    strikePrice = 100.0
    face = 100.0

    timeSteps = range(100, 400, 100)
    strikePrice = 100.0

    testCases.header("TIME","N","PUT_JAM","PUT_TREE","CALL_JAM","CALL_TREE")

    for numTimeSteps in timeSteps:

        sigma = 0.05
        a = 0.1

        start = time.time()
        optionType = FinOptionTypes.EUROPEAN_PUT

        bondOption1 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlement_date, discount_curve, model1)

        bondOption2 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bondOption2.value(settlement_date, discount_curve, model2)

        optionType = FinOptionTypes.EUROPEAN_CALL

        bondOption1 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlement_date, discount_curve, model1)

        bondOption2 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps,  FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bondOption2.value(settlement_date, discount_curve, model2)

        end = time.time()
        period = end - start
        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################


def test_FinBondOptionAmericanConvergenceONE():

    # Build discount curve
    settlement_date = Date(1, 12, 2019)
    discount_curve = FinDiscountCurveFlat(settlement_date, 0.05)

    # Bond details
    issue_date = Date(1, 9, 2014)
    maturity_date = Date(1, 9, 2025)
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    # Option Details
    expiry_date = Date(1, 12, 2020)
    strikePrice = 100.0
    face = 100.0

    testCases.header("TIME", "N", "PUT_AMER", "PUT_EUR",
                     "CALL_AME", "CALL_EUR")

    timeSteps = range(100, 500, 100)

    for numTimeSteps in timeSteps:

        sigma = 0.05
        a = 0.1

        start = time.time()

        optionType = FinOptionTypes.AMERICAN_PUT
        bondOption1 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlement_date, discount_curve, model1)

        optionType = FinOptionTypes.EUROPEAN_PUT
        bondOption2 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bondOption2.value(settlement_date, discount_curve, model2)

        optionType = FinOptionTypes.AMERICAN_CALL
        bondOption1 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlement_date, discount_curve, model1)

        optionType = FinOptionTypes.EUROPEAN_CALL
        bondOption2 = FinBondOption(bond, expiry_date, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bondOption2.value(settlement_date, discount_curve, model2)

        end = time.time()

        period = end - start

        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################


def test_FinBondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlement_date = Date(1, 12, 2019)
    discount_curve = FinDiscountCurveFlat(settlement_date, 0.05)

    # Bond details
    issue_date = Date(1, 12, 2015)
    maturity_date = settlement_date.addTenor("10Y")
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)
    expiry_date = settlement_date.addTenor("18m")
    face = 100.0

    spotValue = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)
    testCases.header("LABEL", "VALUE")
    testCases.print("BOND PRICE", spotValue)

    testCases.header("TIME", "N", "EUR_CALL", "AMER_CALL",
                     "EUR_PUT", "AMER_PUT")

    sigma = 0.01
    a = 0.1
    hwModel = FinModelRatesHW(sigma, a)
    K = 102.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    numStepsVector = range(100, 500, 100)

    for numSteps in numStepsVector:
        hwModel = FinModelRatesHW(sigma, a, numSteps)

        start = time.time()

        europeanCallBondOption = FinBondOption(bond, expiry_date, K, face,
                                               FinOptionTypes.EUROPEAN_CALL)

        v_ec = europeanCallBondOption.value(settlement_date, discount_curve,
                                            hwModel)

        americanCallBondOption = FinBondOption(bond, expiry_date, K, face,
                                               FinOptionTypes.AMERICAN_CALL)

        v_ac = americanCallBondOption.value(settlement_date, discount_curve,
                                            hwModel)

        europeanPutBondOption = FinBondOption(bond, expiry_date, K, face,
                                              FinOptionTypes.EUROPEAN_PUT)

        v_ep = europeanPutBondOption.value(settlement_date, discount_curve,
                                           hwModel)

        americanPutBondOption = FinBondOption(bond, expiry_date, K, face,
                                              FinOptionTypes.AMERICAN_PUT)

        v_ap = americanPutBondOption.value(settlement_date, discount_curve,
                                           hwModel)

        end = time.time()
        period = end - start

        testCases.print(period, numSteps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    if plotGraphs:
        plt.figure()
        plt.plot(numStepsVector, vec_ac, label="American Call")
        plt.legend()

        plt.figure()
        plt.plot(numStepsVector, vec_ap, label="American Put")
        plt.legend()


###############################################################################

def test_FinBondOptionZEROVOLConvergence():

    # Build discount curve
    settlement_date = Date(1, 9, 2019)
    rate = 0.05
    discount_curve = FinDiscountCurveFlat(settlement_date, rate, FinFrequencyTypes.ANNUAL)

    # Bond details
    issue_date = Date(1, 9, 2014)
    maturity_date = Date(1, 9, 2025)
    coupon = 0.06
    freq_type = FinFrequencyTypes.ANNUAL
    accrual_type = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issue_date, maturity_date, coupon, freq_type, accrual_type)

    # Option Details
    expiry_date = Date(1, 12, 2021)
    face = 100.0

    dfExpiry = discount_curve.df(expiry_date)
    fwdCleanValue = bond.cleanPriceFromDiscountCurve(expiry_date, discount_curve)
#    fwdFullValue = bond.fullPriceFromDiscountCurve(expiry_date, discount_curve)
#    print("BOND FwdCleanBondPx", fwdCleanValue)
#    print("BOND FwdFullBondPx", fwdFullValue)
#    print("BOND Accrued:", bond._accruedInterest)

    spotCleanValue = bond.cleanPriceFromDiscountCurve(settlement_date, discount_curve)

    testCases.header("STRIKE", "STEPS",
                     "CALL_INT", "CALL_INT_PV", "CALL_EUR", "CALL_AMER",
                     "PUT_INT", "PUT_INT_PV", "PUT_EUR", "PUT_AMER") 

    numTimeSteps = range(100, 1000, 100)
    strikePrices = [90, 100, 110, 120]

    for strikePrice in strikePrices:
        
        callIntrinsic = max(spotCleanValue - strikePrice, 0)
        putIntrinsic = max(strikePrice - spotCleanValue, 0)
        callIntrinsicPV = max(fwdCleanValue - strikePrice, 0) * dfExpiry
        putIntrinsicPV = max(strikePrice - fwdCleanValue, 0) * dfExpiry

        for numSteps in numTimeSteps:

            sigma = 0.0000001
            a = 0.1
            model = FinModelRatesHW(sigma, a, numSteps)
        
            optionType = FinOptionTypes.EUROPEAN_CALL
            bondOption1 = FinBondOption(bond, expiry_date, strikePrice, face, optionType)    
            v1 = bondOption1.value(settlement_date, discount_curve, model)
    
            optionType = FinOptionTypes.AMERICAN_CALL
            bondOption2 = FinBondOption(bond, expiry_date, strikePrice, face, optionType)   
            v2 = bondOption2.value(settlement_date, discount_curve, model)

            optionType = FinOptionTypes.EUROPEAN_PUT
            bondOption3 = FinBondOption(bond, expiry_date, strikePrice, face, optionType)    
            v3 = bondOption3.value(settlement_date, discount_curve, model)
        
            optionType = FinOptionTypes.AMERICAN_PUT
            bondOption4 = FinBondOption(bond, expiry_date, strikePrice, face, optionType)    
            v4 = bondOption4.value(settlement_date, discount_curve, model)
        
            testCases.print(strikePrice, numSteps,
                            callIntrinsic, callIntrinsicPV, v1, v2,
                            putIntrinsic, putIntrinsicPV, v3, v4)

###############################################################################

test_FinBondOptionZEROVOLConvergence()
test_FinBondOption()
test_FinBondOptionEuropeanConvergence()
test_FinBondOptionAmericanConvergenceONE()
test_FinBondOptionAmericanConvergenceTWO()
testCases.compareTestCases()
