###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.products.bonds.FinBondOption import FinBondOption
from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.models.FinModelRatesHW import FinModelRatesHW, FinHWEuropeanCalcType

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_FinBondOption():

    settlementDate = FinDate(1, 12, 2019)
    issueDate = FinDate(1, 12, 2018)
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    times = np.linspace(0, 10.0, 21)
    dfs = np.exp(-0.05*times)
    dates = settlementDate.addYears(times)
    discountCurve = FinDiscountCurve(settlementDate, dates, dfs)

    expiryDate = settlementDate.addTenor("18m")
    strikePrice = 105.0
    face = 100.0

    ###########################################################################

    strikes = [80, 85, 90, 95, 100, 105, 110, 115, 120]

    optionType = FinOptionTypes.EUROPEAN_CALL

    testCases.header("LABEL", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)
    testCases.print("Fixed Income Price:", price)

    numTimeSteps = 100

    testCases.banner("HW EUROPEAN CALL")
    testCases.header("STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_CALL

    price = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    testCases.banner("HW AMERICAN CALL")
    testCases.header("STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.EUROPEAN_PUT
    testCases.banner("HW EUROPEAN PUT")
    testCases.header("STRIKE", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print(strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_PUT
    testCases.banner("HW AMERICAN PUT")
    testCases.header("STRIKE", "VALUE")

    price = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.02
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face,
                                   optionType)

        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print(strikePrice, v)

###############################################################################


def test_FinBondOptionEuropeanConvergence():

    # CONVERGENCE TESTS
    # COMPARE AMERICAN TREE VERSUS JAMSHIDIAN IN EUROPEAN LIMIT TO CHECK THAT
    # TREE HAS BEEN CORRECTLY CONSTRUCTED. FIND VERY GOOD AGREEMENT.

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinDiscountCurveFlat(settlementDate, 0.05,
                                         FinFrequencyTypes.CONTINUOUS)

    # Bond details
    issueDate = FinDate(1, 12, 2015)
    maturityDate = FinDate(1, 12, 2020)
    coupon = 0.05
    freqType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    # Option Details - put expiry in the middle of a coupon period
    expiryDate = FinDate(1, 3, 2020)
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

        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model1)

        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bondOption2.value(settlementDate, discountCurve, model2)

        optionType = FinOptionTypes.EUROPEAN_CALL

        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model1)

        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps,  FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bondOption2.value(settlementDate, discountCurve, model2)

        end = time.time()
        period = end - start
        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################


def test_FinBondOptionAmericanConvergenceONE():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinDiscountCurveFlat(settlementDate, 0.05)

    # Bond details
    issueDate = FinDate(1, 9, 2014)
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    # Option Details
    expiryDate = FinDate(1, 12, 2020)
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
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinOptionTypes.EUROPEAN_PUT
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_ONLY)
        v2put = bondOption2.value(settlementDate, discountCurve, model2)

        optionType = FinOptionTypes.AMERICAN_CALL
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinOptionTypes.EUROPEAN_CALL
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face,
                                    optionType)

        model2 = FinModelRatesHW(sigma, a, numTimeSteps, FinHWEuropeanCalcType.EXPIRY_TREE)
        v2call = bondOption2.value(settlementDate, discountCurve, model2)

        end = time.time()

        period = end - start

        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################


def test_FinBondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinDiscountCurveFlat(settlementDate, 0.05)

    # Bond details
    issueDate = FinDate(1, 12, 2015)
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)
    expiryDate = settlementDate.addTenor("18m")
    face = 100.0

    spotValue = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)
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

        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face,
                                               FinOptionTypes.EUROPEAN_CALL)

        v_ec = europeanCallBondOption.value(settlementDate, discountCurve,
                                            hwModel)

        americanCallBondOption = FinBondOption(bond, expiryDate, K, face,
                                               FinOptionTypes.AMERICAN_CALL)

        v_ac = americanCallBondOption.value(settlementDate, discountCurve,
                                            hwModel)

        europeanPutBondOption = FinBondOption(bond, expiryDate, K, face,
                                              FinOptionTypes.EUROPEAN_PUT)

        v_ep = europeanPutBondOption.value(settlementDate, discountCurve,
                                           hwModel)

        americanPutBondOption = FinBondOption(bond, expiryDate, K, face,
                                              FinOptionTypes.AMERICAN_PUT)

        v_ap = americanPutBondOption.value(settlementDate, discountCurve,
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
    settlementDate = FinDate(1, 9, 2019)
    rate = 0.05
    discountCurve = FinDiscountCurveFlat(settlementDate, rate, FinFrequencyTypes.ANNUAL)

    # Bond details
    issueDate = FinDate(1, 9, 2014)
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.06
    freqType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    # Option Details
    expiryDate = FinDate(1, 12, 2021)
    face = 100.0

    dfExpiry = discountCurve.df(expiryDate)
    fwdCleanValue = bond.cleanPriceFromDiscountCurve(expiryDate, discountCurve)
#    fwdFullValue = bond.fullPriceFromDiscountCurve(expiryDate, discountCurve)
#    print("BOND FwdCleanBondPx", fwdCleanValue)
#    print("BOND FwdFullBondPx", fwdFullValue)
#    print("BOND Accrued:", bond._accruedInterest)

    spotCleanValue = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)

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
            bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)    
            v1 = bondOption1.value(settlementDate, discountCurve, model)
    
            optionType = FinOptionTypes.AMERICAN_CALL
            bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)   
            v2 = bondOption2.value(settlementDate, discountCurve, model)

            optionType = FinOptionTypes.EUROPEAN_PUT
            bondOption3 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)    
            v3 = bondOption3.value(settlementDate, discountCurve, model)
        
            optionType = FinOptionTypes.AMERICAN_PUT
            bondOption4 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)    
            v4 = bondOption4.value(settlementDate, discountCurve, model)
        
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
