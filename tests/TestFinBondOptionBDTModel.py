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
from financepy.models.FinModelRatesBDT import FinModelRatesBDT

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

    tmat = (maturityDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dates = settlementDate.addYears(times)
    dfs = np.exp(-0.05*times)
    discountCurve = FinDiscountCurve(settlementDate, dates, dfs)

    expiryDate = settlementDate.addTenor("18m")
    strikePrice = 105.0
    face = 100.0

    ###########################################################################

    strikes = [80, 90, 100, 110, 120]

    optionType = FinOptionTypes.EUROPEAN_CALL

    testCases.header("LABEL", "VALUE")

    price = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)
    testCases.print("Fixed Income Price:", price)

    numTimeSteps = 100

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN CALL - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN CALL - BK", strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_CALL

    price = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN CALL - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN CALL - BK", strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.EUROPEAN_PUT

    price = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.01

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN PUT - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN PUT - BK", strikePrice, v)

    ###########################################################################

    optionType = FinOptionTypes.AMERICAN_PUT

    price = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.02

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN PUT - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN PUT - BK", strikePrice, v)

###############################################################################


def test_FinBondOptionAmericanConvergenceONE():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinDiscountCurveFlat(settlementDate, 0.05)

    # Bond details
    issueDate = FinDate(1, 9, 2010)
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

    timeSteps = range(30, 100, 1)

    for numTimeSteps in timeSteps:

        sigma = 0.20

        start = time.time()

        optionType = FinOptionTypes.AMERICAN_PUT
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model)

        optionType = FinOptionTypes.EUROPEAN_PUT
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v2put = bondOption2.value(settlementDate, discountCurve, model)

        optionType = FinOptionTypes.AMERICAN_CALL
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model)

        optionType = FinOptionTypes.EUROPEAN_CALL
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBDT(sigma, numTimeSteps)
        v2call = bondOption2.value(settlementDate, discountCurve, model)

        end = time.time()

        period = end - start

        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################


def test_FinBondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinDiscountCurveFlat(settlementDate,
                                         0.05,
                                         FinFrequencyTypes.CONTINUOUS)

    # Bond details
    issueDate = FinDate(1, 9, 2014)
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.05
    freqType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)
    expiryDate = settlementDate.addTenor("18m")
    face = 100.0

    spotValue = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("BOND PRICE", spotValue)

    testCases.header("TIME", "N", "EUR_CALL", "AMER_CALL",
                     "EUR_PUT", "AMER_PUT")

    sigma = 0.2
    model = FinModelRatesBDT(sigma)
    K = 101.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    if 1 == 1:
        K = 100.0
        bkModel = FinModelRatesBDT(sigma, 100)
        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face,
                                               FinOptionTypes.EUROPEAN_CALL)

        v_ec = europeanCallBondOption.value(settlementDate, discountCurve,
                                            model)
        testCases.header("LABEL", "VALUE")
        testCases.print("OPTION", v_ec)

    numStepsVector = range(100, 100, 1)  # should be 100-400

    for numSteps in numStepsVector:

        bkModel = FinModelRatesBDT(sigma, numSteps)

        start = time.time()

        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face,
                                               FinOptionTypes.EUROPEAN_CALL)
        v_ec = europeanCallBondOption.value(settlementDate, discountCurve,
                                            bkModel)

        americanCallBondOption = FinBondOption(bond, expiryDate, K, face,
                                               FinOptionTypes.AMERICAN_CALL)
        v_ac = americanCallBondOption.value(settlementDate, discountCurve,
                                            bkModel)

        europeanPutBondOption = FinBondOption(bond, expiryDate, K, face,
                                              FinOptionTypes.EUROPEAN_PUT)
        v_ep = europeanPutBondOption.value(settlementDate, discountCurve,
                                           bkModel)

        americanPutBondOption = FinBondOption(bond, expiryDate, K, face,
                                              FinOptionTypes.AMERICAN_PUT)
        v_ap = americanPutBondOption.value(settlementDate, discountCurve,
                                           bkModel)

        end = time.time()
        period = end - start

        testCases.print(period, numSteps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    if plotGraphs:

        plt.figure()
        plt.plot(numStepsVector, vec_ec, label="European Call")
        plt.legend()

        plt.figure()
        plt.plot(numStepsVector, vec_ac, label="American Call")
        plt.legend()

        plt.figure()
        plt.plot(numStepsVector, vec_ep, label="European Put")
        plt.legend()

        plt.figure()
        plt.plot(numStepsVector, vec_ap, label="American Put")
        plt.legend()

###############################################################################
###############################################################################

def test_FinBondOptionZEROVOLConvergence():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019) # CHANGED
    rate = 0.05
    discountCurve = FinDiscountCurveFlat(settlementDate, rate, FinFrequencyTypes.ANNUAL)

    # Bond details
    issueDate = FinDate(1, 9, 2015)
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.06
    freqType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, freqType, accrualType)

    # Option Details
    expiryDate = settlementDate.addTenor("18m") # FinDate(1, 12, 2021)
#    print("EXPIRY:", expiryDate)
    face = 100.0

    dfExpiry = discountCurve.df(expiryDate)
    spotCleanValue = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)
    fwdCleanValue = bond.cleanPriceFromDiscountCurve(expiryDate, discountCurve)
#    print("BOND SpotCleanBondPx", spotCleanValue)
#    print("BOND FwdCleanBondPx", fwdCleanValue)
#    print("BOND Accrued:", bond._accruedInterest)

    spotCleanValue = bond.cleanPriceFromDiscountCurve(settlementDate, discountCurve)

    testCases.header("STRIKE", "STEPS",
                     "CALL_INT", "CALL_INT_PV", "CALL_EUR", "CALL_AMER",
                     "PUT_INT", "PUT_INT_PV", "PUT_EUR", "PUT_AMER") 

    numTimeSteps = range(100, 1000, 200)
    strikePrices = [90, 100, 110, 120]

    for strikePrice in strikePrices:
        
        callIntrinsic = max(spotCleanValue - strikePrice, 0)
        putIntrinsic = max(strikePrice - spotCleanValue, 0)
        callIntrinsicPV = max(fwdCleanValue - strikePrice, 0) * dfExpiry
        putIntrinsicPV = max(strikePrice - fwdCleanValue, 0) * dfExpiry

        for numSteps in numTimeSteps:

            sigma = 0.0000001
            model = FinModelRatesBDT(sigma, numSteps)
        
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
# test_FinBondOptionAmericanConvergenceONE()
test_FinBondOptionAmericanConvergenceTWO()
testCases.compareTestCases()
