# -*- coding: utf-8 -*-

import numpy as np
import time

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinFlatCurve import FinFlatCurve

from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.products.bonds.FinBondOption import FinBondOption
from financepy.products.bonds.FinBondOption import FinBondOptionTypes
from financepy.models.FinModelRatesBK import FinModelRatesBK

import matplotlib.pyplot as plt

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBondOption():

    settlementDate = FinDate(1, 12, 2019)

    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    tmat = (maturityDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dfs = np.exp(-0.05*times)
    discountCurve = FinDiscountCurve(settlementDate, times, dfs)

    expiryDate = settlementDate.addTenor("18m")
    strikePrice = 105.0
    face = 100.0

    ###########################################################################

    strikes = [80, 85, 90, 95, 100, 105, 110, 115, 120]

    optionType = FinBondOptionTypes.EUROPEAN_CALL

    testCases.header("LABEL", "VALUE")

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    testCases.print("Fixed Income Price:", price)

    numTimeSteps = 100

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.20
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN CALL - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20
        a = 0.05

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN CALL - BK", strikePrice, v)

    ###########################################################################

    optionType = FinBondOptionTypes.AMERICAN_CALL

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("Fixed Income Price:", price)

    testCases.header("OPTION TYPE AND MODEL", "STRIKE", "VALUE")

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN CALL - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20
        a = 0.05

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN CALL - BK", strikePrice, v)

    ###########################################################################

    optionType = FinBondOptionTypes.EUROPEAN_PUT

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN PUT - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20
        a = 0.05

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN PUT - BK", strikePrice, v)

    ###########################################################################

    optionType = FinBondOptionTypes.AMERICAN_PUT

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.02
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN PUT - BK", strikePrice, v)

    for strikePrice in strikes:

        sigma = 0.20
        a = 0.05

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesBK(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN PUT - BK", strikePrice, v)

###############################################################################


def test_FinBondOptionAmericanConvergenceONE():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinFlatCurve(settlementDate, 0.05)

    # Bond details
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    # Option Details
    expiryDate = FinDate(1, 12, 2020)
    strikePrice = 100.0
    face = 100.0

    spotValue = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    texp = (expiryDate - settlementDate) / gDaysInYear
    dfExpiry = discountCurve.df(texp)
    tmat = (maturityDate - settlementDate) / gDaysInYear
    dfMat = discountCurve.df(tmat)

    fwdValue = bond.valueBondUsingDiscountCurve(expiryDate, discountCurve)/dfExpiry

    callPV = max(fwdValue - strikePrice, 0) * dfExpiry
    putPV = max(strikePrice - fwdValue, 0) * dfExpiry

    testCases.header("PERIOD","N","PUT_AMER","PUT_EUR","CALL_AME","CALL_EUR")

    timeSteps = range(10, 100, 1)

    for numTimeSteps in timeSteps:

        sigma = 0.20
        a = 0.1

        start = time.time()

        optionType = FinBondOptionTypes.AMERICAN_PUT
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesBK(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinBondOptionTypes.EUROPEAN_PUT
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesBK(sigma, a, numTimeSteps)
        v2put = bondOption2.value(settlementDate, discountCurve, model2)

        optionType = FinBondOptionTypes.AMERICAN_CALL
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesBK(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinBondOptionTypes.EUROPEAN_CALL
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesBK(sigma, a, numTimeSteps)
        v2call = bondOption2.value(settlementDate, discountCurve, model2)

        end = time.time()

        period = end - start

        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################

def test_FinBondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinFlatCurve(settlementDate, 0.05, -1)

    # Bond details
    maturityDate = FinDate(1, 9, 2025)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)
    expiryDate = settlementDate.addTenor("18m")
    face = 100.0

    spotValue = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("BOND PRICE", spotValue)

    testCases.header("PERIOD","N","EUR_CALL","AMER_CALL","EUR_PUT","AMER_PUT")

    sigma = 0.2
    a = 0.1
    bkModel = FinModelRatesBK(sigma, a)
    K = 101.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    if 1==1:
        K = 100.0
        bkModel = FinModelRatesBK(sigma, a, 100)
        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.EUROPEAN_CALL)
        v_ec = europeanCallBondOption.value(settlementDate, discountCurve, bkModel)
        print("OPTION", v_ec)

    numStepsVector = range(100, 100, 1) # should be 100-400

    for numSteps in numStepsVector:

        bkModel = FinModelRatesBK(sigma, a, numSteps)

        start = time.time()

        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.EUROPEAN_CALL)
        v_ec = europeanCallBondOption.value(settlementDate, discountCurve, bkModel)

        americanCallBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.AMERICAN_CALL)
        v_ac = americanCallBondOption.value(settlementDate, discountCurve, bkModel)

        europeanPutBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.EUROPEAN_PUT)
        v_ep = europeanPutBondOption.value(settlementDate, discountCurve, bkModel)

        americanPutBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.AMERICAN_PUT)
        v_ap = americanPutBondOption.value(settlementDate, discountCurve, bkModel)

        end = time.time()
        period = end - start

        testCases.print(period, numSteps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

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


#test_FinBondOption()
#test_FinBondOptionAmericanConvergenceONE()
test_FinBondOptionAmericanConvergenceTWO()