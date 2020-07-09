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
from financepy.models.FinModelRatesHW import FinModelRatesHW

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

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesHW(sigma, a, numTimeSteps)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN CALL - HW", strikePrice, v)

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
        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN CALL - HW", strikePrice, v)

    ###########################################################################

    optionType = FinBondOptionTypes.EUROPEAN_PUT

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("EUROPEAN PUT - HW", strikePrice, v)

    ###########################################################################

    optionType = FinBondOptionTypes.AMERICAN_PUT

    price = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    for strikePrice in strikes:

        sigma = 0.02
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinModelRatesHW(sigma, a)
        v = bondOption.value(settlementDate, discountCurve, model)
        testCases.print("AMERICAN PUT - HW", strikePrice, v)

###############################################################################

def test_FinBondOptionEuropeanConvergence():

    # CONVERGENCE TESTS
    # COMPARE AMERICAN TREE VERSUS JAMSHIDIAN IN EUROPEAN LIMIT TO CHECK THAT
    # TREE HAS BEEN CORRECTLY CONSTRUCTED. FIND VERY GOOD AGREEMENT.

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinFlatCurve(settlementDate, 0.05,
                                 FinFrequencyTypes.CONTINUOUS)

    # Bond details
    maturityDate = FinDate(1, 12, 2020)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    # Option Details - put expiry in the middle of a coupon period
    expiryDate = FinDate(1, 3, 2020)
    strikePrice = 100.0
    face = 100.0

    spotValue = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)

    timeSteps = range(10, 400, 10)
    strikePrice = 100.0

    texp = (expiryDate - settlementDate) / gDaysInYear
    dfExpiry = discountCurve.df(texp)
    tmat = (maturityDate - settlementDate) / gDaysInYear
    dfMat = discountCurve.df(tmat)

    fwdValue = bond.valueBondUsingDiscountCurve(expiryDate, discountCurve)/dfExpiry

    callPV = max(fwdValue - strikePrice,0) * dfExpiry
    putPV = max(strikePrice - fwdValue,0) * dfExpiry

    testCases.header("PERIOD","N","PUT_JAM","PUT_TREE","CALL_JAM","CALL_TREE")

    for numTimeSteps in timeSteps:

        sigma = 0.05
        a = 0.1

        start = time.time()
        optionType = FinBondOptionTypes.EUROPEAN_PUT

        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model1)

        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesHW(sigma, a, numTimeSteps, False)
        v2put = bondOption2.value(settlementDate, discountCurve, model2)

        optionType = FinBondOptionTypes.EUROPEAN_CALL

        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model1)

        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesHW(sigma, a, numTimeSteps, False)
        v2call = bondOption2.value(settlementDate, discountCurve, model2)

        end = time.time()
        period = end - start
        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

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

        sigma = 0.05
        a = 0.1

        start = time.time()

        optionType = FinBondOptionTypes.AMERICAN_PUT
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1put = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinBondOptionTypes.EUROPEAN_PUT
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesHW(sigma, a, numTimeSteps, False)
        v2put = bondOption2.value(settlementDate, discountCurve, model2)

        optionType = FinBondOptionTypes.AMERICAN_CALL
        bondOption1 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model1 = FinModelRatesHW(sigma, a, numTimeSteps)
        v1call = bondOption1.value(settlementDate, discountCurve, model1)

        optionType = FinBondOptionTypes.EUROPEAN_CALL
        bondOption2 = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model2 = FinModelRatesHW(sigma, a, numTimeSteps, False)
        v2call = bondOption2.value(settlementDate, discountCurve, model2)

        end = time.time()

        period = end - start

        testCases.print(period, numTimeSteps, v1put, v2put, v1call, v2call)

###############################################################################

def test_FinBondOptionAmericanConvergenceTWO():

    # Build discount curve
    settlementDate = FinDate(1, 12, 2019)
    discountCurve = FinFlatCurve(settlementDate, 0.05)

    # Bond details
    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)
    expiryDate = settlementDate.addTenor("18m")
    face = 100.0

    spotValue = bond.valueBondUsingDiscountCurve(settlementDate, discountCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("BOND PRICE", spotValue)

    testCases.header("PERIOD","N","EUR_CALL","AMER_CALL","EUR_PUT","AMER_PUT")

    sigma = 0.01
    a = 0.1
    hwModel = FinModelRatesHW(sigma, a)
    K = 102.0

    vec_ec = []
    vec_ac = []
    vec_ep = []
    vec_ap = []

    numStepsVector = range(10,500,10)

    for numSteps in numStepsVector:
        hwModel = FinModelRatesHW(sigma, a, numSteps)

        start = time.time()

        europeanCallBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.EUROPEAN_CALL)
        v_ec = europeanCallBondOption.value(settlementDate, discountCurve, hwModel)

        americanCallBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.AMERICAN_CALL)
        v_ac = americanCallBondOption.value(settlementDate, discountCurve, hwModel)

        europeanPutBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.EUROPEAN_PUT)
        v_ep = europeanPutBondOption.value(settlementDate, discountCurve, hwModel)

        americanPutBondOption = FinBondOption(bond, expiryDate, K, face, FinBondOptionTypes.AMERICAN_PUT)
        v_ap = americanPutBondOption.value(settlementDate, discountCurve, hwModel)

        end = time.time()
        period = end - start

        testCases.print(period, numSteps, v_ec, v_ac, v_ep, v_ap)

        vec_ec.append(v_ec)
        vec_ac.append(v_ac)
        vec_ep.append(v_ep)
        vec_ap.append(v_ap)

    plt.figure()
    plt.plot(numStepsVector, vec_ac, label="American Call")
    plt.legend()

    plt.figure()
    plt.plot(numStepsVector, vec_ap, label="American Put")
    plt.legend()

###############################################################################


test_FinBondOption()
test_FinBondOptionEuropeanConvergence()
test_FinBondOptionAmericanConvergenceONE()
test_FinBondOptionAmericanConvergenceTWO()