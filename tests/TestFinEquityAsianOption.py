###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityAsianOption import FinEquityAsianOption
from financepy.products.equity.FinEquityAsianOption import FinAsianOptionValuationMethods
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

testConvergence = False
testTimeEvolution = False
testMCTimings = True

###############################################################################


def testConvergence():

    valueDate = FinDate(1, 1, 2014)
    startAveragingDate = FinDate(1, 6, 2014)
    expiryDate = FinDate(1, 1, 2015)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 120  # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield)

    asianOption = FinEquityAsianOption(startAveragingDate,
                                       expiryDate,
                                       K,
                                       FinOptionTypes.EUROPEAN_CALL,
                                       numObservations)

    testCases.header(
        "K",
        "Geometric",
        "Turnbull_Wakeman",
        "Curran",
        "FastMC",
        "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    numPathsList = [5000]

    for numPaths in numPathsList:

        accruedAverage = stockPrice * 1.1

        valueMC_fast = asianOption._valueMC_fast(valueDate,
                                                 stockPrice,
                                                 discountCurve,
                                                 dividendCurve,
                                                 model,
                                                 numPaths,
                                                 seed,
                                                 accruedAverage)

        valueMC_CV = asianOption.valueMC(valueDate,
                                         stockPrice,
                                         discountCurve,
                                         dividendCurve,
                                         model,
                                         numPaths,
                                         seed,
                                         accruedAverage)

        valueGeometric = asianOption.value(valueDate,
                                           stockPrice,
                                           discountCurve,
                                           dividendCurve,
                                           model,
                                           FinAsianOptionValuationMethods.GEOMETRIC,
                                           accruedAverage)

        valueTurnbullWakeman = asianOption.value(valueDate,
                                                 stockPrice,
                                                 discountCurve,
                                                 dividendCurve,
                                                 model,
                                                 FinAsianOptionValuationMethods.TURNBULL_WAKEMAN,
                                                 accruedAverage)

        valueCurran = asianOption.value(valueDate,
                                        stockPrice,
                                        discountCurve,
                                        dividendCurve,
                                        model,
                                        FinAsianOptionValuationMethods.CURRAN,
                                        accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(valueMC_fast)
        valuesMC_CV.append(valueMC_CV)

        testCases.print(
            numPaths,
            valueGeometric,
            valueTurnbullWakeman,
            valueCurran,
            valueMC_fast,
            valueMC_CV)

#    import matplotlib.pyplot as plt
#    x = numPathsList
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesGeometric,label="Geometric")
#    plt.plot(x,valuesTurnbull,label="Turbull_Wakeman")
#    plt.plot(x,valuesCurran,label="Curran")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_CV,label="MC_CV")
#    plt.legend()
#    plt.xlabel("Number of Paths")
#    plt.show()

###############################################################################


def testTimeEvolution():

    startAveragingDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 100  # weekly as we have a year
    accruedAverage = None
    K = 100
    seed = 1976

    model = FinModelBlackScholes(volatility)

    asianOption = FinEquityAsianOption(startAveragingDate,
                                       expiryDate,
                                       K,
                                       FinOptionTypes.EUROPEAN_CALL,
                                       numObservations)

    testCases.header(
        "Date",
        "Geometric",
        "Turnbull_Wakeman",
        "Curran",
        "FastMC",
        "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    valueDates = []
    valueDates.append(FinDate(1, 4, 2014))
    valueDates.append(FinDate(1, 6, 2014))
    valueDates.append(FinDate(1, 8, 2014))
    valueDates.append(FinDate(1, 2, 2015))
    valueDates.append(FinDate(1, 4, 2015))
    valueDates.append(FinDate(1, 6, 2015))
    valueDates.append(FinDate(1, 8, 2015))

    numPaths = 10000

    for valueDate in valueDates:

        accruedAverage = stockPrice * 0.9

        discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
        dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield)

        valueMC_fast = asianOption._valueMC_fast(valueDate,
                                                 stockPrice,
                                                 discountCurve,
                                                 dividendCurve,
                                                 model,
                                                 numPaths,
                                                 seed,
                                                 accruedAverage)

        valueMC_CV = asianOption.valueMC(valueDate,
                                         stockPrice,
                                         discountCurve,
                                         dividendCurve,
                                         model,
                                         numPaths,
                                         seed,
                                         accruedAverage)

        valueGeometric = asianOption.value(valueDate,
                                           stockPrice,
                                           discountCurve,
                                           dividendCurve,
                                           model,
                                           FinAsianOptionValuationMethods.GEOMETRIC,
                                           accruedAverage)

        valueTurnbullWakeman = asianOption.value(valueDate,
                                                 stockPrice,
                                                 discountCurve,
                                                 dividendCurve,
                                                 model,
                                                 FinAsianOptionValuationMethods.TURNBULL_WAKEMAN,
                                                 accruedAverage)

        valueCurran = asianOption.value(valueDate,
                                        stockPrice,
                                        discountCurve,
                                        dividendCurve,
                                        model,
                                        FinAsianOptionValuationMethods.CURRAN,
                                        accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(valueMC_fast)
        valuesMC_CV.append(valueMC_CV)

        testCases.print(
            str(valueDate),
            valueGeometric,
            valueTurnbullWakeman,
            valueCurran,
            valueMC_fast,
            valueMC_CV)

#    import matplotlib.pyplot as plt
#    x = [ dt.date() for dt in valueDates]
#
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesGeometric,label="Geometric")
#    plt.plot(x,valuesTurnbull,label="Turbull_Wakeman")
#    plt.plot(x,valuesCurran,label="Curran")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_CV,label="MC_CV")
#    plt.legend()
#    plt.xlabel("Valuation Date")
#    plt.show()

##########################################################################


def testMCTimings():

    valueDate = FinDate(1, 1, 2014)
    startAveragingDate = FinDate(1, 6, 2014)
    expiryDate = FinDate(1, 1, 2015)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 120  # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield)

    asianOption = FinEquityAsianOption(startAveragingDate,
                                       expiryDate,
                                       K,
                                       FinOptionTypes.EUROPEAN_CALL,
                                       numObservations)

    testCases.header(
        "NUMPATHS",
        "VALUE",
        "TIME",
        "VALUE_MC",
        "TIME",
        "VALUE_MC_CV",
        "TIME")

    valuesMC = []
    valuesMC_fast = []
    valuesMC_fast_CV = []

    tvaluesMC = []
    tvaluesMC_fast = []
    tvaluesMC_fast_CV = []

    numPathsList = [5000]

    for numPaths in numPathsList:

        accruedAverage = stockPrice * 1.1

        start = time.time()
        valueMC = asianOption.valueMC(valueDate,
                                      stockPrice,
                                      discountCurve,
                                      dividendCurve,
                                      model,
                                      numPaths,
                                      seed,
                                      accruedAverage)

        end = time.time()
        t_MC = end - start

        start = time.time()
        valueMC_fast = asianOption._valueMC_fast(valueDate,
                                                 stockPrice,
                                                 discountCurve,
                                                 dividendCurve,
                                                 model,
                                                 numPaths,
                                                 seed,
                                                 accruedAverage)

        end = time.time()
        t_MC_fast = end - start

        start = time.time()
        valueMC_fast_CV = asianOption.valueMC(valueDate,
                                              stockPrice,
                                              discountCurve,
                                              dividendCurve,
                                              model,
                                              numPaths,
                                              seed,
                                              accruedAverage)

        end = time.time()
        t_MC_fast_CV = end - start

        valuesMC.append(valueMC)
        valuesMC_fast.append(valueMC_fast)
        valuesMC_fast_CV.append(valueMC_fast_CV)

        tvaluesMC.append(t_MC)
        tvaluesMC_fast.append(t_MC_fast)
        tvaluesMC_fast_CV.append(t_MC_fast_CV)

        testCases.print(
            numPaths,
            valueMC,
            t_MC,
            valueMC_fast,
            t_MC_fast,
            valueMC_fast_CV,
            t_MC_fast_CV)

#    import matplotlib.pyplot as plt
#    x = numPathsList
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesMC,label="Basic MC")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_fast_CV,label="MC_Fast CV")
#    plt.legend()
#    plt.xlabel("Number of Paths")
#    plt.show()


testConvergence()
testMCTimings()
testTimeEvolution()
testCases.compareTestCases()
