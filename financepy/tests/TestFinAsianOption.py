# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""
import sys
sys.path.append("..//..")

from financepy.finutils.FinDate import FinDate
from financepy.products.equities.FinAsianOption import FinAsianOption
from financepy.products.equities.FinOption import FinOptionTypes
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

testCases = FinTestCases(__file__,globalTestCaseMode)
import time

###############################################################################

testConvergence = False
testTimeEvolution = False
testMCTimings = True

###############################################################################

def testConvergence():

    valueDate = FinDate(2014,1,1)
    startAveragingDate = FinDate(2014,6,1)
    expiryDate = FinDate(2015,1,1)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 120 # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    asianOption = FinAsianOption(startAveragingDate,
                                 expiryDate,
                                 K,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 numObservations)

    testCases.header("K","Geometric", "Turnbull_Wakeman", "Curran", "FastMC", "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    numPathsList = [5000]

    for numPaths in numPathsList:

        accruedAverage = stockPrice * 1.1
        
        valueMC_fast = asianOption.valueMC_fast(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        valueMC_CV = asianOption.valueMC_fast_CV(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        valueGeometric = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "GEOMETRIC",
                                 accruedAverage)

        valueTurnbullWakeman = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "TURNBULL_WAKEMAN",
                                 accruedAverage)

        valueCurran = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "CURRAN",
                                 accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(valueMC_fast)
        valuesMC_CV.append(valueMC_CV)

        testCases.print(numPaths,valueGeometric,valueTurnbullWakeman, valueCurran, valueMC_fast,valueMC_CV)

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

    startAveragingDate = FinDate(2015,1,1)
    expiryDate = FinDate(2016,1,1)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 100 # weekly as we have a year
    accruedAverage = None
    K = 100
    seed = 1976

    asianOption = FinAsianOption(startAveragingDate,
                                 expiryDate,
                                 K,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 numObservations)

    testCases.header("Date","Geometric", "Turnbull_Wakeman", "Curran", "FastMC", "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    valueDates = []
    valueDates.append(FinDate(2014,4,1))
    valueDates.append(FinDate(2014,6,1))
    valueDates.append(FinDate(2014,8,1))
    valueDates.append(FinDate(2015,2,1))
    valueDates.append(FinDate(2015,4,1))
    valueDates.append(FinDate(2015,6,1))
    valueDates.append(FinDate(2015,8,1))
    
    numPaths = 10000
    
    for valueDate in valueDates:

        accruedAverage = stockPrice * 0.9
        
        valueMC_fast = asianOption.valueMC_fast(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        valueMC_CV = asianOption.valueMC_fast_CV(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        valueGeometric = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "GEOMETRIC",
                                 accruedAverage)

        valueTurnbullWakeman = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "TURNBULL_WAKEMAN",
                                 accruedAverage)

        valueCurran = asianOption.value(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 "CURRAN",
                                 accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(valueMC_fast)
        valuesMC_CV.append(valueMC_CV)

        testCases.print(str(valueDate),valueGeometric,valueTurnbullWakeman, valueCurran, valueMC_fast,valueMC_CV)

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
    
################################################################################

def testMCTimings():

    valueDate = FinDate(2014,1,1)
    startAveragingDate = FinDate(2014,6,1)
    expiryDate = FinDate(2015,1,1)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.30
    dividendYield = 0.10
    numObservations = 120 # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    asianOption = FinAsianOption(startAveragingDate,
                                 expiryDate,
                                 K,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 numObservations)

    testCases.header("NUMPATHS","VALUE","TIME","VALUE_MC", "TIME", "VALUE_MC_CV","TIME")

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
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        end = time.time()
        t_MC = end - start
        
        start = time.time()
        valueMC_fast = asianOption.valueMC_fast(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
                                 numPaths,
                                 seed,
                                 accruedAverage)

        end = time.time()
        t_MC_fast = end - start

        start = time.time()
        valueMC_fast_CV = asianOption.valueMC_fast_CV(valueDate,
                                 stockPrice,
                                 dividendYield,
                                 volatility,
                                 interestRate,
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

        testCases.print(numPaths,valueMC,t_MC,valueMC_fast,t_MC_fast,valueMC_fast_CV,t_MC_fast_CV)

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
