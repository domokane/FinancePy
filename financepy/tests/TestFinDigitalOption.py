# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import sys
sys.path.append("..//..")

from financepy.finutils.FinDate import FinDate
from financepy.products.equities.FinDigitalOption import FinDigitalOption
from financepy.products.equities.FinOption import FinOptionTypes

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

################################################################################

def test_FinDigitalOption():

    valueDate = FinDate(2015,1,1)
    expiryDate = FinDate(2016,1,1)
    stockPrice = 100.0
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01

    import time

    callOptionValues = []
    callOptionValuesMC = []
    numPathsList = [10000,20000,40000,80000,160000,320000,640000,1280000,2560000]

    testCases.header("NumLoops","ValueBS","ValueMC","TIME")
    
    for numPaths in numPathsList:

        callOption = FinDigitalOption(expiryDate, 100.0, FinOptionTypes.DIGITAL_CALL)
        value = callOption.value(valueDate,stockPrice, dividendYield, volatility, interestRate)
        start = time.time()
        valueMC = callOption.valueMC(valueDate,stockPrice, dividendYield, volatility, interestRate, numPaths)
        end = time.time()
        duration = end-start
        testCases.print(numPaths,value,valueMC,duration)

        callOptionValues.append(value)
        callOptionValuesMC.append(valueMC)

#    plt.figure(figsize=(10,8))
#    plt.plot(numPathsList, callOptionValues, color = 'b', label="Call Option")
#    plt.plot(numPathsList, callOptionValuesMC, color = 'r', label = "Call Option MC")
#    plt.xlabel("Num Loops")
#    plt.legend(loc='best')
    
################################################################################

    stockPrices = range(50,150)
    callOptionValues = []
    callOptionDeltas = []
    callOptionVegas = []
    callOptionThetas = []

    for stockPrice in stockPrices:
        callOption = FinDigitalOption(expiryDate, 100.0, FinOptionTypes.DIGITAL_CALL)
        value = callOption.value(valueDate,stockPrice, dividendYield, volatility, interestRate)
        delta = callOption.delta(valueDate,stockPrice, dividendYield, volatility, interestRate)
        vega = callOption.vega(valueDate,stockPrice, dividendYield, volatility, interestRate)
        theta = callOption.theta(valueDate,stockPrice, dividendYield, volatility, interestRate)
        callOptionValues.append(value)
        callOptionDeltas.append(delta)
        callOptionVegas.append(vega)
        callOptionThetas.append(theta)

    putOptionValues = []
    putOptionDeltas = []
    putOptionVegas = []
    putOptionThetas = []

    for stockPrice in stockPrices:
        putOption = FinDigitalOption(expiryDate, 100.0, FinOptionTypes.DIGITAL_PUT)
        value = putOption.value(valueDate,stockPrice, dividendYield, volatility, interestRate)
        delta = putOption.delta(valueDate,stockPrice, dividendYield, volatility, interestRate)
        vega = putOption.vega(valueDate,stockPrice, dividendYield, volatility, interestRate)
        theta = putOption.theta(valueDate,stockPrice, dividendYield, volatility, interestRate)
        putOptionValues.append(value)
        putOptionDeltas.append(delta)
        putOptionVegas.append(vega)
        putOptionThetas.append(theta)

################################################################################

test_FinDigitalOption()
testCases.compareTestCases()