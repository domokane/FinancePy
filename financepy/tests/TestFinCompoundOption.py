# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate

from financepy.products.equities.FinOption import FinOptionTypes
from financepy.products.equities.FinCompoundOption import FinCompoundOption

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

################################################################################                

def test_FinCompoundOption():

    valueDate = FinDate(2015,1,1)
    expiryDate1 = FinDate(2017,1,1)
    expiryDate2 = FinDate(2018,1,1)
    k1 = 5.0
    k2 = 95.0
    stockPrice = 85.0
    volatility = 0.15
    interestRate = 0.035
    dividendYield = 0.01

    optionType1 = FinOptionTypes.EUROPEAN_CALL
    optionType2 = FinOptionTypes.EUROPEAN_PUT

    numStepsList = [200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000]

    cmpdOption = FinCompoundOption(expiryDate1, expiryDate2, k1, k2, optionType1, optionType2)
    stockPrice = 85.0       

    testCases.header("TYPE1","TYPE2","K1","K2","S","Exact","TreeSteps","TreeValue")                   
    for numSteps in numStepsList:
    
        value = cmpdOption.value(valueDate, stockPrice, interestRate, dividendYield, volatility)
        values = cmpdOption.valueTree(valueDate, stockPrice, interestRate, dividendYield, volatility, numSteps)
        testCases.print(optionType1,optionType2,k1,k2,stockPrice,value,numSteps,values[0])

    testCases.header("TYPE1","TYPE2","K1","K2","S","Exact","TreeSteps","TreeValue","Diff","DELTA","GAMMA","THETA")

    for optionType1 in [FinOptionTypes.EUROPEAN_CALL,FinOptionTypes.EUROPEAN_PUT]:
        for optionType2 in [FinOptionTypes.EUROPEAN_CALL,FinOptionTypes.EUROPEAN_PUT]:
                    
            cmpdOption = FinCompoundOption(expiryDate1, expiryDate2, k1, k2, optionType1, optionType2)
            stockPrices = range(70,100)
        
            for stockPrice in stockPrices:
                value = cmpdOption.value(valueDate, stockPrice, interestRate, dividendYield, volatility)
                delta = cmpdOption.delta(valueDate, stockPrice, interestRate, dividendYield, volatility)
                vega = cmpdOption.vega(valueDate, stockPrice, interestRate, dividendYield, volatility)
                theta = cmpdOption.theta(valueDate, stockPrice, interestRate, dividendYield, volatility)
                values = cmpdOption.valueTree(valueDate, stockPrice, interestRate, dividendYield, volatility)
                diff = value - values[0]
                testCases.print(optionType1,optionType2,k1,k2,stockPrice,value,numSteps,values[0],diff,delta,vega,theta)
                
################################################################################                

test_FinCompoundOption()
testCases.compareTestCases()
