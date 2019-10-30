# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""
import sys
sys.path.append("..//..")

from financepy.finutils.FinDate import FinDate
from financepy.products.equities.FinBarrierOption import FinBarrierOption
from financepy.products.equities.FinBarrierOption import FinBarrierTypes
from financepy.models.FinProcessSimulator import FinGBMNumericalScheme
from financepy.models.FinProcessSimulator import FinProcessTypes

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinBarrierOption():

    valueDate = FinDate(2015,1,1)
    expiryDate = FinDate(2016,1,1)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    optionType = FinBarrierTypes.DOWN_AND_OUT_CALL

    drift = interestRate - dividendYield
    scheme = FinGBMNumericalScheme.NORMAL
    processType = FinProcessTypes.GBM

    #######################################################################

    import time
    start = time.time()
    numObservationsPerYear = 100

    testCases.header("Type","K","B","S:","Value:","ValueMC","Diff","TIME")   
    
    for optionType in FinBarrierTypes: 
                            
        for stockPrice in range(80,120,10):

            B = 110.0
            K = 100.0
            
            option = FinBarrierOption(expiryDate, K, optionType, B, numObservationsPerYear)
            value = option.value(valueDate,stockPrice, interestRate, dividendYield, volatility)
            start = time.time()
            modelParams = (stockPrice,drift,volatility,scheme)
            valueMC = option.valueMC(valueDate, stockPrice, interestRate, processType, modelParams)

            end = time.time()
            timeElapsed = round(end-start,3)
            diff = valueMC - value
            
            testCases.print(optionType,K,B,stockPrice,value,valueMC,diff,timeElapsed)   

    testCases.header("Type","K","B","S:","Value:","ValueMC","Diff","TIME")   

    for stockPrice in range(80,120,10):

            B = 100.0
            K = 110.0
            
            option = FinBarrierOption(expiryDate, K, optionType, B, numObservationsPerYear)
            value = option.value(valueDate,stockPrice, interestRate, dividendYield, volatility)
            start = time.time() 
            modelParams = (stockPrice,drift,volatility,scheme)
            valueMC = option.valueMC(valueDate, stockPrice, interestRate, processType, modelParams)
            end = time.time()
            timeElapsed = round(end-start,3)
            diff = valueMC - value
            
            testCases.print(optionType,K,B,stockPrice,value,valueMC,diff,timeElapsed)   

    end = time.time()


################################################################################

    stockPrices = range(50,150,10)
    B = 105.0

    testCases.header("Type","K","B","S:","Value","Delta","Vega","Theta")   
    
    for optionType in FinBarrierTypes:

        for stockPrice in stockPrices:

            barrierOption = FinBarrierOption(expiryDate, 100.0, optionType, B, numObservationsPerYear)

            value = barrierOption.value(valueDate,stockPrice, interestRate, dividendYield, volatility)
            delta = barrierOption.delta(valueDate,stockPrice, interestRate, dividendYield, volatility)
            vega  = barrierOption.vega(valueDate,stockPrice, interestRate, dividendYield, volatility)
            theta = barrierOption.theta(valueDate,stockPrice, interestRate, dividendYield, volatility)

            testCases.print(optionType,K,B,stockPrice,value,delta,vega,theta)   
            
test_FinBarrierOption()
testCases.compareTestCases()