# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""
import time
import numpy as np

from financepy.finutils.FinDate import FinDate
from financepy.products.equities.FinOption import FinOptionTypes, FinOptionModelTypes
from financepy.products.equities.FinVanillaOption import FinVanillaOption
from financepy.products.equities.FinBinomialTree import FinBinomialTree, FinTreeExerciseTypes, FinTreePayoffTypes
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)


def test_FinBinomialTree():

    stockPrice = 50.0
    riskFreeRate = 0.06
    dividendYield = 0.04
    volatility = 0.40
                
    strikePrice = 50.0
    valueDate = FinDate(2016,1,1)
    expiryDate = FinDate(2017,1,1)
    
    numStepsList = [100,500,1000,2000,5000]
    
    modelType = FinOptionModelTypes.BLACKSCHOLES
    modelParams = (volatility)

    testCases.banner("================== EUROPEAN PUT =======================")
    
    putOption = FinVanillaOption(expiryDate, strikePrice, FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(valueDate, stockPrice,  riskFreeRate, dividendYield, modelType, modelParams)
    delta = putOption.delta(valueDate, stockPrice,  riskFreeRate, dividendYield, modelType, modelParams)
    gamma = putOption.gamma(valueDate, stockPrice,  riskFreeRate, dividendYield, modelType, modelParams)
    theta = putOption.theta(valueDate, stockPrice,  riskFreeRate, dividendYield, modelType, modelParams)
    testCases.header("BS Value","BS Delta","BS Gamma","BS Theta")
    testCases.print(value,delta,gamma,theta)
    
    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.EUROPEAN
    params = np.array([-1,strikePrice])

    testCases.header("NumSteps","Results","TIME")
    
    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
        results = tree.value(stockPrice, riskFreeRate, dividendYield, volatility,
                             numSteps, valueDate, payoff, expiryDate, payoff, exercise, params)    
        end = time.time()
        duration = end-start
        testCases.print(numSteps,results,duration)
    
    testCases.banner("================== AMERICAN PUT =======================")
    
    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.AMERICAN
    params = np.array([-1,strikePrice])
    
    testCases.header("NumSteps","Results","TIME")

    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
        results = tree.value(stockPrice, riskFreeRate, dividendYield, volatility,
                             numSteps, valueDate, payoff, expiryDate, payoff, exercise, params)    
        end = time.time()
        duration = end-start
        testCases.print(numSteps,results,duration)
        
    testCases.banner("================== EUROPEAN CALL =======================")
    
    callOption = FinVanillaOption(expiryDate, strikePrice, FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valueDate, stockPrice, riskFreeRate, dividendYield, modelType, modelParams)
    delta = callOption.delta(valueDate, stockPrice, riskFreeRate, dividendYield, modelType, modelParams)
    gamma = callOption.gamma(valueDate, stockPrice, riskFreeRate, dividendYield, modelType, modelParams)
    theta = callOption.theta(valueDate, stockPrice, riskFreeRate, dividendYield, modelType, modelParams)
    testCases.header("BS Value","BS Delta","BS Gamma","BS Theta")
    testCases.print(value,delta,gamma,theta)
    
    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.EUROPEAN
    params = np.array([1.0,strikePrice])
    
    testCases.header("NumSteps","Results","TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
    
        results = tree.value(stockPrice, riskFreeRate, dividendYield, volatility,
                             numSteps, valueDate, payoff, expiryDate, payoff, exercise, params)    
        
        end = time.time()
        duration = end-start
        testCases.print(numSteps,results,duration)
    
    testCases.banner("================== AMERICAN CALL =======================")
    
    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.AMERICAN
    params = np.array([1.0,strikePrice])
    
    testCases.header("NumSteps","Results","TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
    
        results = tree.value(stockPrice, riskFreeRate, dividendYield, volatility,
                             numSteps, valueDate, payoff, expiryDate, payoff, exercise, params)    
        
        end = time.time()
        duration = end-start
        testCases.print(numSteps,results,duration)
 
test_FinBinomialTree()
testCases.compareTestCases()
