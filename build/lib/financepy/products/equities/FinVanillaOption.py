# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from math import exp, log, sqrt
from scipy import optimize

from ...finutils.FinDate import FinDate
from ...finutils.FinMath import N, nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...products.equities.FinOption import FinOption, FinOptionTypes, FinOptionModelTypes

###############################################################################

def f(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    interestRate = args[3]
    dividendYield = args[4]
    modelType = args[5]    
    price = args[6]    
    modelParams = (volatility)

    objFn = self.value(valueDate,
                       stockPrice,
                       interestRate,
                       dividendYield, 
                       modelType, 
                       modelParams) - price
    
#    print(volatility, price, objFn)
    return objFn

###############################################################################

def fvega(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    interestRate = args[3]
    dividendYield = args[4]
    modelType = args[5]

    fprime = self.vega(valueDate,stockPrice,interestRate,dividendYield,modelType,volatility)
    return fprime

###############################################################################

class FinVanillaOption(FinOption):

    def __init__(self,
                 expiryDate,
                 strikePrice,
                 optionType):

        if optionType != FinOptionTypes.EUROPEAN_CALL and optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type", optionType)
            
        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType

################################################################################
            
    def value(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType,
              modelParams):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if stockPrice <= 0.0:
            raise FinError("Stock price must be greater than zero.")

        t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear
        lnS0k = log(stockPrice/ self._strikePrice)
        sqrtT = sqrt(t)

        if modelType == FinOptionModelTypes.BLACKSCHOLES:

            (volatility) = modelParams

            if abs(volatility) < 1e-5:
                volatility = 1e-5

            den = volatility * sqrtT
            v2 = volatility * volatility
            mu = interestRate - dividendYield
            d1 = (lnS0k + (mu + v2 / 2.0) * t)/den
            d2 = (lnS0k + (mu - v2 / 2.0) * t)/den
    
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = stockPrice * exp(-dividendYield * t) * N(d1)
                v = v - self._strikePrice * exp(-interestRate * t) * N(d2)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = self._strikePrice * exp(-interestRate * t) * N(-d2)
                v = v - stockPrice * exp(-dividendYield * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        return v
        
###############################################################################

    def delta(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType,
              modelParams):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if stockPrice <= 0.0:
            raise FinError("Stock price must be greater than zero.")

        if valueDate == self._expiryDate:
            t = 1e-10
        else:
            t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear

        if modelType == FinOptionModelTypes.BLACKSCHOLES:

            (volatility) = modelParams

            if abs(volatility) < 1e-5:
                volatility = 1e-5
    
            lnS0k = log(float(stockPrice) / self._strikePrice)
            sqrtT = sqrt(t)
            mu = interestRate - dividendYield
            d1 = lnS0k + (mu + volatility * volatility / 2.0) * t
            d1 = d1 / volatility / sqrtT
    
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                delta = exp(-dividendYield * t) * N(d1)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                delta = -exp(-dividendYield * t) * N(-d1)
            else:
                print("Unknown option type")
                return None

        return delta

###############################################################################

    def gamma(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType,
              modelParams):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if stockPrice <= 0.0:
            raise FinError("Stock price must be greater than zero.")

        if valueDate == self._expiryDate:
            t = 1e-10
        else:
            t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear

        if modelType == FinOptionModelTypes.BLACKSCHOLES:

            (volatility) = modelParams    

            if abs(volatility) < 1e-5:
                volatility = 1e-5

            lnS0k = log(float(stockPrice) / self._strikePrice)
            mu = interestRate - dividendYield
            d1 = (lnS0k + (mu + volatility*volatility / 2.0) * t)/ volatility / sqrt(t)
            gamma = exp(-dividendYield*t) * nprime(d1) / stockPrice / volatility / sqrt(t)
            return gamma

        else:
            
            return None

###############################################################################

    def vega(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType,
              modelParams):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if stockPrice <= 0.0:
            raise FinError("Stock price must be greater than zero.")

        if valueDate == self._expiryDate:
            t = 1e-10
        else:
            t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear

        if modelType == FinOptionModelTypes.BLACKSCHOLES:

            (volatility) = modelParams    

            if abs(volatility) < 1e-5:
                volatility = 1e-5

            lnS0k = log(float(stockPrice) / self._strikePrice)
            d1 = lnS0k + (interestRate - dividendYield + volatility * volatility / 2.0) * t
            d1 = d1 / (volatility * sqrt(t))
            vega = stockPrice * sqrt(t) * exp(-dividendYield * t) * nprime(d1)
            return vega

        else:
            
            return None

###############################################################################

    def theta(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType,
              modelParams):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if stockPrice <= 0.0:
            raise FinError("Stock price must be greater than zero.")

        if valueDate == self._expiryDate:
            t = 1e-10
        else:
            t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear

        if modelType == FinOptionModelTypes.BLACKSCHOLES:

            (volatility) = modelParams    

            if abs(volatility) < 1e-5:
                volatility = 1e-5
    
            lnS0k = log(float(stockPrice) / self._strikePrice)
            den = volatility * sqrt(t)
            v2 = volatility * volatility
            mu = interestRate - dividendYield
            d1 = (lnS0k + (mu + v2 / 2.0) * t)/den
            d2 = (lnS0k + (mu - v2 / 2.0) * t)/den
            v = 0.0

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = - stockPrice * exp(-dividendYield * t) * nprime(d1) * volatility / 2.0 / sqrt(t)
                v = v - interestRate * self._strikePrice * exp(-interestRate * t) * N(d2)
                v = v + dividendYield * stockPrice * exp(-dividendYield * t) * N(d1)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = - stockPrice * exp(-dividendYield * t) * nprime(d1) * volatility / 2.0 / sqrt(t)
                v = v + interestRate * self._strikePrice * exp(-interestRate * t) * N(-d2)
                v = v - dividendYield * stockPrice * exp(-dividendYield * t) * N(-d1)
            else:
                raise FinError("Unknown option type")
                return 0.0

        return v

###############################################################################

    def impliedVolatility(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              price):

        modelType = FinOptionModelTypes.BLACKSCHOLES
        argtuple = (self,valueDate,stockPrice,interestRate, dividendYield, 
                    modelType, price)

        sigma = optimize.newton(f, x0=0.2, fprime=fvega, args=argtuple, 
                                tol=1e-5, maxiter=50, fprime2=None)
        return sigma

################################################################################
        
    def valueMC(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              modelType, 
              modelParams,
              numPaths = 10000,
              seed = 4242):

        if modelType == FinOptionModelTypes.BLACKSCHOLES:
            (volatility) = modelParams
        else:
            raise FinError("Model Type invalid")
            
        np.random.seed(seed)
        t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear
        mu = interestRate - dividendYield
        v2 = volatility**2
        K = self._strikePrice
        sqrtdt = np.sqrt(t)
        
        # Use Antithetic variables
        g = np.random.normal(0.0,1.0,size=(1,numPaths))
        s = stockPrice * np.exp((mu - v2/2.0) * t)
        m = np.exp(g * sqrtdt * volatility )
        s_1 = s * m
        s_2 = s / m

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff_a_1 = np.maximum(s_1-K,0)
            payoff_a_2 = np.maximum(s_2-K,0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff_a_1 = np.maximum(K-s_1,0)
            payoff_a_2 = np.maximum(K-s_2,0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * exp(-interestRate * t)/2.0
        return v

################################################################################

    def value_MC_OLD(self,
              valueDate,
              stockPrice,
              interestRate,
              dividendYield,
              terminalS,
              seed = 4242):

        self.validate(valueDate, stockPrice, interestRate, dividendYield, 0.1)

        t = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear
        K = self._strikePrice
        
        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            path_payoff = np.maximum(terminalS-K,0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            path_payoff = np.maximum(K-terminalS,0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(path_payoff)
        v = payoff * exp(-interestRate * t)
        return v

################################################################################
