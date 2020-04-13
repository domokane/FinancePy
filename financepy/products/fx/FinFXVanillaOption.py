# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm

#from ...finutils.FinMath import N, nprime
from ...finutils.FinDate import FinDate
from ...finutils.FinMath import nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...products.fx.FinFXOption import FinFXOption, FinFXOptionTypes
from ...products.fx.FinFXModelTypes import FinFXModel
from ...products.fx.FinFXModelTypes import FinFXModelBlackScholes

N = norm.cdf

###############################################################################


def f(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    discountCurve = args[3]
    dividendYield = args[4]
    price = args[5]

    model = FinFXModelBlackScholes(volatility)

    objFn = self.value(valueDate,
                       stockPrice,
                       discountCurve,
                       dividendYield,
                       model) - price

#    print(volatility, price, objFn)
    return objFn

###############################################################################


def fvega(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    discountCurve = args[3]
    dividendYield = args[4]

    model = FinFXModelBlackScholes(volatility)
    fprime = self.vega(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    return fprime

###############################################################################


class FinFXVanillaOption(FinFXOption):

    def __init__(self,
                 expiryDate,
                 strikeFXRate, # value of a unit of foreign in domestic currency
                 optionType):

        if optionType != FinFXOptionTypes.EUROPEAN_CALL and \
         optionType != FinFXOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type", optionType)

        self._expiryDate = expiryDate
        self._strikeFXRate = strikeFXRate
        self._optionType = optionType

###############################################################################

    def value(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinEquityModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        domInterestRate = -np.log(domDf)/t

        forDf = forDiscountCurve.df(t)
        forInterestRate = -np.log(forDf)/t

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(spotFXRate / self._strikeFXRate)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = domInterestRate - forInterestRate
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinFXOptionTypes.EUROPEAN_CALL:
                v = spotFXRate * np.exp(-forInterestRate * t) * N(d1)
                v = v - self._strikeFXRate * np.exp(-domInterestRate*t) * N(d2)
            elif self._optionType == FinFXOptionTypes.EUROPEAN_PUT:
                v = self._strikeFXRate * np.exp(-domInterestRate * t) * N(-d2)
                v = v - spotFXRate * np.exp(-forInterestRate * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v

###############################################################################

    def xdelta(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        domInterestRate = -np.log(domDf)/t

        forDf = forDiscountCurve.df(t)
        forInterestRate = -np.log(forDf)/t

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(spotFXRate / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = domInterestRate - forInterestRate
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den

            if self._optionType == FinFXOptionTypes.EUROPEAN_CALL:
                delta = np.exp(-domInterestRate * t) * N(d1)
            elif self._optionType == FinFXOptionTypes.EUROPEAN_PUT:
                delta = -np.exp(-domInterestRate * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        return delta

###############################################################################

    def xgamma(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        domInterestRate = -np.log(domDf)/t

        forDf = forDiscountCurve.df(t)
        forInterestRate = -np.log(forDf)/t

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(spotFXRate / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = domInterestRate - forInterestRate
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            gamma = np.exp(-forInterestRate * t) * nprime(d1) / spotFXRate / den

        else:
            raise FinError("Unknown Model Type")

        return gamma

###############################################################################

    def xvega(self,
             valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinEquityModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        domInterestRate = -np.log(domDf)/t

        forDf = forDiscountCurve.df(t)
        forInterestRate = -np.log(forDf)/t

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(spotFXRate / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = domInterestRate - forInterestRate
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            vega = spotFXRate * sqrtT * np.exp(-forInterestRate * t) * nprime(d1)
        else:
            raise FinError("Unknown Model type")

        return vega

###############################################################################

    def xtheta(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve,
              model):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinEquityModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        domInterestRate = -np.log(domDf)/t

        forDf = forDiscountCurve.df(t)
        forInterestRate = -np.log(forDf)/t

        if type(model) == FinFXModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(spotFXRate / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = domInterestRate - forInterestRate
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinFXOptionTypes.EUROPEAN_CALL:
                v = - spotFXRate * np.exp(-forInterestRate * t) * \
                    nprime(d1) * volatility / 2.0 / sqrtT
                v = v - domInterestRate * self._strikePrice * \
                    np.exp(-domInterestRate * t) * N(d2)
                v = v + forInterestRate * spotFXRate * \
                    np.exp(-forInterestRate * t) * N(d1)
            elif self._optionType == FinFXOptionTypes.EUROPEAN_PUT:
                v = - spotFXRate * np.exp(-forInterestRate * t) * \
                    nprime(d1) * volatility / 2.0 / sqrtT
                v = v + domInterestRate * self._strikePrice * \
                    np.exp(-domInterestRate * t) * N(-d2)
                v = v - forInterestRate * spotFXRate * \
                    np.exp(-forInterestRate * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v

###############################################################################

    def impliedVolatility(self,
                          valueDate,
                          stockPrice,
                          discountCurve,
                          dividendYield,
                          price):

        argtuple = (self, valueDate, stockPrice,
                    discountCurve, dividendYield, price)

        sigma = optimize.newton(f, x0=0.2, fprime=fvega, args=argtuple,
                                tol=1e-5, maxiter=50, fprime2=None)
        return sigma

###############################################################################

    def valueMC(self,
                valueDate,
                spotFXRate,
                domDiscountCurve,
                forDiscountCurve,
                model,
                numPaths=10000,
                seed=4242):

        if model._parentType == FinFXModel:
            volatility = model._volatility
        else:
            raise FinError("Model Type invalid")

        np.random.seed(seed)
        t = (self._expiryDate - valueDate) / gDaysInYear

        domDF = domDiscountCurve.df(self._expiryDate)
        forDF = forDiscountCurve.df(self._expiryDate)

        domRate = -np.log(domDF)/t
        forRate = -np.log(forDF)/t

        mu = domRate - forRate
        v2 = volatility**2
        K = self._strikeFXRate
        sqrtdt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, numPaths))
        s = spotFXRate * np.exp((mu - v2 / 2.0) * t)
        m = np.exp(g * sqrtdt * volatility)
        s_1 = s * m
        s_2 = s / m

        if self._optionType == FinFXOptionTypes.EUROPEAN_CALL:
            payoff_a_1 = np.maximum(s_1 - K, 0)
            payoff_a_2 = np.maximum(s_2 - K, 0)
        elif self._optionType == FinFXOptionTypes.EUROPEAN_PUT:
            payoff_a_1 = np.maximum(K - s_1, 0)
            payoff_a_2 = np.maximum(K - s_2, 0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * np.exp(-domRate * t) / 2.0
        return v

###############################################################################

    def value_MC_OLD(self,
                     valueDate,
                     stockPrice,
                     discountCurve,
                     dividendYield,
                     terminalS,
                     seed=4242):

        self.validate(valueDate, stockPrice, discountCurve, dividendYield, 0.1)

        t = (self._expiryDate - valueDate) / gDaysInYear

        df = discountCurve.df(t)
        r = -np.log(df)/t

        K = self._strikePrice

        if self._optionType == FinFXOptionTypes.EUROPEAN_CALL:
            path_payoff = np.maximum(terminalS - K, 0)
        elif self._optionType == FinFXOptionTypes.EUROPEAN_PUT:
            path_payoff = np.maximum(K - terminalS, 0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(path_payoff)
        v = payoff * np.exp(-r * t)
        return v

###############################################################################
