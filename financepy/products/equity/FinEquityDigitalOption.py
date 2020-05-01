# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from math import exp, log, sqrt
import numpy as np

from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...products.equity.FinEquityOption import FinEquityOption, FinEquityOptionTypes

##########################################################################
##########################################################################


class FinEquityDigitalOption(FinEquityOption):

    def __init__(self,
                 expiryDate,
                 strikePrice,
                 optionType):

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType

##########################################################################

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if valueDate == self._expiryDate:
            if self._optionType == FinOptionTypes.DIGITAL_CALL:
                return np.heaviside(stockPrice - self._strikePrice, 0)
            elif self._optionType == FinOptionTypes.DIGITAL_PUT:
                return np.heaviside(self._strikePrice - stockPrice, 0)
            else:
                raise FinError("Unknown option type.")

        t = (self._expiryDate - valueDate) / gDaysInYear
        lnS0k = log(float(stockPrice) / self._strikePrice)
        sqrtT = sqrt(t)

        df = discountCurve.df(t)
        r = -log(df)/t

        volatility = model._volatility

        if abs(volatility) < gSmall:
            volatility = gSmall

        den = volatility * sqrtT
        v2 = volatility * volatility
        mu = r - dividendYield
        d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

        if self._optionType == FinEquityOptionTypes.DIGITAL_CALL:
            v = exp(-r * t) * N(d2)
        elif self._optionType == FinEquityOptionTypes.DIGITAL_PUT:
            v = exp(-r * t) * N(-d2)
        else:
            raise FinError("Unknown option type")

        return v

##########################################################################

    def valueMC(self,
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model,
                numPaths=10000,
                seed=4242):

        np.random.seed(seed)
        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(t)
        r = -np.log(df)/t

        volatility = model._volatility

        mu = r - dividendYield
        v2 = volatility**2
        K = self._strikePrice
        sqrtdt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, numPaths))
        s = stockPrice * np.exp((mu - v2 / 2.0) * t)
        m = np.exp(g * sqrtdt * volatility)
        s_1 = s * m
        s_2 = s / m

        if self._optionType == FinEquityOptionTypes.DIGITAL_CALL:
            payoff_a_1 = np.heaviside(s_1 - K, 0)
            payoff_a_2 = np.heaviside(s_2 - K, 0)
        elif self._optionType == FinEquityOptionTypes.DIGITAL_PUT:
            payoff_a_1 = np.heaviside(K - s_1, 0)
            payoff_a_2 = np.heaviside(K - s_2, 0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * exp(-r * t) / 2.0
        return v

##########################################################################
