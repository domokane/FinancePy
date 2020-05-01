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
from ...products.equities.FinOption import FinOption, FinOptionTypes

##########################################################################
##########################################################################


class FinFXDigitalOption(FinOption):

    def __init__(self,
                 expiryDate,
                 strikePrice,  # ONE UNIT OF FOREIGN IN DOMESTIC CC
                 currencyPair,  # FORDOM
                 optionType
                 notional,
                 notionalCurrency):
        ''' Create the FX Digital Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency. '''

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._currencyPair = currencyPair
        self._optionType = optionType

        self._currencyPair = currencyPair
        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        if notionalCurrency != self._domName and notionalCurrency != self._forName:
            raise FinError("Notional currency not in currency pair.")

###############################################################################
###############################################################################

    def value(self,
              valueDate,
              spotFXRate,  # ONE UNIT OF FOREIGN IN DOMESTIC CCY
              domDiscountCurve,
              forDiscountCurve,
              model):
        ''' Valuation of a digital option using Black-Scholes model. This
        allows for 4 cases - first upper barriers that when crossed pay out
        cash (calls) and lower barriers than when crossed from above cause a
        cash payout (puts) PLUS the fact that the cash payment can be in
        domestic or foreign currency. '''

        if type(valueDate) == FinDate:
            spotDate = valueDate.addWorkDays(self._spotDays)
            tdel = (self._deliveryDate - spotDate) / gDaysInYear
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            tdel = valueDate
            texp = tdel

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if model._parentType != FinFXModel:
            raise FinError("Model is not inherited off type FinFXModel.")

        if np.any(tdel < 0.0):
            raise FinError("Option time to maturity is less than zero.")

        tdel = np.maximum(tdel, 1e-10)

        domDF = domDiscountCurve.df(tdel)
        rd = -np.log(domDF)/tdel

        forDF = forDiscountCurve.df(tdel)
        rf = -np.log(forDF)/tdel

        S0 = spotFXRate
        K = self._strikeFXRate

        if type(model) == FinFXModelBlackScholes \

            volatility = model._volatility
            lnS0k = log(S0 / K)
            den = volatility * sqrt(texp)
            v2 = volatility * volatility
            mu = rd - rf
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinOptionTypes.DIGITAL_CALL
                and self._forName = self._notionalCurrency:
                    v = S0 * exp(-rf * tdel) * N(d2)
            elif self._optionType == FinOptionTypes.DIGITAL_PUT:
                and self._forName = self._notionalCurrency:
                    v = S0 * exp(-rf * tdel) * N(-d2)
            if self._optionType == FinOptionTypes.DIGITAL_CALL
                and self._domName = self._notionalCurrency:
                    v = exp(-rd * tdel) * N(d2)
            elif self._optionType == FinOptionTypes.DIGITAL_PUT:
                and self._domName = self._notionalCurrency:
                    v = exp(-rd * tdel) * N(-d2)
            else:
                raise FinError("Unknown option type")

            v = v * notional

        return v

##########################################################################
##########################################################################
