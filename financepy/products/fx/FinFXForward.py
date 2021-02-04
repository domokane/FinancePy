##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EURUSD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
###############################################################################


class FinFXForward():
    ''' Contract to buy or sell currency at a forward rate decided today. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikeFXRate: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
                 currencyPair: str,  # FORDOM
                 notional: float,
                 notionalCurrency: str,  # must be FOR or DOM
                 spotDays: int = 0):
        ''' Creates a FinFXForward which allows the owner to buy the FOR
        against the DOM currency at the strikeFXRate and to pay it in the
        notional currency. '''

        checkArgumentTypes(self.__init__, locals())

        deliveryDate = expiryDate.addWeekDays(spotDays)

        ''' The FX rate is the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1) '''

        if deliveryDate < expiryDate:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currencyPair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._expiryDate = expiryDate
        self._deliveryDate = deliveryDate
        self._strikeFXRate = strikeFXRate

        self._currencyPair = currencyPair
        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        if notionalCurrency != self._domName and notionalCurrency != self._forName:
            raise FinError("Notional currency not in currency pair.")

        self._notional = notional
        self._notionalCurrency = notionalCurrency
        self._spotDays = spotDays

###############################################################################

    def value(self,
              valueDate,
              spotFXRate,  # 1 unit of foreign in domestic
              domDiscountCurve,
              forDiscountCurve):
        ''' Calculate the value of an FX forward contract where the current
        FX rate is the spotFXRate. '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        newFwdFXRate = self.forward(valueDate,
                                    spotFXRate,
                                    domDiscountCurve,
                                    forDiscountCurve)

        domDF = domDiscountCurve._df(t)

        if self._notionalCurrency == self._domName:
            self._notional_dom = self._notional
            self._notional_for = self._notional / self._strikeFXRate
        elif self._notionalCurrency == self._forName:
            self._notional_dom = self._notional * self._strikeFXRate
            self._notional_for = self._notional
        else:
            raise FinError("Invalid notional currency.")

        if self._notionalCurrency == self._forName:
            v = (newFwdFXRate - self._strikeFXRate)
            v = v * self._notional * domDF
        elif self._notionalCurrency == self._domName:
            v = (newFwdFXRate - self._strikeFXRate)
            v = v * self._notional * domDF * newFwdFXRate

        self._cash_dom = v * self._notional_dom / self._strikeFXRate
        self._cash_for = v * self._notional_for / spotFXRate

        return {"value": v,
                "cash_dom": self._cash_dom,
                "cash_for": self._cash_for,
                "not_dom": self._notional_dom,
                "not_for": self._notional_for,
                "ccy_dom": self._domName,
                "ccy_for": self._forName}

###############################################################################

    def forward(self,
                valueDate,
                spotFXRate,  # 1 unit of foreign in domestic
                domDiscountCurve,
                forDiscountCurve):
        ''' Calculate the FX Forward rate that makes the value of the FX
        contract equal to zero. '''

        if type(valueDate) == FinDate:
            t = (self._deliveryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        forDF = forDiscountCurve._df(t)
        domDF = domDiscountCurve._df(t)

        fwdFXRate = spotFXRate * forDF / domDF
        return fwdFXRate

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE FX RATE", self._strikeFXRate)
        s += labelToString("CURRENCY PAIR", self._currencyPair)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("NOTIONAL CCY", self._notionalCurrency)
        s += labelToString("SPOT DAYS", self._spotDays, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
