##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.Date import Date
from ...utils.FinGlobalVariables import gDaysInYear
from ...utils.FinGlobalTypes import FinLongShort
from ...utils.FinError import FinError
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################
# ADD START DATE TO CLASS ?
###############################################################################


class FinEquityForward():
    """ Contract to buy or sell a stock in future at a price agreed today. """

    def __init__(self,
                 expiry_date: Date,
                 forwardPrice: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
                 notional: float,
                 longShort: FinLongShort = FinLongShort.LONG):
        """ Creates a FinEquityForward which allows the owner to buy the stock
        at a price agreed today. Need to specify if LONG or SHORT."""

        checkArgumentTypes(self.__init__, locals())

        self._expiry_date = expiry_date
        self._forwardPrice = forwardPrice
        self._notional = notional
        self._longShort = longShort
        
###############################################################################

    def value(self,
              valuation_date,
              stockPrice,  # Current stock price
              discount_curve,
              dividendCurve):
        """ Calculate the value of an equity forward contract from the stock
        price and discound and dividend curves. """

        if type(valuation_date) == Date:
            t = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        fwdStockPrice = self.forward(valuation_date,
                                     stockPrice,
                                     discount_curve,
                                     dividendCurve)

        discountDF = discount_curve._df(t)

        v = (fwdStockPrice - self._forwardPrice)
        v = v * self._notional * discountDF
        
        if self._longShort == FinLongShort.SHORT:
            v = v * (-1.0)

        return v

###############################################################################

    def forward(self,
                valuation_date,
                stockPrice,  # Current stock price
                discount_curve,
                dividendCurve):
        """ Calculate the forward price of the equity forward contract. """

        if type(valuation_date) == Date:
            t = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(stockPrice <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        discountDF = discount_curve._df(t)
        dividendDF = dividendCurve._df(t)

        fwdStockPrice = stockPrice * dividendDF / discountDF
        return fwdStockPrice

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("FORWARD PRICE", self._forwardPrice)
        s += labelToString("LONG OR SHORT", self._longShort)
        s += labelToString("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################

