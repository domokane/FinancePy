##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.global_types import FinLongShort
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

###############################################################################
# ADD START DATE TO CLASS ?
###############################################################################


class EquityForward():
    """ Contract to buy or sell a stock in future at a price agreed today. """

    def __init__(self,
                 expiry_date: Date,
                 forward_price: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
                 notional: float,
                 long_short: FinLongShort = FinLongShort.LONG):
        """ Creates a EquityForward which allows the owner to buy the stock
        at a price agreed today. Need to specify if LONG or SHORT."""

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._forward_price = forward_price
        self._notional = notional
        self._long_short = long_short

###############################################################################

    def value(self,
              valuation_date,
              stock_price,  # Current stock price
              discount_curve,
              dividend_curve):
        """ Calculate the value of an equity forward contract from the stock
        price and discount and dividend discount. """

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date")

        if type(valuation_date) == Date:
            t = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        fwdStockPrice = self.forward(valuation_date,
                                     stock_price,
                                     discount_curve,
                                     dividend_curve)

        discountDF = discount_curve._df(t)

        v = (fwdStockPrice - self._forward_price)
        v = v * self._notional * discountDF

        if self._long_short == FinLongShort.SHORT:
            v = v * (-1.0)

        return v

###############################################################################

    def forward(self,
                valuation_date,
                stock_price,  # Current stock price
                discount_curve,
                dividend_curve):
        """ Calculate the forward price of the equity forward contract. """

        if type(valuation_date) == Date:
            t = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        discountDF = discount_curve._df(t)
        dividendDF = dividend_curve._df(t)

        fwdStockPrice = stock_price * dividendDF / discountDF
        return fwdStockPrice

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("FORWARD PRICE", self._forward_price)
        s += label_to_string("LONG OR SHORT", self._long_short)
        s += label_to_string("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
