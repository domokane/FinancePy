##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
from ...utils.global_types import FinLongShort
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

###############################################################################
# ADD START DATE TO CLASS ?
###############################################################################


class EquityForward:
    """Contract to buy or sell a stock in future at a price agreed today."""

    def __init__(
        self,
        expiry_dt: Date,
        forward_price: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
        notional: float,
        long_short: FinLongShort = FinLongShort.LONG,
    ):
        """Creates a EquityForward which allows the owner to buy the stock
        at a price agreed today. Need to specify if LONG or SHORT."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.forward_price = forward_price
        self.notional = notional
        self.long_short = long_short

    ###########################################################################

    def value(
        self,
        value_dt,
        stock_price,  # Current stock price
        discount_curve,
        dividend_curve,
    ):
        """Calculate the value of an equity forward contract from the stock
        price and discount and dividend discount."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        fwd_stock_price = self.forward(
            value_dt, stock_price, discount_curve, dividend_curve
        )

        discount_df = discount_curve.df_t(t)

        v = fwd_stock_price - self.forward_price
        v = v * self.notional * discount_df

        if self.long_short == FinLongShort.SHORT:
            v = v * (-1.0)

        return v

    ###########################################################################

    def forward(
        self,
        value_dt,
        stock_price,  # Current stock price
        discount_curve,
        dividend_curve,
    ):
        """Calculate the forward price of the equity forward contract."""

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        discount_df = discount_curve.df_t(t)
        dividend_df = dividend_curve.df_t(t)

        fwd_stock_price = stock_price * dividend_df / discount_df
        return fwd_stock_price

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("FORWARD PRICE", self.forward_price)
        s += label_to_string("LONG OR SHORT", self.long_short)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
