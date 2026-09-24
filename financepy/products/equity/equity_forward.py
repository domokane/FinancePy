##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.date import Date
from ...utils.global_types import LongShortTypes
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price

########################################################################################
# ADD START DATE TO CLASS ?
########################################################################################


class EquityForward:
    """Contract to buy or sell a stock in future at a price agreed today."""

    def __init__(
        self,
        expiry_dt: Date,
        forward_price: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
        notional: float,
        long_short: LongShortTypes = LongShortTypes.LONG,
    ) -> None:
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
        value_dt: Date,
        stock_price: float,  # Current stock price
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
    ):
        """Calculate the value of an equity forward contract from the stock
        price and discount and dividend discount."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        discount_df = discount_curve.df(self.expiry_dt)
        dividend_df = dividend_curve.df(self.expiry_dt)

        mkt_fwd_stock_price = stock_price * dividend_df / discount_df

        v = mkt_fwd_stock_price - self.forward_price
        v = v * self.notional * discount_df

        if self.long_short == LongShortTypes.SHORT:
            v = v * (-1.0)

        return v

    ###########################################################################

    def forward(
        self,
        value_dt: Date,
        stock_price: float,  # Current stock price
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
    ):
        """Calculate the value of an equity forward contract from the stock
        price and discount and dividend discount."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        discount_df = discount_curve.df(self.expiry_dt)
        dividend_df = dividend_curve.df(self.expiry_dt)

        mkt_fwd_stock_price = stock_price * dividend_df / discount_df

        return mkt_fwd_stock_price

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("FORWARD PRICE", self.forward_price)
        s += label_to_string("LONG OR SHORT", self.long_short)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
