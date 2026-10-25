##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

import math

from ...utils.date import Date
from ...utils.global_types import LongShortTypes
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.check_values import check_curve_dt, check_stock_price
from ...market.curves.discount_curve import DiscountCurve


class EquityQuantoForward:
    """Quanto equity forward with a fixed currency conversion rate.

    Assumes:
        - Stock price is denominated in foreign currency.
        - FX is quoted as domestic currency per foreign currency.
        - Equity and FX volatilities are annualized.
        - Correlation is between equity and FX log returns.
        - Volatilities and correlation are constant.
        - Dividend and foreign discount curves are foreign-currency
          discount-factor curves.
        - Notional includes the fixed quanto conversion factor
          when that factor differs from 1.
    """

    def __init__(
        self,
        expiry_dt: Date,
        forward_price: float,
        notional: float,
        long_short: LongShortTypes = LongShortTypes.LONG,
    ) -> None:

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.forward_price = forward_price
        self.notional = notional
        self.long_short = long_short

    ###########################################################################

    def _quanto_forward(
        self,
        value_dt: Date,
        stock_price: float,
        foreign_discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        equity_volatility: float,
        fx_volatility: float,
        correlation: float,
    ) -> float:
        """Compute the quanto-adjusted fair forward price."""

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date.")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, foreign_discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        if not math.isfinite(equity_volatility) or equity_volatility < 0.0:
            raise FinError("Equity volatility must be finite and non-negative.")

        if not math.isfinite(fx_volatility) or fx_volatility < 0.0:
            raise FinError("FX volatility must be finite and non-negative.")

        if not math.isfinite(correlation) or not -1.0 <= correlation <= 1.0:
            raise FinError("Correlation must be between -1 and 1.")

        # ACT/365 annualization convention.
        tau = (self.expiry_dt - value_dt) / 365.0

        foreign_df = foreign_discount_curve.df(self.expiry_dt)
        dividend_df = dividend_curve.df(self.expiry_dt)

        if foreign_df <= 0.0 or dividend_df <= 0.0:
            raise FinError("Discount factors must be positive.")

        # Ordinary equity forward under the foreign risk-neutral measure.
        equity_forward = stock_price * dividend_df / foreign_df

        # Quanto drift correction for FX = domestic / foreign.
        quanto_adjustment = math.exp(
            -correlation * equity_volatility * fx_volatility * tau
        )

        return equity_forward * quanto_adjustment

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_discount_curve: DiscountCurve,
        foreign_discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        equity_volatility: float,
        fx_volatility: float,
        correlation: float,
    ) -> float:
        """Present value of the quanto equity forward in domestic currency.

        Notional represents shares multiplied by the fixed quanto
        conversion rate. The forward price is the strike in the
        corresponding payoff units.
        """

        check_curve_dt(value_dt, domestic_discount_curve)

        quanto_forward_price = self._quanto_forward(
            value_dt,
            stock_price,
            foreign_discount_curve,
            dividend_curve,
            equity_volatility,
            fx_volatility,
            correlation,
        )

        domestic_df = domestic_discount_curve.df(self.expiry_dt)

        if domestic_df <= 0.0:
            raise FinError("Domestic discount factor must be positive.")

        v = (
            (quanto_forward_price - self.forward_price)
            * self.notional
            * domestic_df
        )

        if self.long_short == LongShortTypes.SHORT:
            v = -v

        return v

    ###########################################################################

    def forward(
        self,
        value_dt: Date,
        stock_price: float,
        foreign_discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        equity_volatility: float,
        fx_volatility: float,
        correlation: float,
    ) -> float:
        """Calculate the fair quanto-adjusted forward price."""

        return self._quanto_forward(
            value_dt,
            stock_price,
            foreign_discount_curve,
            dividend_curve,
            equity_volatility,
            fx_volatility,
            correlation,
        )

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
