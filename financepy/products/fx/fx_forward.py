##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.check_values import check_curve_dt
from ...market.curves.discount_curve import DiscountCurve
from ...utils.check_values import check_currency_pair
from ...utils.fx_conventions import fx_spot_date

########################################################################################


def fx_forward_rate(
    value_dt: Date,
    settlement_dt: Date,
    spot_fx_rate: np.ndarray | float,
    currency_pair: str,
    domestic_curve: DiscountCurve,
    foreign_curve: DiscountCurve,
) -> float:
    """Calculate the fair forward exchange rate.

    The spot FX rate is assumed to settle on the
    conventional spot settlement date determined
    from the valuation date and currency pair.
    """
    check_currency_pair(currency_pair)

    s0 = spot_fx_rate

    check_curve_dt(value_dt, domestic_curve, foreign_curve)

    spot_settle_dt = fx_spot_date(currency_pair, value_dt)

    if settlement_dt < spot_settle_dt:
        raise FinError("Forward settlement date before spot settlement date.")

    # Discount factors to the contractual settlement date.
    dom_df = domestic_curve.df(settlement_dt)
    for_df = foreign_curve.df(settlement_dt)

    # Discount factors to the conventional spot settlement date.
    dom_df_spot = domestic_curve.df(spot_settle_dt)
    for_df_spot = foreign_curve.df(spot_settle_dt)

    fwd_rate = (
        s0
        * for_df
        / for_df_spot
        * dom_df_spot
        / dom_df
    )

    return fwd_rate

########################################################################################


class FXForward:
    """Deliverable FX forward transaction.

    The transaction exchanges two currencies at an agreed
    exchange rate on a specified future settlement date.

    The buy_currency determines which currency is received
    and which is paid.

    FX rates are expressed in units of the second currency
    per unit of the first currency.
    """

    def __init__(
        self,
        settlement_dt: Date,
        strike_fx_rate: float,
        currency_pair: str,
        notional: float,
        notional_ccy: str,
        buy_ccy: str,
    ) -> None:

        check_argument_types(self.__init__, locals())

        notional_ccy = notional_ccy.upper()
        buy_ccy = buy_ccy.upper()

        check_currency_pair(currency_pair)

        if not np.isfinite(notional) or notional <= 0.0:
            raise FinError("Notional must be positive and finite.")

        self.for_ccy = currency_pair[0:3]
        self.dom_ccy = currency_pair[3:]

        if self.for_ccy == self.dom_ccy:
            raise FinError("Currencies in pair must be different.")

        if notional_ccy not in (self.for_ccy, self.dom_ccy):
            raise FinError("Notional currency not in currency pair.")

        if buy_ccy not in (self.for_ccy, self.dom_ccy):
            raise FinError("Buy currency not in currency pair.")

        self.settlement_dt = settlement_dt
        self.strike_fx_rate = strike_fx_rate
        self.currency_pair = currency_pair

        self.notional = notional
        self.notional_ccy = notional_ccy
        self.buy_ccy = buy_ccy

        if notional_ccy == self.for_ccy:
            self.notional_for = notional
            self.notional_dom = notional * strike_fx_rate
        else:
            self.notional_dom = notional
            self.notional_for = notional / strike_fx_rate

        # +1: buy first currency, sell second currency.
        # -1: sell first currency, buy second currency.
        self.direction = 1 if buy_ccy == self.for_ccy else -1

        self.cashflow_for = self.direction * self.notional_for
        self.cashflow_dom = -self.direction * self.notional_dom

        self.spot_dt = None

    ##########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
    ) -> dict[str, float]:
        """Calculate the present value of the FX forward in both currencies."""

        if value_dt > self.settlement_dt:
            raise FinError("Valuation date after settlement date.")

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        self.spot_dt = fx_spot_date(self.currency_pair, value_dt)

        if self.settlement_dt < self.spot_dt:
            raise FinError("Settlement date before conventional spot date.")

        fair_forward = fx_forward_rate(
            value_dt,
            self.settlement_dt,
            spot_fx_rate,
            self.currency_pair,
            domestic_curve,
            foreign_curve,
        )

        dom_df = domestic_curve.df(self.settlement_dt)

        npv_dom = dom_df * (self.cashflow_for * fair_forward + self.cashflow_dom)

        # Convert domestic PV using the implied cash FX rate.
        dom_df_spot = domestic_curve.df(self.spot_dt)
        for_df_spot = foreign_curve.df(self.spot_dt)

        cash_fx_rate = spot_fx_rate * dom_df_spot / for_df_spot

        npv_for = npv_dom / cash_fx_rate

        return {
            self.dom_ccy: npv_dom,
            self.for_ccy: npv_for,
        }

    ##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("SPOT DATE", self.spot_dt)
        s += label_to_string("SETTLEMENT DATE", self.settlement_dt)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("NOTIONAL CCY", self.notional_ccy)
        s += label_to_string("BUY CURRENCY", self.buy_ccy)
        return s
