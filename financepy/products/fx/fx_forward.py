##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

###############################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EUR USD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
###############################################################################


class FXForward:
    """Contract to buy or sell currency at a forward rate decided today."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_fx_rate: float,  # PRICE OF 1 UNIT OF FOR IN DOM CCY
        currency_pair: str,  # FOR DOM
        notional: float,
        notional_currency: str,  # must be FOR or DOM
        spot_days: int = 0,
    ):
        """Creates a FinFXForward which allows the owner to buy the FOR
        against the DOM currency at the strike_fx_rate and to pay it in the
        notional currency."""

        check_argument_types(self.__init__, locals())

        delivery_dt = expiry_dt.add_weekdays(spot_days)

        """ The FX rate is the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1) """

        if delivery_dt < expiry_dt:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.expiry_dt = expiry_dt
        self.delivery_dt = delivery_dt
        self.strike_fx_rate = strike_fx_rate

        self.currency_pair = currency_pair
        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        if (
            notional_currency != self.dom_name
            and notional_currency != self.for_name
        ):
            raise FinError("Notional currency not in currency pair.")

        self.notional = notional
        self.notional_currency = notional_currency
        self.spot_days = spot_days
        self.notional_dom = None
        self.notional_for = None
        self.cash_dom = None
        self.cash_for = None

    ###########################################################################

    def value(
        self,
        value_dt,
        spot_fx_rate,  # 1 unit of foreign in domestic
        domestic_curve,
        foreign_curve,
    ):
        """Calculate the value of an FX forward contract where the current
        FX rate is the spot_fx_rate."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as option value date"
            )

        if foreign_curve.value_dt != value_dt:
            raise FinError(
                "Foreign Curve valuation date not same as option value date"
            )

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        newfwd_fx_rate = self.forward(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve
        )

        dom_df = domestic_curve.df_t(t)

        if self.notional_currency == self.dom_name:
            self.notional_dom = self.notional
            self.notional_for = self.notional / self.strike_fx_rate
        elif self.notional_currency == self.for_name:
            self.notional_dom = self.notional * self.strike_fx_rate
            self.notional_for = self.notional
        else:
            raise FinError("Invalid notional currency.")

        if self.notional_currency == self.for_name:
            v = newfwd_fx_rate - self.strike_fx_rate
            v = v * self.notional * dom_df
        elif self.notional_currency == self.dom_name:
            v = newfwd_fx_rate - self.strike_fx_rate
            v = v * self.notional * dom_df * newfwd_fx_rate

        self.cash_dom = v * self.notional_dom / self.strike_fx_rate
        self.cash_for = v * self.notional_for / spot_fx_rate

        return {
            "value": v,
            "cash_dom": self.cash_dom,
            "cash_for": self.cash_for,
            "not_dom": self.notional_dom,
            "not_for": self.notional_for,
            "ccy_dom": self.dom_name,
            "ccy_for": self.for_name,
        }

    ###########################################################################

    def forward(
        self,
        value_dt,
        spot_fx_rate,  # 1 unit of foreign in domestic
        domestic_curve,
        foreign_curve,
    ):
        """Calculate the FX Forward rate that makes the value of the FX
        contract equal to zero."""

        if isinstance(value_dt, Date):
            t = (self.delivery_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        for_df = foreign_curve.df_t(t)
        dom_df = domestic_curve.df_t(t)

        fwd_fx_rate = spot_fx_rate * for_df / dom_df
        return fwd_fx_rate

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("NOTIONAL CCY", self.notional_currency)
        s += label_to_string("SPOT DAYS", self.spot_days, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
