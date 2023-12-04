##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

###############################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EUR USD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
###############################################################################


class FXForward:
    """ Contract to buy or sell currency at a forward rate decided today. """

    def __init__(self,
                 expiry_dt: Date,
                 strike_fx_rate: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
                 currency_pair: str,  # FOR DOM
                 notional: float,
                 notional_currency: str,  # must be FOR or DOM
                 spot_days: int = 0):
        """ Creates a FinFXForward which allows the owner to buy the FOR
        against the DOM currency at the strike_fx_rate and to pay it in the
        notional currency. """

        check_argument_types(self.__init__, locals())

        delivery_dt = expiry_dt.add_weekdays(spot_days)

        """ The FX rate is the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1) """

        if delivery_dt < expiry_dt:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._expiry_dt = expiry_dt
        self._delivery_dt = delivery_dt
        self._strike_fx_rate = strike_fx_rate

        self._currency_pair = currency_pair
        self._forName = self._currency_pair[0:3]
        self._domName = self._currency_pair[3:6]

        if notional_currency != self._domName and notional_currency != self._forName:
            raise FinError("Notional currency not in currency pair.")

        self._notional = notional
        self._notional_currency = notional_currency
        self._spot_days = spot_days

###############################################################################

    def value(self,
              value_dt,
              spot_fx_rate,  # 1 unit of foreign in domestic
              dom_discount_curve,
              for_discount_curve):
        """ Calculate the value of an FX forward contract where the current
        FX rate is the spot_fx_rate. """

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self._expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if dom_discount_curve._value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as option value date")

        if for_discount_curve._value_dt != value_dt:
            raise FinError(
                "Foreign Curve valuation date not same as option value date")

        if isinstance(value_dt, Date):
            t = (self._expiry_dt - value_dt) / gDaysInYear
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        newFwdFXRate = self.forward(value_dt,
                                    spot_fx_rate,
                                    dom_discount_curve,
                                    for_discount_curve)

        dom_df = dom_discount_curve._df(t)

        if self._notional_currency == self._domName:
            self._notional_dom = self._notional
            self._notional_for = self._notional / self._strike_fx_rate
        elif self._notional_currency == self._forName:
            self._notional_dom = self._notional * self._strike_fx_rate
            self._notional_for = self._notional
        else:
            raise FinError("Invalid notional currency.")

        if self._notional_currency == self._forName:
            v = (newFwdFXRate - self._strike_fx_rate)
            v = v * self._notional * dom_df
        elif self._notional_currency == self._domName:
            v = (newFwdFXRate - self._strike_fx_rate)
            v = v * self._notional * dom_df * newFwdFXRate

        self._cash_dom = v * self._notional_dom / self._strike_fx_rate
        self._cash_for = v * self._notional_for / spot_fx_rate

        return {"value": v,
                "cash_dom": self._cash_dom,
                "cash_for": self._cash_for,
                "not_dom": self._notional_dom,
                "not_for": self._notional_for,
                "ccy_dom": self._domName,
                "ccy_for": self._forName}

###############################################################################

    def forward(self,
                value_dt,
                spot_fx_rate,  # 1 unit of foreign in domestic
                dom_discount_curve,
                for_discount_curve):
        """ Calculate the FX Forward rate that makes the value of the FX
        contract equal to zero. """

        if isinstance(value_dt, Date):
            t = (self._delivery_dt - value_dt) / gDaysInYear
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        for_df = for_discount_curve._df(t)
        dom_df = dom_discount_curve._df(t)

        fwdFXRate = spot_fx_rate * for_df / dom_df
        return fwdFXRate

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_dt)
        s += label_to_string("STRIKE FX RATE", self._strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self._currency_pair)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("NOTIONAL CCY", self._notional_currency)
        s += label_to_string("SPOT DAYS", self._spot_days, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
