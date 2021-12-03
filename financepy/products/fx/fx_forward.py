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
                 expiry_date: Date,
                 strike_fx_rate: float,  # PRICE OF 1 UNIT OF FOREIGN IN DOM CCY
                 currency_pair: str,  # FOR DOM
                 notional: float,
                 notional_currency: str,  # must be FOR or DOM
                 spot_days: int = 0):
        """ Creates a FinFXForward which allows the owner to buy the FOR
        against the DOM currency at the strike_fx_rate and to pay it in the
        notional currency. """

        check_argument_types(self.__init__, locals())

        delivery_date = expiry_date.add_weekdays(spot_days)

        """ The FX rate is the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1) """

        if delivery_date < expiry_date:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._expiry_date = expiry_date
        self._delivery_date = delivery_date
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
              valuation_date,
              spot_fx_rate,  # 1 unit of foreign in domestic
              dom_discount_curve,
              for_discount_curve):
        """ Calculate the value of an FX forward contract where the current
        FX rate is the spot_fx_rate. """

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if dom_discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Domestic Curve valuation date not same as option valuation date")

        if for_discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Foreign Curve valuation date not same as option valuation date")

        if type(valuation_date) == Date:
            t = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        newFwdFXRate = self.forward(valuation_date,
                                    spot_fx_rate,
                                    dom_discount_curve,
                                    for_discount_curve)

        domDF = dom_discount_curve._df(t)

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
            v = v * self._notional * domDF
        elif self._notional_currency == self._domName:
            v = (newFwdFXRate - self._strike_fx_rate)
            v = v * self._notional * domDF * newFwdFXRate

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
                valuation_date,
                spot_fx_rate,  # 1 unit of foreign in domestic
                dom_discount_curve,
                for_discount_curve):
        """ Calculate the FX Forward rate that makes the value of the FX
        contract equal to zero. """

        if type(valuation_date) == Date:
            t = (self._delivery_date - valuation_date) / gDaysInYear
        else:
            t = valuation_date

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        forDF = for_discount_curve._df(t)
        domDF = dom_discount_curve._df(t)

        fwdFXRate = spot_fx_rate * forDF / domDF
        return fwdFXRate

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
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
