##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.math import n_vect  # n_prime_vect

from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
# from ...products.equity.EquityOption import FinOption
from ...utils.date import Date
# from ...products.fx.FinFXModelTypes import FinFXModel
from ...models.black_scholes import BlackScholes
from ...utils.helpers import check_argument_types
from ...utils.global_types import OptionTypes


class FXDoubleDigitalOption:

    def __init__(self,
                 expiry_date: Date,
                 upper_strike: (float, np.ndarray),
                 lower_strike: (float, np.ndarray),
                 currency_pair: str,  # FORDOM
                 notional: float,
                 prem_currency: str,
                 spot_days: int = 0):
        """ Create the FX Double Digital Option object. Inputs include
        expiry date, upper strike, lower strike, currency pair,
        option type notional and the currency of the notional.
        An adjustment for spot days is enabled. All currency rates
        must be entered in the price in domestic currency of one unit
        of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency. """

        check_argument_types(self.__init__, locals())

        delivery_date = expiry_date.add_weekdays(spot_days)

        if delivery_date < expiry_date:
            raise FinError("Delivery date must be on or after expiry date.")

        self._expiry_date = expiry_date
        self._delivery_date = delivery_date

        if np.any(upper_strike < 0.0):
            raise FinError("Negative upper strike.")

        if np.any(lower_strike < 0.0):
            raise FinError("Negative lower strike.")

        self._upper_strike = upper_strike
        self._lower_strike = lower_strike

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._currency_pair = currency_pair
        self._forName = self._currency_pair[0:3]
        self._domName = self._currency_pair[3:6]

        if prem_currency != self._domName and prem_currency != self._forName:
            raise FinError("Notional currency not in currency pair.")

        self._prem_currency = prem_currency

        self._notional = notional

        self._spot_days = spot_days

###############################################################################

    def value(self,
              value_date,
              spot_fx_rate,  # 1 unit of foreign in domestic
              dom_discount_curve,
              for_discount_curve,
              model):
        """ Valuation of a double digital option using Black-Scholes model.
        The option pays out the notional in the premium currency if the
        fx rate is between the upper and lower strike at maturity. The
        valuation is equivalent to the valuation of the difference of
        the value of two digital puts, one with the upper and the other
        with the lower strike """

        if isinstance(value_date, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if dom_discount_curve._value_date != value_date:
            raise FinError(
                "Domestic Curve valuation date not same as valuation date")

        if for_discount_curve._value_date != value_date:
            raise FinError(
                "Foreign Curve valuation date not same as valuation date")

        if isinstance(value_date, Date):
            spot_date = value_date.add_weekdays(self._spot_days)
            tdel = (self._delivery_date - spot_date) / gDaysInYear
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        else:
            tdel = value_date
            t_exp = tdel

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(tdel < 0.0):
            raise FinError("Option time to maturity is less than zero.")

        tdel = np.maximum(tdel, 1e-10)

        # TODO RESOLVE TDEL versus TEXP
        domDF = dom_discount_curve._df(tdel)
        forDF = for_discount_curve._df(tdel)

        r_d = -np.log(domDF) / tdel
        r_f = -np.log(forDF) / tdel

        S0 = spot_fx_rate
        K1 = self._lower_strike
        K2 = self._upper_strike

        if type(model) == BlackScholes:

            volatility = model._volatility
            lnS0k1 = np.log(S0 / K1)
            lnS0k2 = np.log(S0 / K2)
            den = volatility * np.sqrt(t_exp)
            v2 = volatility * volatility
            mu = r_d - r_f
            lower_d2 = (lnS0k1 + (mu - v2 / 2.0) * tdel) / den
            upper_d2 = (lnS0k2 + (mu - v2 / 2.0) * tdel) / den

            if self._prem_currency == self._forName:
                lower_digital = S0 * np.exp(-r_f * tdel) * n_vect(-lower_d2)
                upper_digital = S0 * np.exp(-r_f * tdel) * n_vect(-upper_d2)
            elif self._prem_currency == self._domName:
                lower_digital = np.exp(-r_f * tdel) * n_vect(-lower_d2)
                upper_digital = np.exp(-r_f * tdel) * n_vect(-upper_d2)

            v = (upper_digital - lower_digital) * self._notional

        return v

###############################################################################
