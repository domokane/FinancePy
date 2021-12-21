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

###############################################################################


class FXDigitalOption:

    def __init__(self,
                 expiry_date: Date,
                 strike_fx_rate: (float, np.ndarray),
                 currency_pair: str,  # FORDOM
                 option_type: (OptionTypes, list),
                 notional: float,
                 prem_currency: str,
                 spot_days: int = 0):
        """ Create the FX Digital Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency. """

        check_argument_types(self.__init__, locals())

        delivery_date = expiry_date.add_weekdays(spot_days)

        if delivery_date < expiry_date:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._expiry_date = expiry_date
        self._delivery_date = delivery_date

        if np.any(strike_fx_rate < 0.0):
            raise FinError("Negative strike.")

        self._strike_fx_rate = strike_fx_rate

        self._currency_pair = currency_pair
        self._forName = self._currency_pair[0:3]
        self._domName = self._currency_pair[3:6]

        if prem_currency != self._domName and prem_currency != self._forName:
            raise FinError("Notional currency not in currency pair.")

        self._prem_currency = prem_currency

        self._notional = notional

        if option_type != OptionTypes.DIGITAL_CALL and\
           option_type != OptionTypes.DIGITAL_PUT:
            raise FinError("Unknown Digital Option Type:" + option_type)

        self._option_type = option_type
        self._spot_days = spot_days

###############################################################################

    def value(self,
              valuation_date,
              spot_fx_rate,  # 1 unit of foreign in domestic
              dom_discount_curve,
              for_discount_curve,
              model):
        """ Valuation of a digital option using Black-Scholes model. This
        allows for 4 cases - first upper barriers that when crossed pay out
        cash (calls) and lower barriers than when crossed from above cause a
        cash payout (puts) PLUS the fact that the cash payment can be in
        domestic or foreign currency. """

        if isinstance(valuation_date, Date) is False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if dom_discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Domestic Curve valuation date not same as valuation date")

        if for_discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Foreign Curve valuation date not same as valuation date")

        if type(valuation_date) == Date:
            spot_date = valuation_date.add_weekdays(self._spot_days)
            tdel = (self._delivery_date - spot_date) / gDaysInYear
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            tdel = valuation_date
            texp = tdel

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(tdel < 0.0):
            raise FinError("Option time to maturity is less than zero.")

        tdel = np.maximum(tdel, 1e-10)

        # TODO RESOLVE TDEL versus TEXP
        domDF = dom_discount_curve._df(tdel)
        forDF = for_discount_curve._df(tdel)

        rd = -np.log(domDF) / tdel
        rf = -np.log(forDF) / tdel

        S0 = spot_fx_rate
        K = self._strike_fx_rate

        if type(model) == BlackScholes:

            volatility = model._volatility
            lnS0k = np.log(S0 / K)
            den = volatility * np.sqrt(texp)
            v2 = volatility * volatility
            mu = rd - rf
            d2 = (lnS0k + (mu - v2 / 2.0) * tdel) / den

            if self._option_type == OptionTypes.DIGITAL_CALL and \
                    self._forName == self._prem_currency:
                v = S0 * np.exp(-rf * tdel) * n_vect(d2)
            elif self._option_type == OptionTypes.DIGITAL_PUT and \
                    self._forName == self._prem_currency:
                v = S0 * np.exp(-rf * tdel) * n_vect(-d2)
            elif self._option_type == OptionTypes.DIGITAL_CALL and \
                    self._domName == self._prem_currency:
                v = np.exp(-rd * tdel) * n_vect(d2)
            elif self._option_type == OptionTypes.DIGITAL_PUT and \
                    self._domName == self._prem_currency:
                v = np.exp(-rd * tdel) * n_vect(-d2)
            else:
                raise FinError("Unknown option type")

            v = v * self._notional

        return v

###############################################################################
