##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...utils.fin_math import N
from ...utils.global_variables import gDaysInYear
from ...utils.FinError import FinError
# from ...products.equity.FinEquityOption import FinOption
from ...utils.date import Date
#from ...products.fx.FinFXModelTypes import FinFXModel
from ...models.black_scholes import FinModelBlackScholes
from ...utils.helper_functions import labelToString, check_argument_types
from ...utils.FinGlobalTypes import FinOptionTypes

###############################################################################


class FinFXDigitalOption():

    def __init__(self,
                 expiry_date: Date,
                 strikePrice: float,  # 1 unit of foreign in domestic
                 currencyPair: str,  # FORDOM
                 optionType: FinOptionTypes,
                 notional: float,
                 premCurrency: str):
        """ Create the FX Digital Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strikePrice = float(strikePrice)
        self._currencyPair = currencyPair
        self._optionType = optionType

        self._currencyPair = currencyPair
        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        if premCurrency != self._domName and premCurrency != self._forName:
            raise FinError("Notional currency not in currency pair.")

###############################################################################

    def value(self,
              valuation_date,
              spotFXRate,  # 1 unit of foreign in domestic
              domDiscountCurve,
              forDiscountCurve,
              model):
        """ Valuation of a digital option using Black-Scholes model. This
        allows for 4 cases - first upper barriers that when crossed pay out
        cash (calls) and lower barriers than when crossed from above cause a
        cash payout (puts) PLUS the fact that the cash payment can be in
        domestic or foreign currency. """

        if type(valuation_date) == Date:
            spotDate = valuation_date.addWeekDays(self._spotDays)
            tdel = (self._deliveryDate - spotDate) / gDaysInYear
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            tdel = valuation_date
            texp = tdel

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(tdel < 0.0):
            raise FinError("Option time to maturity is less than zero.")

        tdel = np.maximum(tdel, 1e-10)

        domDF = domDiscountCurve.df(tdel)
        rd = -np.log(domDF)/tdel

        forDF = forDiscountCurve.df(tdel)
        rf = -np.log(forDF)/tdel

        S0 = spotFXRate
        K = self._strikeFXRate

        if type(model) == FinModelBlackScholes:

            volatility = model._volatility
            lnS0k = log(S0 / K)
            den = volatility * sqrt(texp)
            v2 = volatility * volatility
            mu = rd - rf
            d2 = (lnS0k + (mu - v2 / 2.0) * tdel) / den

            if self._optionType == FinOptionTypes.DIGITAL_CALL and \
                self._forName == self._premCurrency:
                    v = S0 * exp(-rf * tdel) * N(d2)
            elif self._optionType == FinOptionTypes.DIGITAL_PUT and \
                self._forName == self._premCurrency:
                    v = S0 * exp(-rf * tdel) * N(-d2)
            if self._optionType == FinOptionTypes.DIGITAL_CALL and \
                self._domName == self._premCurrency:
                    v = exp(-rd * tdel) * N(d2)
            elif self._optionType == FinOptionTypes.DIGITAL_PUT and \
                self._domName == self._premCurrency:
                    v = exp(-rd * tdel) * N(-d2)
            else:
                raise FinError("Unknown option type")

            v = v * self.premNotional

        return v

###############################################################################
