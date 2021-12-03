##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from enum import Enum


from ...utils.global_vars import gDaysInYear, gSmall
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...products.equity.equity_option import EquityOption
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve

from ...utils.math import n_vect


###############################################################################


class FinDigitalOptionTypes(Enum):
    CASH_OR_NOTHING = 1,
    ASSET_OR_NOTHING = 2

###############################################################################


class EquityDigitalOption(EquityOption):
    """ A EquityDigitalOption is an option in which the buyer receives some
    payment if the stock price has crossed a barrier ONLY at expiry and zero
    otherwise. There are two types: cash-or-nothing and the asset-or-nothing
    option. We do not care whether the stock price has crossed the barrier
    today, we only care about the barrier at option expiry. For a continuously-
    monitored barrier, use the EquityOneTouchOption class. """

    def __init__(self,
                 expiry_date: Date,
                 barrier_price: float,
                 option_type: OptionTypes,
                 underlying_type: FinDigitalOptionTypes):
        """ Create the digital option by specifying the expiry date, the
        barrier price and the type of option which is either a EUROPEAN_CALL
        or a EUROPEAN_PUT or an AMERICAN_CALL or AMERICAN_PUT. There are two
        types of underlying - cash or nothing and asset or nothing. """

        check_argument_types(self.__init__, locals())

        if option_type != OptionTypes.EUROPEAN_CALL and option_type != OptionTypes.EUROPEAN_PUT:
            raise FinError("Option type must be EUROPEAN CALL or EUROPEAN PUT")

        self._expiry_date = expiry_date
        self._barrier_price = float(barrier_price)
        self._option_type = option_type
        self._underlying_type = underlying_type

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: (float, np.ndarray),
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Digital Option valuation using the Black-Scholes model assuming a
        barrier at expiry. Handles both cash-or-nothing and asset-or-nothing
        options."""

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date")

        t = (self._expiry_date - valuation_date) / gDaysInYear
        t = max(t, 1e-6)

        S0 = stock_price
        X = self._barrier_price
        lnS0k = np.log(S0 / X)

        sqrtT = np.sqrt(t)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/t

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/t

        volatility = model._volatility

        if abs(volatility) < gSmall:
            volatility = gSmall

        d1 = (lnS0k + (r - q + volatility*volatility / 2.0) * t)
        d1 = d1 / volatility / sqrtT
        d2 = d1 - volatility * sqrtT

        if self._underlying_type == FinDigitalOptionTypes.CASH_OR_NOTHING:
            if self._option_type == OptionTypes.EUROPEAN_CALL:
                v = np.exp(-r * t) * n_vect(d2)
            elif self._option_type == OptionTypes.EUROPEAN_PUT:
                v = np.exp(-r * t) * n_vect(-d2)
        elif self._underlying_type == FinDigitalOptionTypes.ASSET_OR_NOTHING:
            if self._option_type == OptionTypes.EUROPEAN_CALL:
                v = S0 * np.exp(-q * t) * n_vect(d1)
            elif self._option_type == OptionTypes.EUROPEAN_PUT:
                v = S0 * np.exp(-q * t) * n_vect(-d1)
        else:
            raise FinError("Unknown underlying type.")

        return v

###############################################################################

    def value_mc(self,
                 valuation_date: Date,
                 stock_price: float,
                 discount_curve: DiscountCurve,
                 dividend_curve: DiscountCurve,
                 model,
                 num_paths: int = 10000,
                 seed: int = 4242):
        """ Digital Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Product assumes a barrier only at expiry. Monte Carlo
        handles both a cash-or-nothing and an asset-or-nothing option."""

        np.random.seed(seed)
        t = (self._expiry_date - valuation_date) / gDaysInYear
        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/t

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/t

        volatility = model._volatility
        K = self._barrier_price
        sqrt_dt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, num_paths))
        s = stock_price * np.exp((r - q - volatility * volatility / 2.0) * t)
        m = np.exp(g * sqrt_dt * volatility)

        s_1 = s * m
        s_2 = s / m

        if self._underlying_type == FinDigitalOptionTypes.CASH_OR_NOTHING:
            if self._option_type == OptionTypes.EUROPEAN_CALL:
                payoff_a_1 = np.heaviside(s_1 - K, 0.0)
                payoff_a_2 = np.heaviside(s_2 - K, 0.0)
            elif self._option_type == OptionTypes.EUROPEAN_PUT:
                payoff_a_1 = np.heaviside(K - s_1, 0.0)
                payoff_a_2 = np.heaviside(K - s_2, 0.0)
        elif self._underlying_type == FinDigitalOptionTypes.ASSET_OR_NOTHING:
            if self._option_type == OptionTypes.EUROPEAN_CALL:
                payoff_a_1 = s_1 * np.heaviside(s_1 - K, 0.0)
                payoff_a_2 = s_2 * np.heaviside(s_2 - K, 0.0)
            elif self._option_type == OptionTypes.EUROPEAN_PUT:
                payoff_a_1 = s_1 * np.heaviside(K - s_1, 0.0)
                payoff_a_2 = s_2 * np.heaviside(K - s_2, 0.0)

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * df / 2.0
        return v

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("BARRIER LEVEL", self._barrier_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("UNDERLYING TYPE", self._underlying_type, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
