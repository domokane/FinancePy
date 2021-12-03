##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from enum import Enum

from ...utils.math import N
from ...utils.global_vars import gDaysInYear, gSmall
from ...utils.error import FinError
from ...models.gbm_process_simulator import FinGBMProcess
from ...products.fx.fx_option import FXOption
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.global_types import OptionTypes
from ...market.curves.discount_curve import DiscountCurve


##########################################################################
# TODO: Attempt control variate adjustment to monte carlo
# TODO: Sobol for Monte Carlo
# TODO: TIGHTEN UP LIMIT FOR W FROM 100
# TODO: Vectorise the analytical pricing formula
##########################################################################


##########################################################################
# FLOAT STRIKE LOOKBACK CALL PAYS MAX(S(T)-SMIN,0)
# FLOAT STRIKE LOOKBACK PUT PAYS MAX(SMAX-S(T),0)
##########################################################################


class FXFloatLookbackOption(FXOption):
    """ This is an FX option in which the strike of the option is not fixed
    but is set at expiry to equal the minimum fx rate in the case of a call
    or the maximum fx rate in the case of a put. """

    def __init__(self,
                 expiry_date: Date,
                 option_type: OptionTypes):
        """ Create the FX Float Look Back Option by specifying the expiry
        date and the option type. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._option_type = option_type

    ##########################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              domestic_curve: DiscountCurve,
              foreign_curve: DiscountCurve,
              volatility: float,
              stock_min_max: float):
        """ Valuation of the Floating Lookback option using Black-Scholes
        using the formulae derived by Goldman, Sosin and Gatto (1979). """

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve._valuation_date != valuation_date:
            raise FinError(
                "Domestic Curve valuation date not same as option valuation date")

        if foreign_curve._valuation_date != valuation_date:
            raise FinError(
                "Foreign Curve valuation date not same as option valuation date")

        t = (self._expiry_date - valuation_date) / gDaysInYear

        df = domestic_curve._df(t)
        r = -np.log(df) / t

        dq = foreign_curve._df(t)
        q = -np.log(dq) / t

        v = volatility
        s0 = stock_price
        smin = 0.0
        smax = 0.0

        if self._option_type == OptionTypes.EUROPEAN_CALL:
            smin = stock_min_max
            if smin > s0:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._option_type == OptionTypes.EUROPEAN_PUT:
            smax = stock_min_max
            if smax < s0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        if abs(r - q) < gSmall:
            q = r + gSmall

        dq = np.exp(-q * t)
        df = np.exp(-r * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = np.exp(b * t)

        # Taken from Haug Page 142
        if self._option_type == OptionTypes.EUROPEAN_CALL:

            a1 = (np.log(s0 / smin) + (b + (v ** 2) / 2.0) * t) / v / np.sqrt(t)
            a2 = a1 - v * np.sqrt(t)

            if smin == s0:
                term = N(-a1 + 2.0 * b * np.sqrt(t) / v) - expbt * N(-a1)
            elif s0 < smin and w < -100:
                term = - expbt * N(-a1)
            else:
                term = ((s0 / smin) ** (-w)) * N(-a1 + 2.0 *
                                                 b * np.sqrt(t) / v) - expbt * N(-a1)

            v = s0 * dq * N(a1) - smin * df * N(a2) + s0 * df * u * term

        elif self._option_type == OptionTypes.EUROPEAN_PUT:

            b1 = (np.log(s0 / smax) + (b + (v ** 2) / 2.0) * t) / v / np.sqrt(t)
            b2 = b1 - v * np.sqrt(t)

            if smax == s0:
                term = -N(b1 - 2.0 * b * np.sqrt(t) / v) + expbt * N(b1)
            elif s0 < smax and w > 100:
                term = expbt * N(b1)
            else:
                term = (-(s0 / smax) ** (-w)) * \
                    N(b1 - 2.0 * b * np.sqrt(t) / v) + expbt * N(b1)

            v = smax * df * N(-b2) - s0 * dq * N(-b1) + s0 * df * u * term

        else:
            raise FinError("Unknown lookback option type:" +
                           str(self._option_type))

        return v

    ##########################################################################

    def value_mc(self,
                 valuation_date,
                 stock_price,
                 domestic_curve,
                 foreign_curve,
                 volatility,
                 stock_min_max,
                 num_paths=10000,
                 num_steps_per_year=252,
                 seed=4242):

        t = (self._expiry_date - valuation_date) / gDaysInYear
        df = domestic_curve._df(t)
        r = -np.log(df) / t

        dq = domestic_curve._df(t)
        q = -np.log(dq) / t

        num_time_steps = int(t * num_steps_per_year)
        mu = r - q

        option_type = self._option_type
        smin = 0.0
        smax = 0.0

        if self._option_type == OptionTypes.EUROPEAN_CALL:
            smin = stock_min_max
            if smin > stock_price:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._option_type == OptionTypes.EUROPEAN_PUT:
            smax = stock_min_max
            if smax < stock_price:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        model = FinGBMProcess()
        Sall = model.get_paths(
            num_paths,
            num_time_steps,
            t,
            mu,
            stock_price,
            volatility,
            seed)

        # Due to anti-thetics we have doubled the number of paths
        num_paths = 2 * num_paths
        payoff = np.zeros(num_paths)

        if option_type == OptionTypes.EUROPEAN_CALL:
            SMin = np.min(Sall, axis=1)
            SMin = np.minimum(SMin, smin)
            payoff = np.maximum(Sall[:, -1] - SMin, 0.0)
        elif option_type == OptionTypes.EUROPEAN_PUT:
            SMax = np.max(Sall, axis=1)
            SMax = np.maximum(SMax, smax)
            payoff = np.maximum(SMax - Sall[:, -1], 0.0)
        else:
            raise FinError("Unknown lookback option type:" + str(option_type))

        v = payoff.mean() * np.exp(-r * t)
        return v

##########################################################################
