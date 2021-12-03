##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...utils.math import N
from ...utils.global_vars import gDaysInYear, gSmall
from ...utils.error import FinError
from ...models.gbm_process_simulator import FinGBMProcess
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
# FIXED STRIKE LOOKBACK CALL PAYS MAX(SMAX-K,0)
# FIXED STRIKE LOOKBACK PUT PAYS MAX(K-SMIN,0)
##########################################################################


class FXFixedLookbackOption:
    """ The Class for FX Fixed Strike Lookback options. """

    def __init__(self,
                 expiry_date: Date,
                 option_type: OptionTypes,
                 optionStrike: float):
        """ Create option with expiry date, option type and the option strike
        """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._option_type = option_type
        self._optionStrike = optionStrike

##########################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              dom_discount_curve: DiscountCurve,
              for_discount_curve: DiscountCurve,
              volatility: float,
              stock_min_max: float):
        """ Value FX Fixed Lookback Option using Black Scholes model and
        analytical formulae. """

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

        t = (self._expiry_date - valuation_date) / gDaysInYear

        df = dom_discount_curve.df(self._expiry_date)
        r = -np.log(df)/t

        dq = for_discount_curve.df(self._expiry_date)
        q = -np.log(dq)/t

        v = volatility
        s0 = stock_price
        k = self._optionStrike
        smin = 0.0
        smax = 0.0

        if self._option_type == OptionTypes.EUROPEAN_CALL:
            smax = stock_min_max
            if smax < s0:
                raise FinError(
                    "The Smax value must be >= the stock price.")
        elif self._option_type == OptionTypes.EUROPEAN_PUT:
            smin = stock_min_max
            if smin > s0:
                raise FinError(
                    "The Smin value must be <= the stock price.")

        # There is a risk of an overflow in the limit of q=r which
        # we remove by adjusting the value of the dividend
        if abs(r - q) < gSmall:
            q = r + gSmall

        df = exp(-r * t)
        dq = exp(-q * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = exp(b * t)

        # Taken from Hull Page 536 (6th edition) and Haug Page 143
        if self._option_type == OptionTypes.EUROPEAN_CALL:

            if k > smax:
                d1 = (log(s0/k) + (b+v*v/2.0)*t)/v/sqrt(t)
                d2 = d1 - v * sqrt(t)

                if s0 == k:
                    term = -N(d1 - 2.0 * b * sqrt(t) / v) + expbt * N(d1)
                elif s0 < k and w > 100:
                    term = expbt * N(d1)
                else:
                    term = (-(s0 / k)**(-w)) * N(d1 - 2.0 *
                                                 b * sqrt(t) / v) + expbt * N(d1)

                v = s0 * dq * N(d1) - k * df * N(d2) + s0 * df * u * term

            else:
                e1 = (log(s0/smax) + (b+v*v/2.0)*t) / v / sqrt(t)
                e2 = e1 - v * sqrt(t)

                if s0 == smax:
                    term = -N(e1 - 2.0 * b * sqrt(t) / v) + expbt * N(e1)
                elif s0 < smax and w > 100:
                    term = expbt * N(e1)
                else:
                    term = (-(s0 / smax)**(-w)) * \
                        N(e1 - 2.0 * b * sqrt(t) / v) + expbt * N(e1)

                v = df * (smax - k) + s0 * dq * N(e1) - \
                    smax * df * N(e2) + s0 * df * u * term

        elif self._option_type == OptionTypes.EUROPEAN_PUT:

            if k >= smin:
                f1 = (log(s0 / smin) + (b + v * v / 2.0) * t) / v / sqrt(t)
                f2 = f1 - v * sqrt(t)

                if s0 == smin:
                    term = N(-f1 + 2.0 * b * sqrt(t) / v) - expbt * N(-f1)
                elif s0 > smin and w < -100:
                    term = -expbt * N(-f1)
                else:
                    term = ((s0 / smin)**(-w)) * N(-f1 + 2.0 *
                                                   b * sqrt(t) / v) - expbt * N(-f1)

                v = df * (k - smin) - s0 * dq * N(-f1) + \
                    smin * df * N(-f2) + s0 * df * u * term

            else:
                d1 = (log(s0 / k) + (b + v * v / 2) * t) / v / sqrt(t)
                d2 = d1 - v * sqrt(t)
                if s0 == k:
                    term = N(-d1 + 2.0 * b * sqrt(t) / v) - expbt * N(-d1)
                elif s0 > k and w < -100:
                    term = -expbt * N(-d1)
                else:
                    term = ((s0 / k)**(-w)) * N(-d1 + 2.0 *
                                                b * sqrt(t) / v) - expbt * N(-d1)

                v = k * df * N(-d2) - s0 * dq * N(-d1) + s0 * df * u * term

        else:
            raise FinError("Unknown lookback option type:" +
                           str(self._option_type))

        return v

###############################################################################

    def value_mc(self,
                 valuation_date: Date,
                 spot_fx_rate: float,  # FORDOM
                 domestic_curve: DiscountCurve,
                 foreign_curve: DiscountCurve,
                 volatility: float,
                 spot_fx_rateMinMax: float,
                 num_paths: int = 10000,
                 num_steps_per_year: int = 252,
                 seed: int = 4242):
        """ Value FX Fixed Lookback option using Monte Carlo. """

        t = (self._expiry_date - valuation_date) / gDaysInYear
        S0 = spot_fx_rate

        df = domestic_curve._df(t)
        rd = -np.log(df)/t

        dq = foreign_curve._df(t)
        rf = -np.log(dq)/t

        mu = rd - rf

        num_time_steps = int(t * num_steps_per_year)

        option_type = self._option_type
        k = self._optionStrike

        smin = 0.0
        smax = 0.0

        if self._option_type == OptionTypes.EUROPEAN_CALL:
            smax = spot_fx_rateMinMax
            if smax < S0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")
        elif self._option_type == OptionTypes.EUROPEAN_PUT:
            smin = spot_fx_rateMinMax
            if smin > S0:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")

        model = FinGBMProcess()
        Sall = model.get_paths(
            num_paths,
            num_time_steps,
            t,
            mu,
            S0,
            volatility,
            seed)

        # Due to antithetics we have doubled the number of paths
        num_paths = 2 * num_paths
        payoff = np.zeros(num_paths)

        if option_type == OptionTypes.EUROPEAN_CALL:
            SMax = np.max(Sall, axis=1)
            smaxs = np.ones(num_paths) * smax
            payoff = np.maximum(SMax - k, 0.0)
            payoff = np.maximum(payoff, smaxs - k)
        elif option_type == OptionTypes.EUROPEAN_PUT:
            SMin = np.min(Sall, axis=1)
            smins = np.ones(num_paths) * smin
            payoff = np.maximum(k - SMin, 0.0)
            payoff = np.maximum(payoff, k - smins)
        else:
            raise FinError("Unknown lookback option type:" + str(option_type))

        v = payoff.mean() * exp(-rd*t)
        return v
