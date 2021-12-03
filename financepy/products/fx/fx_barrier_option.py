##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np
from enum import Enum

from ...utils.error import FinError
from ...utils.math import N
from ...utils.global_vars import gDaysInYear
from ...products.fx.fx_option import FXOption
from ...models.process_simulator import FinProcessSimulator
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date


###############################################################################


class FinFXBarrierTypes(Enum):
    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8


###############################################################################


class FXBarrierOption(FXOption):

    def __init__(self,
                 expiry_date: Date,
                 strike_fx_rate: float,  # 1 unit of foreign in domestic
                 currency_pair: str,  # FORDOM
                 option_type: FinFXBarrierTypes,
                 barrier_level: float,
                 num_observations_per_year: int,
                 notional: float,
                 notional_currency: str):
        """ Create FX Barrier option product. This is an option that cancels if
        the FX rate crosses a barrier during the life of the option. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strike_fx_rate = float(strike_fx_rate)
        self._currency_pair = currency_pair
        self._barrier_level = float(barrier_level)
        self._num_observations_per_year = int(num_observations_per_year)
        self._option_type = option_type
        self._notional = notional
        self._notional_currency = notional_currency

    ##########################################################################

    def value(self,
              valuation_date,
              spot_fx_rate,
              dom_discount_curve,
              for_discount_curve,
              model):
        """ Value FX Barrier Option using Black-Scholes model with closed-form
        analytical models. """

        # This prices the option using the formulae given in the paper
        # by Clewlow, Llanos and Strickland December 1994 which can be found at
        # https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf

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

        K = self._strike_fx_rate
        S0 = spot_fx_rate
        h = self._barrier_level

        t = (self._expiry_date - valuation_date) / gDaysInYear
        lnS0k = log(float(S0) / K)
        sqrtT = sqrt(t)

        dq = for_discount_curve._df(t)
        df = dom_discount_curve._df(t)
        rd = -log(df) / t
        rf = -log(dq) / t

        volatility = model._volatility
        sigmaRootT = volatility * sqrtT
        v2 = volatility * volatility
        mu = rd - rf
        d1 = (lnS0k + (mu + v2 / 2.0) * t) / sigmaRootT
        d2 = (lnS0k + (mu - v2 / 2.0) * t) / sigmaRootT

        c = S0 * dq * N(d1) - K * df * N(d2)
        p = K * df * N(-d2) - S0 * dq * N(-d1)
        #        print("CALL:",c,"PUT:",p)

        if self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL and S0 <= h:
            return 0.0
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_CALL and S0 >= h:
            return 0.0
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_PUT and S0 >= h:
            return 0.0
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT and S0 <= h:
            return 0.0
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL and S0 <= h:
            return c
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_CALL and S0 >= h:
            return c
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_PUT and S0 >= h:
            return p
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT and S0 <= h:
            return p

        num_observations = t * self._num_observations_per_year

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        if self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        ll = (mu + v2 / 2.0) / v2
        y = log(h * h / (S0 * K)) / sigmaRootT + ll * sigmaRootT
        x1 = log(S0 / h) / sigmaRootT + ll * sigmaRootT
        y1 = log(h / S0) / sigmaRootT + ll * sigmaRootT
        hOverS = h / S0

        if self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            if h >= K:
                c_do = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * N(y1) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y1 - sigmaRootT)
                price = c_do
            else:
                c_di = S0 * dq * pow(hOverS, 2.0 * ll) * N(y) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y - sigmaRootT)
                price = c - c_di
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            if h <= K:
                c_di = S0 * dq * pow(hOverS, 2.0 * ll) * N(y) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y - sigmaRootT)
                price = c_di
            else:
                c_do = S0 * dq * N(x1) \
                    - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * N(y1) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y1 - sigmaRootT)
                price = c - c_do
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            if h >= K:
                c_ui = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c_ui
            else:
                price = c
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            if h > K:
                c_ui = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c - c_ui
            else:
                price = 0.0
        elif self._option_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            if h > K:
                p_ui = -S0 * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigmaRootT)
                price = p_ui
            else:
                p_uo = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                       - K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    N(-y1 + sigmaRootT)
                price = p - p_uo
        elif self._option_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            if h >= K:
                p_ui = -S0 * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigmaRootT)
                price = p - p_ui
            else:
                p_uo = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                       - K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    N(-y1 + sigmaRootT)
                price = p_uo
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            if h >= K:
                price = 0.0
            else:
                p_di = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                       - K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p - p_di
        elif self._option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            if h >= K:
                price = p
            else:
                p_di = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                       - K * df * pow(hOverS, 2.0 * ll - 2.0) * \
                    (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p_di
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        return price

    ###############################################################################

    def value_mc(self,
                 valuation_date,
                 spot_fx_rate,
                 dom_interest_rate,
                 process_type,
                 model_params,
                 num_ann_steps=552,
                 num_paths=5000,
                 seed=4242):
        """ Value the FX Barrier Option using Monte Carlo. """

        t = (self._expiry_date - valuation_date) / gDaysInYear
        num_time_steps = int(t * num_ann_steps)
        K = self._strike_fx_rate
        B = self._barrier_level
        S0 = spot_fx_rate
        option_type = self._option_type

        process = FinProcessSimulator()

        rd = dom_interest_rate

        #######################################################################

        if option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL and S0 <= B:
            return 0.0
        elif option_type == FinFXBarrierTypes.UP_AND_OUT_CALL and S0 >= B:
            return 0.0
        elif option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT and S0 <= B:
            return 0.0
        elif option_type == FinFXBarrierTypes.UP_AND_OUT_PUT and S0 >= B:
            return 0.0

        #######################################################################

        simple_call = False
        simple_put = False

        if option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL and S0 <= B:
            simple_call = True
        elif option_type == FinFXBarrierTypes.UP_AND_IN_CALL and S0 >= B:
            simple_call = True
        elif option_type == FinFXBarrierTypes.UP_AND_IN_PUT and S0 >= B:
            simple_put = True
        elif option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT and S0 <= B:
            simple_put = True

        if simple_put or simple_call:
            Sall = process.get_process(
                process_type, t, model_params, 1, num_paths, seed)

        if simple_call:
            sT = Sall[:, -1]
            c = (np.maximum(sT - K, 0.0)).mean()
            c = c * exp(-rd * t)
            return c

        if simple_put:
            sT = Sall[:, -1]
            p = (np.maximum(K - sT, 0.0)).mean()
            p = p * exp(-rd * t)
            return p

        # Get full set of paths
        Sall = process.get_process(process_type,
                                   t,
                                   model_params,
                                   num_time_steps,
                                   num_paths,
                                   seed)

        (num_paths, num_time_steps) = Sall.shape

        if option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL or \
                option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL or \
                option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT or \
                option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:

            barrierCrossedFromAbove = [False] * num_paths

            for p in range(0, num_paths):
                barrierCrossedFromAbove[p] = np.any(Sall[p] <= B)

        if option_type == FinFXBarrierTypes.UP_AND_IN_CALL or \
                option_type == FinFXBarrierTypes.UP_AND_OUT_CALL or \
                option_type == FinFXBarrierTypes.UP_AND_IN_PUT or \
                option_type == FinFXBarrierTypes.UP_AND_OUT_PUT:

            barrierCrossedFromBelow = [False] * num_paths
            for p in range(0, num_paths):
                barrierCrossedFromBelow[p] = np.any(Sall[p] >= B)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if option_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif option_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromAbove
        elif option_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromBelow
        elif option_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif option_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromBelow
        elif option_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif option_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif option_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromAbove
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        v = payoff.mean() * exp(-rd * t)

        return v

    ###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE FX RATE", self._strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self._currency_pair)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("BARRIER LEVEL", self._barrier_level)
        s += label_to_string("NUM OBSERVATIONS",
                             self._num_observations_per_year)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("NOTIONAL CURRENCY", self._notional_currency, "")
        return s

    ###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
