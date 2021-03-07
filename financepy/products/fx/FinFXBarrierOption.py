##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np
from enum import Enum


from ...utils.FinError import FinError
from ...utils.math import N
from ...utils.global_vars import gDaysInYear
from ...products.fx.FinFXOption import FinFXOption
from ...models.process_simulator import FinProcessSimulator
from ...utils.helpers import labelToString, check_argument_types
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


class FinFXBarrierOption(FinFXOption):

    def __init__(self,
                 expiry_date: Date,
                 strikeFXRate: float,  # 1 unit of foreign in domestic
                 currencyPair: str,  # FORDOM
                 optionType: FinFXBarrierTypes,
                 barrierLevel: float,
                 numObservationsPerYear: int,
                 notional: float,
                 notionalCurrency: str):
        """ Create FX Barrier option product. This is an option that cancels if
        the FX rate crosses a barrier during the life of the option. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strikeFXRate = float(strikeFXRate)
        self._currencyPair = currencyPair
        self._barrierLevel = float(barrierLevel)
        self._numObservationsPerYear = int(numObservationsPerYear)
        self._optionType = optionType
        self._notional = notional
        self._notionalCurrency = notionalCurrency

##########################################################################

    def value(self,
              valuation_date,
              spotFXRate,
              domDiscountCurve,
              forDiscountCurve,
              model):
        """ Value FX Barrier Option using Black-Scholes model with closed-form
        analytical models. """

        # This prices the option using the formulae given in the paper
        # by Clewlow, Llanos and Strickland December 1994 which can be found at
        # https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf

        K = self._strikeFXRate
        S0 = spotFXRate
        h = self._barrierLevel

        t = (self._expiry_date - valuation_date) / gDaysInYear
        lnS0k = log(float(S0)/K)
        sqrtT = sqrt(t)

        dq = forDiscountCurve._df(t)
        df = domDiscountCurve._df(t)
        rd = -log(df)/t
        rf = -log(dq)/t

        volatility = model._volatility
        sigmaRootT = volatility * sqrtT
        v2 = volatility * volatility
        mu = rd - rf
        d1 = (lnS0k + (mu + v2 / 2.0) * t) / sigmaRootT
        d2 = (lnS0k + (mu - v2 / 2.0) * t) / sigmaRootT

        c = S0 * dq * N(d1) - K * df * N(d2)
        p = K * df * N(-d2) - S0 * dq * N(-d1)
#        print("CALL:",c,"PUT:",p)

        if self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL and S0 <= h:
            return 0.0
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_CALL and S0 >= h:
            return 0.0
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_PUT and S0 >= h:
            return 0.0
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT and S0 <= h:
            return 0.0
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL and S0 <= h:
            return c
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_CALL and S0 >= h:
            return c
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_PUT and S0 >= h:
            return p
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT and S0 <= h:
            return p

        numObservations = t * self._numObservationsPerYear

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        if self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        ll = (mu + v2 / 2.0) / v2
        y = log(h * h / (S0 * K)) / sigmaRootT + ll * sigmaRootT
        x1 = log(S0 / h) / sigmaRootT + ll * sigmaRootT
        y1 = log(h / S0) / sigmaRootT + ll * sigmaRootT
        hOverS = h / S0

        if self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            if h >= K:
                c_do = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * N(y1) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y1 - sigmaRootT)
                price = c_do
            else:
                c_di = S0 * dq * pow(hOverS, 2.0 * ll) * N(y) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * N(y - sigmaRootT)
                price = c - c_di
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL:
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
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_CALL:
            if h >= K:
                c_ui = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                    - S0 * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) *(N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c_ui
            else:
                price = c
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_CALL:
            if h > K:
                c_ui = S0 * dq * N(x1) - K * df * N(x1 - sigmaRootT) \
                     - S0 * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                     + K * df * pow(hOverS, 2.0 * ll - 2.0) * (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c - c_ui
            else:
                price = 0.0
        elif self._optionType == FinFXBarrierTypes.UP_AND_IN_PUT:
            if h > K:
                p_ui = -S0 * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigmaRootT)
                price = p_ui
            else:
                p_uo = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y1 + sigmaRootT)
                price = p - p_uo
        elif self._optionType == FinFXBarrierTypes.UP_AND_OUT_PUT:
            if h >= K:
                p_ui = -S0 * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                    + K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigmaRootT)
                price = p - p_ui
            else:
                p_uo = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y1 + sigmaRootT)
                price = p_uo
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            if h >= K:
                price = 0.0
            else:
                p_di = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p - p_di
        elif self._optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            if h >= K:
                price = p
            else:
                p_di = -S0 * dq * N(-x1) \
                    + K * df * N(-x1 + sigmaRootT) \
                    + S0 * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                    - K * df * pow(hOverS, 2.0 * ll - 2.0) * (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p_di
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        return price

###############################################################################

    def valueMC(self,
                valuation_date,
                spotFXRate,
                domInterestRate,
                processType,
                modelParams,
                numAnnSteps=552,
                num_paths=5000,
                seed=4242):
        """ Value the FX Barrier Option using Monte Carlo. """

        t = (self._expiry_date - valuation_date) / gDaysInYear
        numTimeSteps = int(t * numAnnSteps)
        K = self._strikeFXRate
        B = self._barrierLevel
        S0 = spotFXRate
        optionType = self._optionType

        process = FinProcessSimulator()

        rd = domInterestRate

        #######################################################################

        if optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL and S0 <= B:
            return 0.0
        elif optionType == FinFXBarrierTypes.UP_AND_OUT_CALL and S0 >= B:
            return 0.0
        elif optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT and S0 <= B:
            return 0.0
        elif optionType == FinFXBarrierTypes.UP_AND_OUT_PUT and S0 >= B:
            return 0.0

        #######################################################################

        simpleCall = False
        simplePut = False

        if optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL and S0 <= B:
            simpleCall = True
        elif optionType == FinFXBarrierTypes.UP_AND_IN_CALL and S0 >= B:
            simpleCall = True
        elif optionType == FinFXBarrierTypes.UP_AND_IN_PUT and S0 >= B:
            simplePut = True
        elif optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT and S0 <= B:
            simplePut = True

        if simplePut or simpleCall:
            Sall = process.getProcess(
                processType, t, modelParams, 1, num_paths, seed)

        if simpleCall:
            sT = Sall[:, -1]
            c = (np.maximum(sT - K, 0.0)).mean()
            c = c * exp(-rd * t)
            return c

        if simplePut:
            sT = Sall[:, -1]
            p = (np.maximum(K - sT, 0.0)).mean()
            p = p * exp(-rd * t)
            return p

        # Get full set of paths
        Sall = process.getProcess(processType,
                                  t,
                                  modelParams,
                                  numTimeSteps,
                                  num_paths,
                                  seed)

        (num_paths, numTimeSteps) = Sall.shape

        if optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL or \
           optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL or \
           optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT or \
           optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT:

            barrierCrossedFromAbove = [False] * num_paths

            for p in range(0, num_paths):
                barrierCrossedFromAbove[p] = np.any(Sall[p] <= B)

        if optionType == FinFXBarrierTypes.UP_AND_IN_CALL or \
           optionType == FinFXBarrierTypes.UP_AND_OUT_CALL or \
           optionType == FinFXBarrierTypes.UP_AND_IN_PUT or \
           optionType == FinFXBarrierTypes.UP_AND_OUT_PUT:

            barrierCrossedFromBelow = [False] * num_paths
            for p in range(0, num_paths):
                barrierCrossedFromBelow[p] = np.any(Sall[p] >= B)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if optionType == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromAbove
        elif optionType == FinFXBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromBelow
        elif optionType == FinFXBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinFXBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromBelow
        elif optionType == FinFXBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromAbove
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        v = payoff.mean() * exp(-rd*t)

        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("STRIKE FX RATE", self._strikeFXRate)
        s += labelToString("CURRENCY PAIR", self._currencyPair)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("BARRIER LEVEL", self._barrierLevel)
        s += labelToString("NUM OBSERVATIONS", self._numObservationsPerYear)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("NOTIONAL CURRENCY", self._notionalCurrency, "")
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)


###############################################################################

