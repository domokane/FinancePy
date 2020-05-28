##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
import numpy as np
from enum import Enum

from ...finutils.FinError import FinError
from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear
from ...products.equity.FinEquityOption import FinEquityOption
from ...models.FinProcessSimulator import FinProcessSimulator

from ...finutils.FinHelperFunctions import labelToString

class FinEquityBarrierTypes(Enum):
    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8


###############################################################################
##########################################################################

class FinEquityBarrierOption(FinEquityOption):
    ''' Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum FinEquityBarrierTypes. '''

    def __init__(self,
                 expiryDate,
                 strikePrice,
                 optionType,
                 barrierLevel,
                 numObservationsPerYear,
                 notional=1.0):

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._barrierLevel = float(barrierLevel)
        self._numObservationsPerYear = int(numObservationsPerYear)

        if optionType not in FinEquityBarrierTypes:
            raise FinError("Option Type ", optionType, " unknown.")

        self._optionType = optionType
        self._notional = notional

##########################################################################

    def value(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):

        # This prices the option using the formulae given in the paper
        # by Clewlow, Llanos and Strickland December 1994 which can be found at
        # https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf

        t = (self._expiryDate - valueDate) / gDaysInYear
        lnS0k = log(float(stockPrice) / self._strikePrice)
        sqrtT = sqrt(t)

        df = discountCurve.df(t)
        r = -np.log(df)/t

        k = self._strikePrice
        s = stockPrice
        h = self._barrierLevel

        volatility = model._volatility
        sigmaRootT = volatility * sqrtT
        v2 = volatility * volatility
        mu = r - dividendYield
        d1 = (lnS0k + (mu + v2 / 2.0) * t) / sigmaRootT
        d2 = (lnS0k + (mu - v2 / 2.0) * t) / sigmaRootT
        df = exp(-r * t)
        dq = exp(-dividendYield * t)

        c = s * dq * N(d1) - k * df * N(d2)
        p = k * df * N(-d2) - s * dq * N(-d1)
#        print("CALL:",c,"PUT:",p)

        if self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL and s <= h:
            return 0.0
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL and s >= h:
            return 0.0
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT and s >= h:
            return 0.0
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT and s <= h:
            return 0.0
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL and s <= h:
            return c
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_CALL and s >= h:
            return c
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_PUT and s >= h:
            return p
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT and s <= h:
            return p

        numObservations = t * self._numObservationsPerYear

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        if self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / numObservations))
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        l = (mu + v2 / 2.0) / v2
        y = log(h * h / (s * k)) / sigmaRootT + l * sigmaRootT
        x1 = log(s / h) / sigmaRootT + l * sigmaRootT
        y1 = log(h / s) / sigmaRootT + l * sigmaRootT
        hOverS = h / s

        if self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL:
            if h >= k:
                c_do = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigmaRootT)
                price = c_do
            else:
                c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigmaRootT)
                price = c - c_di
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL:
            if h <= k:
                c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigmaRootT)
                price = c_di
            else:
                c_do = s * dq * N(x1) \
                    - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigmaRootT)
                price = c - c_do
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_CALL:
            if h >= k:
                c_ui = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c_ui
            else:
                price = c
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL:
            if h > k:
                c_ui = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                     - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                     + k * df * pow(hOverS, 2.0 * l - 2.0) * (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c - c_ui
            else:
                price = 0.0
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_PUT:
            if h > k:
                p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigmaRootT)
                price = p_ui
            else:
                p_uo = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y1 + sigmaRootT)
                price = p - p_uo
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:
            if h >= k:
                p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigmaRootT)
                price = p - p_ui
            else:
                p_uo = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y1 + sigmaRootT)
                price = p_uo
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:
            if h >= k:
                price = 0.0
            else:
                p_di = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p - p_di
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT:
            if h >= k:
                price = p
            else:
                p_di = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p_di
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        v = price * self._notional
        return v

###############################################################################

    def valueMC(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            processType,
            modelParams,
            numAnnSteps=252,
            numPaths=10000,
            seed=4242):

        t = (self._expiryDate - valueDate) / gDaysInYear
        numTimeSteps = int(t * numAnnSteps)
        K = self._strikePrice
        B = self._barrierLevel
        optionType = self._optionType

        process = FinProcessSimulator()

        df = discountCurve.df(t)
        r = -np.log(df)/t

        #######################################################################

        if optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL and stockPrice <= B:
            return 0.0
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL and stockPrice >= B:
            return 0.0
        elif optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT and stockPrice <= B:
            return 0.0
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT and stockPrice >= B:
            return 0.0

        #######################################################################

        simpleCall = False
        simplePut = False

        if optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL and stockPrice <= B:
            simpleCall = True
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_CALL and stockPrice >= B:
            simpleCall = True
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_PUT and stockPrice >= B:
            simplePut = True
        elif optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT and stockPrice <= B:
            simplePut = True

        if simplePut or simpleCall:
            Sall = process.getProcess(
                processType, t, modelParams, 1, numPaths, seed)

        if simpleCall:
            c = (np.maximum(Sall[:, -1] - K, 0)).mean()
            c = c * exp(-r * t)
            return c

        if simplePut:
            p = (np.maximum(K - Sall[:, -1], 0)).mean()
            p = p * exp(-r * t)
            return p

        # Get full set of paths
        Sall = process.getProcess(
            processType,
            t,
            modelParams,
            numTimeSteps,
            numPaths,
            seed)
        (numPaths, numTimeSteps) = Sall.shape

        if optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL or \
           optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL or \
           optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT or \
           optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:

            barrierCrossedFromAbove = [False] * numPaths

            for p in range(0, numPaths):
                barrierCrossedFromAbove[p] = np.any(Sall[p] <= B)

        if optionType == FinEquityBarrierTypes.UP_AND_IN_CALL or \
           optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL or \
           optionType == FinEquityBarrierTypes.UP_AND_IN_PUT or \
           optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:

            barrierCrossedFromBelow = [False] * numPaths
            for p in range(0, numPaths):
                barrierCrossedFromBelow[p] = np.any(Sall[p] >= B)

        payoff = np.zeros(numPaths)
        ones = np.ones(numPaths)

        if optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0) * barrierCrossedFromAbove
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0) * barrierCrossedFromBelow
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0) * barrierCrossedFromBelow
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0) * barrierCrossedFromAbove
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        v = payoff.mean() * exp(- r * t)

        return v * self._notional

##########################################################################
