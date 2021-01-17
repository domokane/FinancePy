###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from enum import Enum

from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear
from ...products.equity.FinEquityOption import FinEquityOption
from ...models.FinProcessSimulator import FinProcessSimulator
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate


from ...finutils.FinMath import N

# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

###############################################################################


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


class FinEquityBarrierOption(FinEquityOption):
    ''' Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum FinEquityBarrierTypes. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikePrice: float,
                 optionType: FinEquityBarrierTypes,
                 barrierLevel: float,
                 numObservationsPerYear: (int, float) = 252,
                 notional: float = 1.0):
        ''' Create the FinEquityBarrierOption by specifying the expiry date,
        strike price, option type, barrier level, the number of observations
        per year and the notional. '''

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._barrierLevel = float(barrierLevel)
        self._numObservationsPerYear = int(numObservationsPerYear)

        if optionType not in FinEquityBarrierTypes:
            raise FinError("Option Type " + str(optionType) + " unknown.")

        self._optionType = optionType
        self._notional = notional

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: (float, np.ndarray),
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              model):
        ''' This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        '''

        if isinstance(stockPrice, int):
            stockPrice = float(stockPrice)

        if isinstance(stockPrice, float):
            stockPrices = [stockPrice]
        else:
            stockPrices = stockPrice

        values = []
        for s in stockPrices:
            v = self._valueOne(valueDate, s, discountCurve,
                               dividendCurve, model)
            values.append(v)

        if isinstance(stockPrice, float):
            return values[0]
        else:
            return np.array(values)

###############################################################################

    def _valueOne(self,
                  valueDate: FinDate,
                  stockPrice: (float, np.ndarray),
                  discountCurve: FinDiscountCurve,
                  dividendCurve: FinDiscountCurve,
                  model):
        ''' This values a single option. Because of its structure it cannot
        easily be vectorised which is why it has been wrapped. '''

        texp = (self._expiryDate - valueDate) / gDaysInYear

        if texp < 0:
            raise FinError("Option expires before value date.")

        texp = max(texp, 1e-6)

        lnS0k = np.log(stockPrice / self._strikePrice)
        sqrtT = np.sqrt(texp)

        r = discountCurve.ccRate(self._expiryDate)
        q = dividendCurve.ccRate(self._expiryDate)

        k = self._strikePrice
        s = stockPrice
        h = self._barrierLevel

        volatility = model._volatility
        sigmaRootT = volatility * sqrtT
        v2 = volatility * volatility
        mu = r - q
        d1 = (lnS0k + (mu + v2 / 2.0) * texp) / sigmaRootT
        d2 = (lnS0k + (mu - v2 / 2.0) * texp) / sigmaRootT
        df = np.exp(-r * texp)
        dq = np.exp(-q * texp)

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

        numObservations = 1 + texp * self._numObservationsPerYear

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        t = texp / numObservations

        if self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        l = (mu + v2 / 2.0) / v2
        y = np.log(h * h / (s * k)) / sigmaRootT + l * sigmaRootT
        x1 = np.log(s / h) / sigmaRootT + l * sigmaRootT
        y1 = np.log(h / s) / sigmaRootT + l * sigmaRootT
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

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendCurve: FinDiscountCurve,
                processType,
                modelParams,
                numAnnObs: int = 252,
                numPaths: int = 10000,
                seed: int = 4242):
        ''' A Monte-Carlo based valuation of the barrier option which simulates
        the evolution of the stock price of at a specified number of annual
        observation times until expiry to examine if the barrier has been
        crossed and the corresponding value of the final payoff, if any. It
        assumes a GBM model for the stock price. '''

        texp = (self._expiryDate - valueDate) / gDaysInYear
        numTimeSteps = int(texp * numAnnObs)
        K = self._strikePrice
        B = self._barrierLevel
        optionType = self._optionType

        process = FinProcessSimulator()

        r = discountCurve.zeroRate(self._expiryDate)
        
        # TODO - NEED TO DECIDE IF THIS IS PART OF MODEL PARAMS OR NOT ??????????????

        r = discountCurve.ccRate(self._expiryDate)
        q = dividendCurve.ccRate(self._expiryDate)
        
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
                processType, texp, modelParams, 1, numPaths, seed)

        if simpleCall:
            c = (np.maximum(Sall[:, -1] - K, 0.0)).mean()
            c = c * np.exp(-r * texp)
            return c

        if simplePut:
            p = (np.maximum(K - Sall[:, -1], 0.0)).mean()
            p = p * np.exp(-r * texp)
            return p

        # Get full set of paths
        Sall = process.getProcess(processType, texp, modelParams, numTimeSteps,
                                  numPaths, seed)

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
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromAbove
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromBelow
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinEquityBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromBelow
        elif optionType == FinEquityBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif optionType == FinEquityBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromAbove
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._optionType))

        v = payoff.mean() * np.exp(- r * texp)

        return v * self._notional

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("BARRIER LEVEL", self._barrierLevel)
        s += labelToString("NUM OBSERVATIONS", self._numObservationsPerYear)
        s += labelToString("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
