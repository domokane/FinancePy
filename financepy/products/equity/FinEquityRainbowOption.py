##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
import numpy as np
from typing import List

from ...finutils.FinMath import N, M
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...models.FinGBMProcess import FinGBMProcess
from ...products.equity.FinEquityOption import FinEquityOption
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate

from enum import Enum


class FinEquityRainbowOptionTypes(Enum):
    CALL_ON_MAXIMUM = 1
    PUT_ON_MAXIMUM = 2
    CALL_ON_MINIMUM = 3
    PUT_ON_MINIMUM = 4
    CALL_ON_NTH = 5  # MAX(NTH(S1,S2,...,SN)-K,0)
    PUT_ON_NTH = 6  # MAX(K-NTH(S1,S2,...,SN),0)

###############################################################################


def payoffValue(s, payoffTypeValue, payoffParams):

    if payoffTypeValue == FinEquityRainbowOptionTypes.CALL_ON_MINIMUM.value:
        k = payoffParams[0]
        payoff = np.maximum(np.min(s, axis=1) - k, 0.0)
    elif payoffTypeValue == FinEquityRainbowOptionTypes.CALL_ON_MAXIMUM.value:
        k = payoffParams[0]
        payoff = np.maximum(np.max(s, axis=1) - k, 0.0)
    elif payoffTypeValue == FinEquityRainbowOptionTypes.PUT_ON_MINIMUM.value:
        k = payoffParams[0]
        payoff = np.maximum(k - np.min(s, axis=1), 0.0)
    elif payoffTypeValue == FinEquityRainbowOptionTypes.PUT_ON_MAXIMUM.value:
        k = payoffParams[0]
        payoff = np.maximum(k - np.max(s, axis=1), 0.0)
    elif payoffTypeValue == FinEquityRainbowOptionTypes.CALL_ON_NTH.value:
        n = payoffParams[0]
        k = payoffParams[1]
        ssorted = np.sort(s)
        assetn = ssorted[:, -n]
        payoff = np.maximum(assetn - k, 0.0)
    elif payoffTypeValue == FinEquityRainbowOptionTypes.PUT_ON_NTH.value:
        n = payoffParams[0]
        k = payoffParams[1]
        ssorted = np.sort(s)
        assetn = ssorted[:, -n]
        payoff = np.maximum(k - assetn, 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff

###############################################################################


def valueMCFast(t,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                betas,
                numAssets,
                payoffType,
                payoffParams,
                numPaths=10000,
                seed=4242):

    np.random.seed(seed)

    df = discountCurve._df(t)
    r = -log(df)/t

    qs = []
    for curve in dividendCurves:
        dq = curve._df(t)
        q = -np.log(dq)/t
        qs.append(q)

    qs = np.array(qs)
    
    mus = r - qs

    model = FinGBMProcess()

    numTimeSteps = 2
    Sall = model.getPathsAssets(numAssets, numPaths, numTimeSteps,
                                t, mus, stockPrices, volatilities, betas, seed)

    payoff = payoffValue(Sall, payoffType.value, payoffParams)
    payoff = np.mean(payoff)
    v = payoff * exp(-r * t)
    return v

###############################################################################


class FinEquityRainbowOption(FinEquityOption):

    def __init__(self,
                 expiryDate: FinDate,
                 payoffType: FinEquityRainbowOptionTypes,
                 payoffParams: List[float],
                 numAssets: int):

        checkArgumentTypes(self.__init__, locals())

        self._validatePayoff(payoffType, payoffParams, numAssets)

        self._expiryDate = expiryDate
        self._payoffType = payoffType
        self._payoffParams = payoffParams
        self._numAssets = numAssets

###############################################################################

    def _validate(self,
                  stockPrices,
                  dividendCurves,
                  volatilities,
                  betas):

        if len(stockPrices) != self._numAssets:
            raise FinError(
                "Stock prices must be a vector of length "
                + str(self._numAssets))

        if len(dividendCurves) != self._numAssets:
            raise FinError(
                "Dividend curves must be a vector of length "
                + str(self._numAssets))

        if len(volatilities) != self._numAssets:
            raise FinError(
                "Volatilities must be a vector of length "
                + str(self._numAssets))

        if len(betas) != self._numAssets:
            raise FinError("Betas must be a vector of length " +
                           str(self._numAssets))

###############################################################################

    def _validatePayoff(self, payoffType, payoffParams, numAssets):

        numParams = 0

        if payoffType == FinEquityRainbowOptionTypes.CALL_ON_MINIMUM:
            numParams = 1
        elif payoffType == FinEquityRainbowOptionTypes.CALL_ON_MAXIMUM:
            numParams = 1
        elif payoffType == FinEquityRainbowOptionTypes.PUT_ON_MINIMUM:
            numParams = 1
        elif payoffType == FinEquityRainbowOptionTypes.PUT_ON_MAXIMUM:
            numParams = 1
        elif payoffType == FinEquityRainbowOptionTypes.CALL_ON_NTH:
            numParams = 2
        elif payoffType == FinEquityRainbowOptionTypes.PUT_ON_NTH:
            numParams = 2
        else:
            raise FinError("Unknown payoff type")

        if len(payoffParams) != numParams:
            raise FinError(
                "Number of parameters required for " +
                str(payoffType) +
                " must be " +
                str(numParams))

        if payoffType == FinEquityRainbowOptionTypes.CALL_ON_NTH \
           or payoffType == FinEquityRainbowOptionTypes.PUT_ON_NTH:
            n = payoffParams[0]
            if n < 1 or n > numAssets:
                raise FinError("Nth parameter must be 1 to " + str(numAssets))

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrices: np.ndarray,
              discountCurve: FinDiscountCurve,
              dividendCurves: (list),
              volatilities: np.ndarray,
              corrMatrix: np.ndarray):

        if self._numAssets != 2:
            raise FinError("Analytical results for two assets only.")

        if corrMatrix.ndim != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if corrMatrix.shape[0] != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if corrMatrix.shape[1] != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")


        # Use result by Stulz (1982) given by Haug Page 211
        t = (self._expiryDate - valueDate) / gDaysInYear
        r = discountCurve.zeroRate(self._expiryDate)

        q1 = dividendCurves[0].zeroRate(self._expiryDate)
        q2 = dividendCurves[1].zeroRate(self._expiryDate)

        dividendYields = [q1, q2]

        self._validate(stockPrices,
                       dividendYields,
                       volatilities,
                       corrMatrix)

#        q1 = dividendYields[0]
#        q2 = dividendYields[1]

        rho = corrMatrix[0][1]
        s1 = stockPrices[0]
        s2 = stockPrices[1]
        b1 = r - q1
        b2 = r - q2
        v1 = volatilities[0]
        v2 = volatilities[1]
        k = self._payoffParams[0]

        v = sqrt(v1 * v1 + v2 * v2 - 2 * rho * v1 * v2)
        d = (log(s1 / s2) + (b1 - b2 + v * v / 2) * t) / v / sqrt(t)
        y1 = (log(s1 / k) + (b1 + v1 * v1 / 2) * t) / v1 / sqrt(t)
        y2 = (log(s2 / k) + (b2 + v2 * v2 / 2) * t) / v2 / sqrt(t)
        rho1 = (v1 - rho * v2) / v
        rho2 = (v2 - rho * v1) / v
        dq1 = exp(-q1 * t)
        dq2 = exp(-q2 * t)
        df = exp(-r * t)

        if self._payoffType == FinEquityRainbowOptionTypes.CALL_ON_MAXIMUM:
            v = s1 * dq1 * M(y1, d, rho1) + s2 * dq2 * M(y2, -d + v * sqrt(t), rho2) \
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t), -y2 + v2 * sqrt(t), rho))
        elif self._payoffType == FinEquityRainbowOptionTypes.CALL_ON_MINIMUM:
            v = s1 * dq1 * M(y1, -d, -rho1) + s2 * dq2 * M(y2, d - v * sqrt(t), -rho2) \
                - k * df * M(y1 - v1 * sqrt(t), y2 - v2 * sqrt(t), rho)
        elif self._payoffType == FinEquityRainbowOptionTypes.PUT_ON_MAXIMUM:
            cmax1 = s2 * dq2 + s1 * dq1 * N(d) - s2 * dq2 * N(d - v * sqrt(t))
            cmax2 = s1 * dq1 * M(y1, d, rho1) + s2 * dq2 * M(y2, -d + v * sqrt(t), rho2) \
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t), -y2 + v2 * sqrt(t), rho))
            v = k * df - cmax1 + cmax2
        elif self._payoffType == FinEquityRainbowOptionTypes.PUT_ON_MINIMUM:
            cmin1 = s1 * dq1 - s1 * dq1 * N(d) + s2 * dq2 * N(d - v * sqrt(t))
            cmin2 = s1 * dq1 * M(y1, -d, -rho1) + s2 * dq2 * M(y2, d - v * sqrt(
                t), -rho2) - k * df * M(y1 - v1 * sqrt(t), y2 - v2 * sqrt(t), rho)
            v = k * df - cmin1 + cmin2
        else:
            raise FinError("Unsupported Rainbow option type")

        return v

###############################################################################

    def valueMC(self,
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths=10000,
                seed=4242):

        self._validate(stockPrices,
                       dividendCurves,
                       volatilities,
                       corrMatrix)

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        t = (self._expiryDate - valueDate) / gDaysInYear

        v = valueMCFast(t,
                        stockPrices,
                        discountCurve,
                        dividendCurves,
                        volatilities,
                        corrMatrix,
                        self._numAssets,
                        self._payoffType,
                        self._payoffParams,
                        numPaths,
                        seed)

        return v

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("PAYOFF TYPE", self._payoffType)
        s += labelToString("PAYOFF PARAMS", self._payoffParams)
        s += labelToString("NUM ASSETS TYPE", self._numAssets, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
