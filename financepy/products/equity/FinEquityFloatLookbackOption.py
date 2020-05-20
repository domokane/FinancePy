##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from math import exp, log, sqrt
from enum import Enum

from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...models.FinGBMProcess import FinGBMProcess
from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinHelperFunctions import labelToString

##########################################################################
# TODO
# Attempt control variate adjustment to monte carlo
# Adjust for finite sampling in Monte Carlo or Analytic
# TIGHTEN UP LIMIT FOR W FROM 100
##########################################################################


class FinEquityFloatLookbackOptionTypes(Enum):
    FLOATING_CALL = 1
    FLOATING_PUT = 2

##########################################################################
# FLOAT STRIKE LOOKBACK CALL PAYS MAX(S(T)-SMIN,0)
# FLOAT STRIKE LOOKBACK PUT PAYS MAX(SMAX-S(T),0)
##########################################################################


class FinEquityFloatLookbackOption(FinEquityOption):

    def __init__(self,
                 expiryDate,
                 optionType):

        self._expiryDate = expiryDate
        self._optionType = optionType

##########################################################################

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              volatility,
              stockMinMax):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(t)
        r = -np.log(df)/t

        v = volatility
        s0 = stockPrice
        q = dividendYield
        smin = 0.0
        smax = 0.0

        if self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_CALL:
            smin = stockMinMax
            if smin > s0:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_PUT:
            smax = stockMinMax
            if smax < s0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        if abs(r - q) < gSmall:
            q = r + gSmall

        dq = exp(-q * t)
        df = exp(-r * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = exp(b * t)

        # Taken from Haug Page 142
        if self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_CALL:

            a1 = (log(s0 / smin) + (b + (v**2) / 2.0) * t) / v / sqrt(t)
            a2 = a1 - v * sqrt(t)

            if smin == s0:
                term = N(-a1 + 2.0 * b * sqrt(t) / v) - expbt * N(-a1)
            elif s0 < smin and w < -100:
                term = - expbt * N(-a1)
            else:
                term = ((s0 / smin)**(-w)) * N(-a1 + 2.0 *
                                               b * sqrt(t) / v) - expbt * N(-a1)

            v = s0 * dq * N(a1) - smin * df * N(a2) + s0 * df * u * term

        elif self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_PUT:

            b1 = (log(s0 / smax) + (b + (v**2) / 2.0) * t) / v / sqrt(t)
            b2 = b1 - v * sqrt(t)

            if smax == s0:
                term = -N(b1 - 2.0 * b * sqrt(t) / v) + expbt * N(b1)
            elif s0 < smax and w > 100:
                term = expbt * N(b1)
            else:
                term = (-(s0 / smax)**(-w)) * \
                    N(b1 - 2.0 * b * sqrt(t) / v) + expbt * N(b1)

            v = smax * df * N(-b2) - s0 * dq * N(-b1) + s0 * df * u * term

        else:
            raise FinError("Unknown lookback option type:" +
                           str(self._optionType))

        return v

##########################################################################

    def valueMC(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            volatility,
            stockMinMax,
            numPaths=10000,
            numStepsPerYear=252,
            seed=4242):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(t)
        r = -np.log(df)/t

        numTimeSteps = int(t * numStepsPerYear)
        mu = r - dividendYield

        optionType = self._optionType
        smin = 0.0
        smax = 0.0

        if self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_CALL:
            smin = stockMinMax
            if smin > stockPrice:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._optionType == FinEquityFloatLookbackOptionTypes.FLOATING_PUT:
            smax = stockMinMax
            if smax < stockPrice:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        model = FinGBMProcess()
        Sall = model.getPaths(
            numPaths,
            numTimeSteps,
            t,
            mu,
            stockPrice,
            volatility,
            seed)

        # Due to antithetics we have doubled the number of paths
        numPaths = 2 * numPaths
        payoff = np.zeros(numPaths)

        if optionType == FinEquityFloatLookbackOptionTypes.FLOATING_CALL:
            SMin = np.min(Sall, axis=1)
            SMin = np.minimum(SMin, smin)
            payoff = np.maximum(Sall[:, -1] - SMin, 0)
        elif optionType == FinEquityFloatLookbackOptionTypes.FLOATING_PUT:
            SMax = np.max(Sall, axis=1)
            SMax = np.maximum(SMax, smax)
            payoff = np.maximum(SMax - Sall[:, -1], 0)
        else:
            raise FinError("Unknown lookback option type:" + str(optionType))

        v = payoff.mean() * exp(-r * t)
        return v

##########################################################################
