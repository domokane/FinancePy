##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...models.FinGBMProcess import FinGBMProcess
from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate

##########################################################################
# TODO
# Attempt control variate adjustment to monte carlo
# Adjust for finite sampling in Monte Carlo or Analytic
# TIGHTEN UP LIMIT FOR W FROM 100
##########################################################################

from enum import Enum


class FinEquityFixedLookbackOptionTypes(Enum):
    FIXED_CALL = 1
    FIXED_PUT = 2

##########################################################################
# FIXED STRIKE LOOKBACK CALL PAYS MAX(SMAX-K,0)
# FIXED STRIKE LOOKBACK PUT PAYS MAX(K-SMIN,0)
##########################################################################


class FinEquityFixedLookbackOption(FinEquityOption):

    def __init__(self,
                 expiryDate: FinDate,
                 optionType: FinEquityFixedLookbackOptionTypes,
                 optionStrike: float):

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._optionType = optionType
        self._optionStrike = optionStrike

##########################################################################

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              volatility,
              stockMinMax):

        t = (self._expiryDate - valueDate) / gDaysInYear
        r = discountCurve.zeroRate(self._expiryDate)
        
        v = volatility
        s0 = stockPrice
        q = dividendYield
        k = self._optionStrike
        smin = 0.0
        smax = 0.0

        if self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_CALL:
            smax = stockMinMax
            if smax < s0:
                raise FinError(
                    "The Smax value must be >= the stock price.")
        elif self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_PUT:
            smin = stockMinMax
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
        if self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_CALL:

            if k > smax:
                d1 = (log(s0 / k) + (b + v * v / 2) * t) / v / sqrt(t)
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
                e1 = (log(s0 / smax) + (r - q + v * v / 2) * t) / v / sqrt(t)
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

        elif self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_PUT:

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
                           str(self._optionType))

        return v

###############################################################################

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
        r = discountCurve.zeroRate(self._expiryDate)
        mu = r - dividendYield

        numTimeSteps = int(t * numStepsPerYear)

        optionType = self._optionType
        k = self._optionStrike

        smin = 0.0
        smax = 0.0

        if self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_CALL:
            smax = stockMinMax
            if smax < stockPrice:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")
        elif self._optionType == FinEquityFixedLookbackOptionTypes.FIXED_PUT:
            smin = stockMinMax
            if smin > stockPrice:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")

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

        if optionType == FinEquityFixedLookbackOptionTypes.FIXED_CALL:
            SMax = np.max(Sall, axis=1)
            smaxs = np.ones(numPaths) * smax
            payoff = np.maximum(SMax - k, 0)
            payoff = np.maximum(payoff, smaxs - k)
        elif optionType == FinEquityFixedLookbackOptionTypes.FIXED_PUT:
            SMin = np.min(Sall, axis=1)
            smins = np.ones(numPaths) * smin
            payoff = np.maximum(k - SMin, 0)
            payoff = np.maximum(payoff, k - smins)
        else:
            raise FinError("Unknown lookback option type:" + str(optionType))

        v = payoff.mean() * exp(-r * t)
        return v
