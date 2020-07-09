##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...models.FinGBMProcess import FinGBMProcess
#from ...products.equities.FinOption import FinOption
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate

##########################################################################
# TODO
# Attempt control variate adjustment to monte carlo
# Adjust for finite sampling in Monte Carlo or Analytic
# TIGHTEN UP LIMIT FOR W FROM 100
##########################################################################

from enum import Enum


class FinFXFixedLookbackOptionTypes(Enum):
    FIXED_CALL = 1
    FIXED_PUT = 2

##########################################################################
# FIXED STRIKE LOOKBACK CALL PAYS MAX(SMAX-K,0)
# FIXED STRIKE LOOKBACK PUT PAYS MAX(K-SMIN,0)
##########################################################################


class FinFXFixedLookbackOption():
    ''' The Class for FX Fixed Strike Lookback options. '''

    def __init__(self,
                 expiryDate: FinDate,
                 optionType: FinFXFixedLookbackOptionTypes,
                 optionStrike: float):
        ''' Create option with expiry date, option type and the option strike
        '''

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._optionType = optionType
        self._optionStrike = optionStrike

##########################################################################

    def value(self,
              valueDate,
              stockPrice,
              domesticCurve,
              foreignCurve,
              volatility,
              stockMinMax):
        ''' Value FX Fixed Lookback Option using Black Scholes model and
        analytical formulae. '''

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = domesticCurve.df(t)
        r = -np.log(df)/t

        dq = foreignCurve.df(t)
        q = -np.log(dq)/t

        v = volatility
        s0 = stockPrice
        k = self._optionStrike
        smin = 0.0
        smax = 0.0

        if self._optionType == FinFXFixedLookbackOptionTypes.FIXED_CALL:
            smax = stockMinMax
            if smax < s0:
                raise FinError(
                    "The Smax value must be >= the stock price.")
        elif self._optionType == FinFXFixedLookbackOptionTypes.FIXED_PUT:
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
        if self._optionType == FinFXFixedLookbackOptionTypes.FIXED_CALL:

            if k > smax:
                d1 = (log(s0/k) + (b+v*v/2.0)*t)/v/sqrt(t)
                d2 = d1 - v * sqrt(t)

                if s0 == k:
                    term = -N(d1 - 2.0 * b * sqrt(t) / v) + expbt * N(d1)
                elif s0 < k and w > 100:
                    term = expbt * N(d1)
                else:
                    term = (-(s0 / k)**(-w)) * N(d1 - 2.0 *
                                                 b * sqrt(t) / v) + expbt*N(d1)

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
                        N(e1 - 2.0*b*sqrt(t) / v) + expbt * N(e1)

                v = df * (smax - k) + s0 * dq * N(e1) - \
                    smax * df * N(e2) + s0 * df * u * term

        elif self._optionType == FinFXFixedLookbackOptionTypes.FIXED_PUT:

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

    def valueMC(self,
                valueDate,
                spotFXRate,  # FORDOM
                domesticCurve,
                foreignCurve,
                volatility,
                spotFXRateMinMax,
                numPaths=10000,
                numStepsPerYear=252,
                seed=4242):
        ''' Value FX Fixed Lookback option using Monte Carlo. '''

        t = (self._expiryDate - valueDate) / gDaysInYear
        S0 = spotFXRate

        df = domesticCurve.df(t)
        rd = -np.log(df)/t

        dq = foreignCurve.df(t)
        rf = -np.log(dq)/t

        mu = rd - rf

        numTimeSteps = int(t * numStepsPerYear)

        optionType = self._optionType
        k = self._optionStrike

        smin = 0.0
        smax = 0.0

        if self._optionType == FinFXFixedLookbackOptionTypes.FIXED_CALL:
            smax = spotFXRateMinMax
            if smax < S0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")
        elif self._optionType == FinFXFixedLookbackOptionTypes.FIXED_PUT:
            smin = spotFXRateMinMax
            if smin > S0:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")

        model = FinGBMProcess()
        Sall = model.getPaths(
            numPaths,
            numTimeSteps,
            t,
            mu,
            S0,
            volatility,
            seed)

        # Due to antithetics we have doubled the number of paths
        numPaths = 2 * numPaths
        payoff = np.zeros(numPaths)

        if optionType == FinFXFixedLookbackOptionTypes.FIXED_CALL:
            SMax = np.max(Sall, axis=1)
            smaxs = np.ones(numPaths) * smax
            payoff = np.maximum(SMax - k, 0)
            payoff = np.maximum(payoff, smaxs - k)
        elif optionType == FinFXFixedLookbackOptionTypes.FIXED_PUT:
            SMin = np.min(Sall, axis=1)
            smins = np.ones(numPaths) * smin
            payoff = np.maximum(k - SMin, 0)
            payoff = np.maximum(payoff, k - smins)
        else:
            raise FinError("Unknown lookback option type:" + str(optionType))

        v = payoff.mean() * exp(-rd*t)
        return v
