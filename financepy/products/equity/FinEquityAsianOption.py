##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from math import exp, log, sqrt

from numba import njit

from typing import Union

from ...finutils.FinMath import N, covar
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError

from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinOptionTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes 
from ...finutils.FinDate import FinDate

###############################################################################
# An Asian option on an arithmetic average and strike K has a payoff
# Max(SA(T)-K,0) where SA is the arithmetic average
# We define three dates
# - Valuation date for which we want the price
# - Start Averaging Date for when the averaging starts
# - Expiry date for when the payoff is made and the option expires
#
# In the model we have
# tv = is the time now
# t0 = time to the start averaging date in years
# t = time to the expiry date in years
# tau = length of averaging period in years at the start of the option
#
# We can be before the start of the averaging period in which case t0 > 0
# We can be after the start of the averaging period in which case we set t0=0
# and we note that t <= tau
#
# If we are in the averaging period then we need to know the accrued average
# I call this AA and the new average is now given by the accrued average
# The option payoff is now Max( (AA x (tau-t) + SA(t0) x (t-t0))/tau - K,0)
# This simplifies to
#
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#
###############################################################################


@njit(cache=True, fastmath=True)
def valueMC_NUMBA(t0, t, tau, K, n, optionType,
                  stockPrice,
                  interestRate,
                  dividendYield,
                  volatility,
                  numPaths,
                  seed,
                  accruedAverage):

    # Start pricing here
    np.random.seed(seed)

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(
            "In the averaging period you need to enter the accrued average")

        # we adjust the strike to account for the accrued coupon
        K = (K * tau + accruedAverage * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    mu = interestRate - dividendYield
    v2 = volatility**2
    dt = (t - t0) / n

    payoff_a = 0.0

    for p in range(0, numPaths):

        # evolve stock price to start of averaging period
        g = np.random.standard_normal(1)
        s_1 = stockPrice * exp((mu - v2 / 2.0) * t0 +
                               g[0] * sqrt(t0) * volatility)
        s_2 = stockPrice * exp((mu - v2 / 2.0) * t0 -
                               g[0] * sqrt(t0) * volatility)

        # enter averaging period
        s_1_arithmetic = 0.0
        s_2_arithmetic = 0.0

        g = np.random.standard_normal(n)

        for obs in range(0, n):

            s_1 = s_1 * exp((mu - v2 / 2.0) * dt +
                            g[obs] * sqrt(dt) * volatility)
            s_2 = s_2 * exp((mu - v2 / 2.0) * dt -
                            g[obs] * sqrt(dt) * volatility)

            s_1_arithmetic += s_1
            s_2_arithmetic += s_2

        s_1_arithmetic /= n
        s_2_arithmetic /= n

        if optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff_a += max(s_1_arithmetic - K, 0)
            payoff_a += max(s_2_arithmetic - K, 0)
        elif optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff_a += max(K - s_1_arithmetic, 0)
            payoff_a += max(K - s_2_arithmetic, 0)
        else:
            return None

    v_a = payoff_a * exp(-interestRate * t) / numPaths / 2.0
    v_a = v_a * multiplier
    return v_a

###############################################################################


@njit(cache=True, fastmath=True)
def valueMC_fast_NUMBA(t0, t, tau, K, n, optionType,
                       stockPrice,
                       interestRate,
                       dividendYield,
                       volatility,
                       numPaths,
                       seed,
                       accruedAverage):

    np.random.seed(seed)

    mu = interestRate - dividendYield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interestRate

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(
                "In the averaging period you need to enter the accrued average")

        # we adjust the strike to account for the accrued coupon
        K = (K * tau + accruedAverage * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

#    payoff_a = 0.0

    # evolve stock price to start of averaging period
    g = np.random.normal(0.0, 1.0, size=(numPaths))

    s_1 = np.empty(numPaths)
    s_2 = np.empty(numPaths)

    for ip in range(0, numPaths):
        s_1[ip] = stockPrice * \
            exp((mu - v2 / 2.0) * t0 + g[ip] * sqrt(t0) * volatility)
        s_2[ip] = stockPrice * \
            exp((mu - v2 / 2.0) * t0 - g[ip] * sqrt(t0) * volatility)

    s_1_arithmetic = np.zeros(numPaths)
    s_2_arithmetic = np.zeros(numPaths)

    for obs in range(0, n):

        g = np.random.normal(0.0, 1.0, size=(numPaths))

        for ip in range(0, numPaths):
            s_1[ip] = s_1[ip] * exp((mu - v2 / 2.0) *
                                    dt + g[ip] * sqrt(dt) * volatility)
            s_2[ip] = s_2[ip] * exp((mu - v2 / 2.0) *
                                    dt - g[ip] * sqrt(dt) * volatility)

        for ip in range(0, numPaths):
            s_1_arithmetic[ip] += s_1[ip] / n
            s_2_arithmetic[ip] += s_2[ip] / n

    if optionType == FinOptionTypes.EUROPEAN_CALL:
        payoff_a_1 = np.maximum(s_1_arithmetic - K, 0)
        payoff_a_2 = np.maximum(s_2_arithmetic - K, 0)
    elif optionType == FinOptionTypes.EUROPEAN_PUT:
        payoff_a_1 = np.maximum(K - s_1_arithmetic, 0)
        payoff_a_2 = np.maximum(K - s_2_arithmetic, 0)
    else:
        raise FinError("Unknown option type.")

    payoff_a = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    v_a = multiplier * payoff_a * exp(- r * t) / 2.0
    return v_a

##########################################################################


@njit(cache=True, fastmath=True)
def valueMC_fast_CV_NUMBA(t0, t, tau, K, n, optionType,
                          stockPrice,
                          interestRate,
                          dividendYield,
                          volatility,
                          numPaths,
                          seed,
                          accruedAverage,
                          v_g_exact):

    np.random.seed(seed)

    mu = interestRate - dividendYield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interestRate

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(
                "In averaging period you need to enter the accrued average")

        # we adjust the strike to account for the accrued coupon
        K = (K * tau + accruedAverage * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    g = np.random.normal(0.0, 1.0, size=(numPaths))

    s_1 = np.empty(numPaths)
    s_2 = np.empty(numPaths)

    for ip in range(0, numPaths):
        s_1[ip] = stockPrice * \
            exp((mu - v2 / 2.0) * t0 + g[ip] * sqrt(t0) * volatility)
        s_2[ip] = stockPrice * \
            exp((mu - v2 / 2.0) * t0 - g[ip] * sqrt(t0) * volatility)

    s_1_arithmetic = np.zeros(numPaths)
    s_2_arithmetic = np.zeros(numPaths)
    ln_s_1_geometric = np.zeros(numPaths)
    ln_s_2_geometric = np.zeros(numPaths)

    for obs in range(0, n):

        g = np.random.normal(0.0, 1.0, size=(numPaths))
        for ip in range(0, numPaths):
            s_1[ip] = s_1[ip] * exp((mu - v2 / 2.0) *
                                    dt + g[ip] * sqrt(dt) * volatility)
            s_2[ip] = s_2[ip] * exp((mu - v2 / 2.0) *
                                    dt - g[ip] * sqrt(dt) * volatility)

        for ip in range(0, numPaths):
            s_1_arithmetic[ip] += s_1[ip]
            s_2_arithmetic[ip] += s_2[ip]
            ln_s_1_geometric[ip] += np.log(s_1[ip])
            ln_s_2_geometric[ip] += np.log(s_2[ip])

    s_1_geometric = np.empty(numPaths)
    s_2_geometric = np.empty(numPaths)

    for ip in range(0, numPaths):
        s_1_arithmetic[ip] /= n
        s_1_geometric[ip] = np.exp(ln_s_1_geometric[ip] / n)
        s_2_arithmetic[ip] /= n
        s_2_geometric[ip] = np.exp(ln_s_2_geometric[ip] / n)

    if optionType == FinOptionTypes.EUROPEAN_CALL:
        payoff_a_1 = np.maximum(s_1_arithmetic - K, 0.0)
        payoff_g_1 = np.maximum(s_1_geometric - K, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - K, 0.0)
        payoff_g_2 = np.maximum(s_2_geometric - K, 0.0)
    elif optionType == FinOptionTypes.EUROPEAN_PUT:
        payoff_a_1 = np.maximum(K - s_1_arithmetic, 0.0)
        payoff_g_1 = np.maximum(K - s_1_geometric, 0.0)
        payoff_a_2 = np.maximum(K - s_2_arithmetic, 0.0)
        payoff_g_2 = np.maximum(K - s_2_geometric, 0.0)
    else:
        raise FinError("Unknown Option Type")

    payoff_a = np.concatenate((payoff_a_1, payoff_a_2), axis=0)
    payoff_g = np.concatenate((payoff_g_1, payoff_g_2), axis=0)

    m = covar(payoff_a, payoff_a)
    lam = m[0][1] / m[1][1]

    payoff_a_mean = np.mean(payoff_a)
    payoff_g_mean = np.mean(payoff_g)

    v_a = payoff_a_mean * exp(-r * t) * multiplier
    v_g = payoff_g_mean * exp(-r * t) * multiplier

    epsilon = v_g_exact - v_g
    v_a_cv = v_a + lam * epsilon

    return v_a_cv

##########################################################################


class FinEquityAsianOption(FinEquityOption):
    ''' Class to store an Equity Asian Option. '''

    def __init__(self,
                 startAveragingDate: FinDate,
                 expiryDate: FinDate,
                 strikePrice: Union[int, float],
                 optionType: FinOptionTypes,
                 numberOfObservations: int = 0):

        checkArgumentTypes(self.__init__, locals())

        if startAveragingDate > expiryDate:
            raise FinError("Averaging starts after expiry date")

        self._startAveragingDate = startAveragingDate
        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType
        self._numObservations = numberOfObservations

###############################################################################

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model,
              valuationMethod,
              accruedAverage=None):
        ''' Calculate the value of an Asian option. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if valuationMethod == "GEOMETRIC":

            v = self.valueGeometric(valueDate,
                                    stockPrice,
                                    discountCurve,
                                    dividendYield,
                                    model,
                                    accruedAverage)

        elif valuationMethod == "TURNBULL_WAKEMAN":

            v = self.valueTurnbullWakeman(valueDate,
                                          stockPrice,
                                          discountCurve,
                                          dividendYield,
                                          model,
                                          accruedAverage)

        elif valuationMethod == "CURRAN":

            v = self.valueCurran(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividendYield,
                                 model,
                                 accruedAverage)

#        elif valuationMethod == "MILEVSKY_POSNER":
#
#                v = self.valueMilevskyPosner(valueDate,
#                                   stockPrice,
#                                   dividendYield,
#                                   volatility,
#                                   interestRate,
#                                   accruedAverage)
#
#        elif valuationMethod == "LEVY":
#
#                v = self.valueLevy(valueDate,
#                                   stockPrice,
#                                   dividendYield,
#                                   volatility,
#                                   interestRate,
#                                   accruedAverage)

        else:
            raise FinError("Unknown valuation model")

        return v

###############################################################################

    def valueGeometric(self,
                       valueDate,
                       stockPrice,
                       discountCurve,
                       dividendYield,
                       model,
                       accruedAverage):

        # This is based on paper by Kemna and Vorst 1990. It calculates the
        # Geometric Asian option price which is a lower bound on the Arithmetic
        # option price.
        # This should not be used as a valuation model for the Arithmetic Average
        # option but can be used as a control variate for other approaches

        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate)  / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        df = discountCurve.df(t)
        r = -log(df)/t
        volatility = model._volatility

        K = self._strikePrice
        n = self._numObservations
        q = dividendYield
        S0 = stockPrice

        multiplier = 1.0

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(
                "In averaging period you need to enter the accrued average")

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / t
            # the number of options is rescaled also
            multiplier = t / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled
            n = n * t / tau

        sigSq = volatility ** 2
        meanGeo = (r - q - sigSq / 2.0) * (t0 + (t - t0) / 2.0)
        varGeo = sigSq * (t0 + (t - t0) * (2 * n - 1) / (6 * n))
        EG = S0 * exp(meanGeo + varGeo / 2.0)

        d1 = (meanGeo + log(S0 / K) + varGeo) / sqrt(varGeo)
        d2 = d1 - sqrt(varGeo)

        # the Geometric price is the lower bound
        call_g = exp(-r * t) * (EG * N(d1) - K * N(d2))

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = call_g
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            put_g = call_g - (EG - K) * exp(-r * t)
            v = put_g
        else:
            raise FinError("Unknown option type " + str(self._optionType))

        v = v * multiplier
        return v

###############################################################################

    def valueCurran(self,
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividendYield,
                    model,
                    accruedAverage):

        # This is based on paper by Vorst
        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate)  / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        multiplier = 1.0

        df = discountCurve.df(t)
        r = -log(df)/t
        volatility = model._volatility

        S0 = stockPrice
        b = r - dividendYield
        sigma2 = volatility**2
        K = self._strikePrice

        n = self._numObservations

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(
                "In averaging period you need to enter the accrued average")

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / t
            # the number of options is rescaled also
            multiplier = t / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * t / tau + 0.5) + 1

        h = (t - t0) / (n - 1)
        u = (1.0 - exp(b * h * n)) / (1.0 - exp(b * h))
        w = (1.0 - exp((2 * b + sigma2) * h * n)) / \
            (1.0 - exp((2 * b + sigma2) * h))

        FA = (S0 / n) * exp(b * t0) * u
        EA2 = (S0 * S0 / n / n) * exp((2.0 * b + sigma2) * t0)
        EA2 = EA2 * (w + 2.0 / (1 - exp((b + sigma2) * h)) * (u - w))
        sigmaA = sqrt((log(EA2) - 2.0 * log(FA)) / t)

        d1 = (log(FA / K) + sigmaA * sigmaA * t / 2.0) / (sigmaA * sqrt(t))
        d2 = d1 - sigmaA * sqrt(t)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = exp(-r * t) * (FA * N(d1) - K * N(d2))
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            v = exp(-r * t) * (K * N(-d2) - FA * N(-d1))
        else:
            return None

        v = v * multiplier
        return v

##########################################################################

    def valueTurnbullWakeman(self,
                             valueDate,
                             stockPrice,
                             discountCurve,
                             dividendYield,
                             model,
                             accruedAverage):

        # This is based on paper by Turnbull and Wakeman 1991 which uses
        # the edgeworth expansion to find the first two moments of the
        # arithmetic average

        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        multiplier = 1.0
        n = self._numObservations

        df = discountCurve.df(t)
        r = -log(df)/t

        volatility = model._volatility

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(
                "In averaging period you need to enter the accrued average")

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / t
            # the number of options is rescaled also
            multiplier = t / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * t / tau + 0.5) + 1

        # need to handle this
        b = r - dividendYield
        sigma2 = volatility**2
        a1 = b + sigma2
        a2 = 2 * b + sigma2
        S0 = stockPrice

        dt = t - t0

        if b == 0:
            M1 = 1.0
            M2 = 2.0 * exp(sigma2 * t) - 2.0 * \
                exp(sigma2 * t0) * (1.0 + sigma2 * dt)
            M2 = M2 / sigma2 / sigma2 / dt / dt
        else:
            M1 = S0 * (exp(b * t) - exp(b * t0)) / (b * dt)
            M2 = exp(a2 * t) / a1 / a2 / dt / dt + \
                (exp(a2 * t0) / b / dt / dt) * (1.0 / a2 - exp(b * dt) / a1)
            M2 = 2.0 * M2 * S0 * S0

        F0 = M1
        sigma2 = 1.0 / t * log(M2 / M1 / M1)
        sigma = sqrt(sigma2)

        d1 = (log(F0 / K) + sigma2 * t / 2) / sigma / sqrt(t)
        d2 = d1 - sigma * sqrt(t)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            call = exp(-r * t) * (F0 * N(d1) - K * N(d2))
            v = call
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            put = exp(-r * t) * (K * N(-d2) - F0 * N(-d1))
            v = put
        else:
            return None

        v = v * multiplier

        return v

##############################################################################

    def valueMC(self,
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model,
                numPaths,
                seed,
                accruedAverage):

        # Basic validation
        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        if valueDate > self._startAveragingDate and accruedAverage is None:
            raise FinError(
                "In averaging period so need to enter accrued average.")

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate)  / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        df = discountCurve.df(t)
        r = -log(df)/t
        volatility = model._volatility

        K = self._strikePrice
        n = self._numObservations

        v = valueMC_NUMBA(t0, t, tau, K, n, self._optionType,
                          stockPrice,
                          r,
                          dividendYield,
                          volatility,
                          numPaths,
                          seed,
                          accruedAverage)

        return v

###############################################################################

    def valueMC_fast(self,
                     valueDate,
                     stockPrice,
                     discountCurve,
                     dividendYield,
                     model,
                     numPaths,
                     seed,
                     accruedAverage):

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        n = self._numObservations

        df = discountCurve.df(t)
        r = -log(df)/t
        volatility = model._volatility

        v = valueMC_fast_NUMBA(t0, t, tau, K, n, self._optionType,
                               stockPrice,
                               r,
                               dividendYield,
                               volatility,
                               numPaths,
                               seed,
                               accruedAverage)

        return v

##########################################################################

    def valueMC_fast_CV(self,
                        valueDate,
                        stockPrice,
                        discountCurve,
                        dividendYield,
                        model,
                        numPaths,
                        seed,
                        accruedAverage):

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        t = (self._expiryDate - valueDate)  / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        n = self._numObservations

        df = discountCurve.df(t)
        r = -log(df)/t
        volatility = model._volatility

        # For control variate we price a Geometric average option exactly
        v_g_exact = self.valueGeometric(valueDate,
                                        stockPrice,
                                        discountCurve,
                                        dividendYield,
                                        model,
                                        accruedAverage)

        v = valueMC_fast_CV_NUMBA(t0, t, tau, K, n, self._optionType,
                                  stockPrice,
                                  r,
                                  dividendYield,
                                  volatility,
                                  numPaths,
                                  seed,
                                  accruedAverage,
                                  v_g_exact)

        return v

##########################################################################
