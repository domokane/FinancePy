##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit

# TODO: Add perturbatory risk using the analytical methods !!
# TODO: Add Sobol to Monte Carlo

from ...finutils.FinMath import covar
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError

from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...finutils.FinDate import FinDate
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...finutils.FinMath import N


###############################################################################

from enum import Enum


class FinAsianOptionValuationMethods(Enum):
    GEOMETRIC = 1,
    TURNBULL_WAKEMAN = 2,
    CURRAN = 3

###############################################################################


errorStr = "In averaging period so need to enter accrued average."

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
def _valueMC_NUMBA(t0,
                   t,
                   tau,
                   K,
                   n,
                   optionType,
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

    if t0 < 0.0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(errorStr)

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

    for _ in range(0, numPaths):

        # evolve stock price to start of averaging period
        g = np.random.standard_normal(1)
        s_1 = stockPrice * np.exp((mu - v2 / 2.0) * t0 +
                                  g[0] * np.sqrt(t0) * volatility)
        s_2 = stockPrice * np.exp((mu - v2 / 2.0) * t0 -
                                  g[0] * np.sqrt(t0) * volatility)

        # enter averaging period
        s_1_arithmetic = 0.0
        s_2_arithmetic = 0.0

        g = np.random.standard_normal(n)

        for obs in range(0, n):

            s_1 = s_1 * np.exp((mu - v2 / 2.0) * dt +
                               g[obs] * np.sqrt(dt) * volatility)

            s_2 = s_2 * np.exp((mu - v2 / 2.0) * dt -
                               g[obs] * np.sqrt(dt) * volatility)

            s_1_arithmetic += s_1
            s_2_arithmetic += s_2

        s_1_arithmetic /= n
        s_2_arithmetic /= n

        if optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff_a += max(s_1_arithmetic - K, 0.0)
            payoff_a += max(s_2_arithmetic - K, 0.0)
        elif optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff_a += max(K - s_1_arithmetic, 0.0)
            payoff_a += max(K - s_2_arithmetic, 0.0)
        else:
            return None

    v_a = payoff_a * np.exp(-interestRate * t) / numPaths / 2.0
    v_a = v_a * multiplier
    return v_a

###############################################################################


@njit(cache=True, fastmath=True)
def _valueMC_fast_NUMBA(t0: float,
                        t: float,
                        tau: float,
                        K: float,
                        n: int,
                        optionType: FinOptionTypes,
                        stockPrice: float,
                        interestRate: float,
                        dividendYield: float,
                        volatility: float,
                        numPaths: int,
                        seed: int,
                        accruedAverage: float):

    np.random.seed(seed)
    mu = interestRate - dividendYield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interestRate
    numPaths = int(numPaths)

    multiplier = 1.0

    if t0 < 0.0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(errorStr)

        # we adjust the strike to account for the accrued coupon
        K = (K * tau + accruedAverage * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    # g = np.random.normal(0.0, 1.0, size=(numPaths))

    gg = np.empty(numPaths, np.float64)
    for iPath in range(0, numPaths):
        rv = np.random.normal()
        gg[iPath] = rv

    s_1 = np.empty(numPaths)
    s_2 = np.empty(numPaths)

    for ip in range(0, numPaths):
        s_1[ip] = stockPrice * \
            np.exp((mu - v2 / 2.0) * t0 + gg[ip] * np.sqrt(t0) * volatility)
        s_2[ip] = stockPrice * \
            np.exp((mu - v2 / 2.0) * t0 - gg[ip] * np.sqrt(t0) * volatility)

    s_1_arithmetic = np.zeros(numPaths)
    s_2_arithmetic = np.zeros(numPaths)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=(numPaths))

        for ip in range(0, numPaths):
            s_1[ip] = s_1[ip] * np.exp((mu - v2 / 2.0) *
                                       dt + g[ip] * np.sqrt(dt) * volatility)
            s_2[ip] = s_2[ip] * np.exp((mu - v2 / 2.0) *
                                       dt - g[ip] * np.sqrt(dt) * volatility)

        for ip in range(0, numPaths):
            s_1_arithmetic[ip] += s_1[ip] / n
            s_2_arithmetic[ip] += s_2[ip] / n

    if optionType == FinOptionTypes.EUROPEAN_CALL:
        payoff_a_1 = np.maximum(s_1_arithmetic - K, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - K, 0.0)
    elif optionType == FinOptionTypes.EUROPEAN_PUT:
        payoff_a_1 = np.maximum(K - s_1_arithmetic, 0.0)
        payoff_a_2 = np.maximum(K - s_2_arithmetic, 0.0)
    else:
        raise FinError("Unknown option type.")

    payoff_a = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    v_a = multiplier * payoff_a * np.exp(- r * t) / 2.0
    return v_a

###############################################################################


@njit(cache=True, fastmath=True)
def _valueMC_fast_CV_NUMBA(t0, t, tau, K, n, optionType, stockPrice,
                           interestRate, dividendYield, volatility, numPaths,
                           seed, accruedAverage, v_g_exact):

    np.random.seed(seed)
    mu = interestRate - dividendYield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interestRate

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accruedAverage is None:
            raise FinError(errorStr)

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
            np.exp((mu - v2 / 2.0) * t0 + g[ip] * np.sqrt(t0) * volatility)
        s_2[ip] = stockPrice * \
            np.exp((mu - v2 / 2.0) * t0 - g[ip] * np.sqrt(t0) * volatility)

    s_1_arithmetic = np.zeros(numPaths)
    s_2_arithmetic = np.zeros(numPaths)
    ln_s_1_geometric = np.zeros(numPaths)
    ln_s_2_geometric = np.zeros(numPaths)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=(numPaths))
        for ip in range(0, numPaths):
            s_1[ip] = s_1[ip] * np.exp((mu - v2 / 2.0) *
                                       dt + g[ip] * np.sqrt(dt) * volatility)
            s_2[ip] = s_2[ip] * np.exp((mu - v2 / 2.0) *
                                       dt - g[ip] * np.sqrt(dt) * volatility)

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

    # Now we do the control variate adjustment
    m = covar(payoff_a, payoff_a)
    lam = m[0][1] / m[1][1]

    payoff_a_mean = np.mean(payoff_a)
    payoff_g_mean = np.mean(payoff_g)

    v_a = payoff_a_mean * np.exp(-r * t) * multiplier
    v_g = payoff_g_mean * np.exp(-r * t) * multiplier

    epsilon = v_g_exact - v_g
    v_a_cv = v_a + lam * epsilon

    return v_a_cv

###############################################################################


class FinEquityAsianOption():
    ''' Class for an Equity Asian Option. This is an option with a final payoff
    linked to the averaging of the stock price over some specified period
    before the option expires. The valuation is done for both an arithmetic and
    a geometric average but the former can only be done either using an
    analytical approximation of the arithmetic average distribution or by using
    Monte-Carlo simulation.'''

    def __init__(self,
                 startAveragingDate: FinDate,
                 expiryDate: FinDate,
                 strikePrice: float,
                 optionType: FinOptionTypes,
                 numberOfObservations: int = 100):
        ''' Create an FinEquityAsian option object which takes a start date for
        the averaging, an expiry date, a strike price, an option type and a
        number of observations. '''

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
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              model,
              method: FinAsianOptionValuationMethods,
              accruedAverage: float = None):
        ''' Calculate the value of an Asian option using one of the specified
        analytical approximations for an average rate option. These are the
        three enumerated values in the enum FinAsianOptionValuationMethods. The
        choices of approximation are (i) GEOMETRIC - the average is a geometric
        one as in paper by Kenna and Worst (1990), (ii) TURNBULL_WAKEMAN -
        this is a value based on an edgeworth expansion of the moments of the
        arithmetic average, and (iii) CURRAN - another approximative approach
        by Curran based on conditioning on the geometric mean price. Just
        choose the corresponding enumerated value to switch between these
        different approaches.

        Note that the accrued average is only required if the value date is
        inside the averaging period for the option. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        if method == FinAsianOptionValuationMethods.GEOMETRIC:
            v = self._valueGeometric(valueDate,
                                     stockPrice,
                                     discountCurve,
                                     dividendCurve,
                                     model,
                                     accruedAverage)

        elif method == FinAsianOptionValuationMethods.TURNBULL_WAKEMAN:
            v = self._valueTurnbullWakeman(valueDate,
                                           stockPrice,
                                           discountCurve,
                                           dividendCurve,
                                           model,
                                           accruedAverage)

        elif method == FinAsianOptionValuationMethods.CURRAN:
            v = self._valueCurran(valueDate,
                                  stockPrice,
                                  discountCurve,
                                  dividendCurve,
                                  model,
                                  accruedAverage)
        else:
            raise FinError("Unknown valuation model")

        return v

###############################################################################

    def _valueGeometric(self, valueDate, stockPrice, discountCurve,
                        dividendCurve, model, accruedAverage):
        ''' This option valuation is based on paper by Kemna and Vorst 1990. It
        calculates the Geometric Asian option price which is a lower bound on
        the Arithmetic option price. This should not be used as a valuation
        model for the Arithmetic Average option but can be used as a control
        variate for other approaches. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

#        print("r:", r, "q:", q)
        
        volatility = model._volatility

        K = self._strikePrice
        n = self._numObservations
        S0 = stockPrice

        multiplier = 1.0

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(errorStr)

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / texp
            # the number of options is rescaled also
            multiplier = texp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled
            n = n * texp / tau

        sigSq = volatility ** 2
        meanGeo = (r - q - sigSq / 2.0) * (t0 + (texp - t0) / 2.0)
        varGeo = sigSq * (t0 + (texp - t0) * (2 * n - 1) / (6 * n))
        EG = S0 * np.exp(meanGeo + varGeo / 2.0)

        d1 = (meanGeo + np.log(S0 / K) + varGeo) / np.sqrt(varGeo)
        d2 = d1 - np.sqrt(varGeo)

        # the Geometric price is the lower bound
        call_g = np.exp(-r * texp) * (EG * N(d1) - K * N(d2))

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = call_g
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            put_g = call_g - (EG - K) * np.exp(-r * texp)
            v = put_g
        else:
            raise FinError("Unknown option type " + str(self._optionType))

        v = v * multiplier
        return v

###############################################################################

    def _valueCurran(self, valueDate, stockPrice, discountCurve,
                     dividendCurve, model, accruedAverage):
        ''' Valuation of an Asian option using the result by Vorst. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        multiplier = 1.0

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        volatility = model._volatility

        S0 = stockPrice
        b = r - q
        sigma2 = volatility**2
        K = self._strikePrice

        n = self._numObservations

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(errorStr)

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / texp
            # the number of options is rescaled also
            multiplier = texp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * texp / tau + 0.5) + 1

        h = (texp - t0) / (n - 1)
        u = (1.0 - np.exp(b * h * n)) / (1.0 - np.exp(b * h))
        w = (1.0 - np.exp((2 * b + sigma2) * h * n)) / \
            (1.0 - np.exp((2 * b + sigma2) * h))

        FA = (S0 / n) * np.exp(b * t0) * u
        EA2 = (S0 * S0 / n / n) * np.exp((2.0 * b + sigma2) * t0)
        EA2 = EA2 * (w + 2.0 / (1.0 - np.exp((b + sigma2) * h)) * (u - w))
        sigmaA = np.sqrt((np.log(EA2) - 2.0 * np.log(FA)) / texp)

        d1 = (np.log(FA / K) + sigmaA * sigmaA * texp / 2.0) / (sigmaA*np.sqrt(texp))
        d2 = d1 - sigmaA * np.sqrt(texp)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = np.exp(-r * texp) * (FA * N(d1) - K * N(d2))
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            v = np.exp(-r * texp) * (K * N(-d2) - FA * N(-d1))
        else:
            return None

        v = v * multiplier
        return v

###############################################################################

    def _valueTurnbullWakeman(self, valueDate, stockPrice, discountCurve,
                              dividendCurve, model, accruedAverage):
        ''' Asian option valuation based on paper by Turnbull and Wakeman 1991
        which uses the edgeworth expansion to find the first two moments of the
        arithmetic average. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        multiplier = 1.0
        n = self._numObservations

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        volatility = model._volatility

        if t0 < 0:  # we are in the averaging period

            if accruedAverage is None:
                raise FinError(errorStr)

            # we adjust the strike to account for the accrued coupon
            K = (K * tau + accruedAverage * t0) / texp
            # the number of options is rescaled also
            multiplier = texp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * texp / tau + 0.5) + 1

        # need to handle this
        b = r - q
        sigma2 = volatility**2
        a1 = b + sigma2
        a2 = 2 * b + sigma2
        S0 = stockPrice

        dt = texp - t0

        if b == 0:
            M1 = 1.0
            M2 = 2.0 * np.exp(sigma2 * texp) - 2.0 * \
                np.exp(sigma2 * t0) * (1.0 + sigma2 * dt)
            M2 = M2 / sigma2 / sigma2 / dt / dt
        else:
            M1 = S0 * (np.exp(b * texp) - np.exp(b * t0)) / (b * dt)
            M2 = np.exp(a2 * texp) / a1 / a2 / dt / dt + \
                (np.exp(a2 * t0) / b / dt / dt) * (1.0/a2 - np.exp(b*dt) / a1)
            M2 = 2.0 * M2 * S0 * S0

        F0 = M1
        sigma2 = (1.0 / texp) * np.log(M2 / M1 / M1)
        sigma = np.sqrt(sigma2)

        d1 = (np.log(F0 / K) + sigma2 * texp / 2) / sigma / np.sqrt(texp)
        d2 = d1 - sigma * np.sqrt(texp)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            call = np.exp(-r * texp) * (F0 * N(d1) - K * N(d2))
            v = call
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            put = np.exp(-r * texp) * (K * N(-d2) - F0 * N(-d1))
            v = put
        else:
            return None

        v = v * multiplier
        return v

###############################################################################

    def _valueMC(self,
                 valueDate: FinDate,
                 stockPrice: float,
                 discountCurve: FinDiscountCurve,
                 dividendCurve: FinDiscountCurve,
                 model,
                 numPaths: int,
                 seed: int,
                 accruedAverage: float):
        ''' Monte Carlo valuation of the Asian Average option using standard
        Monte Carlo code enhanced by Numba. I have discontinued the use of this
        as it is both slow and has limited variance reduction. '''

        # Basic validation
        if valueDate > self._expiryDate:
            raise FinError("Value date after option expiry date.")

        if valueDate > self._startAveragingDate and accruedAverage is None:
            raise FinError(errorStr)

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        volatility = model._volatility

        K = self._strikePrice
        n = self._numObservations

        v = _valueMC_NUMBA(t0, texp, tau, K, n, self._optionType,
                           stockPrice,
                           r,
                           q,
                           volatility,
                           numPaths,
                           seed,
                           accruedAverage)

        return v

##############################################################################

    def _valueMC_fast(self,
                      valueDate,
                      stockPrice,
                      discountCurve,
                      dividendCurve,  # Yield
                      model,          # Model
                      numPaths,       # Numpaths integer
                      seed,
                      accruedAverage):
        ''' Monte Carlo valuation of the Asian Average option. This method uses
        a lot of Numpy vectorisation. It is also helped by Numba. '''

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        n = self._numObservations

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        volatility = model._volatility

        v = _valueMC_fast_NUMBA(t0, texp, tau,
                                K, n, self._optionType,
                                stockPrice,
                                r,
                                q,
                                volatility,
                                numPaths,
                                seed,
                                accruedAverage)

        return v

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendCurve: FinDiscountCurve,
                model,
                numPaths: int,
                seed: int,
                accruedAverage: float):
        ''' Monte Carlo valuation of the Asian Average option using a control
        variate method that improves accuracy and reduces the variance of the
        price. This uses Numpy and Numba. This is the standard MC pricer. '''

        # the years to the start of the averaging period
        t0 = (self._startAveragingDate - valueDate) / gDaysInYear
        texp = (self._expiryDate - valueDate) / gDaysInYear
        tau = (self._expiryDate - self._startAveragingDate) / gDaysInYear

        K = self._strikePrice
        n = self._numObservations

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        volatility = model._volatility

        # For control variate we price a Geometric average option exactly
        v_g_exact = self._valueGeometric(valueDate,
                                         stockPrice,
                                         discountCurve,
                                         dividendCurve,
                                         model,
                                         accruedAverage)

        v = _valueMC_fast_CV_NUMBA(t0,
                                   texp,
                                   tau,
                                   K,
                                   n,
                                   self._optionType,
                                   stockPrice,
                                   r,
                                   q,
                                   volatility,
                                   numPaths,
                                   seed,
                                   accruedAverage,
                                   v_g_exact)

        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START AVERAGING DATE", self._startAveragingDate)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUM OBSERVATIONS", self._numObservations, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
