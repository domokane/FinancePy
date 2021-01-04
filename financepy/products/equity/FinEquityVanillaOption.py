##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit

# from scipy import optimize
from ...finutils.FinSolvers import newton_secant, bisection, newton

from ...finutils.FinDate import FinDate
from ...finutils.FinMath import nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...models.FinModel import FinModel
from ...models.FinModelBlackScholes import FinModelBlackScholes
from ...models.FinModelBlackScholesAnalytical import bsValue
from ...models.FinModelBlackScholesAnalytical import bsDelta
from ...models.FinModelBlackScholesAnalytical import bsVega
from ...models.FinModelBlackScholesAnalytical import bsGamma
from ...models.FinModelBlackScholesAnalytical import bsRho
from ...models.FinModelBlackScholesAnalytical import bsTheta

from ...models.FinModelBlackScholesMC import _valueMC_NONUMBA_NONUMPY
from ...models.FinModelBlackScholesMC import _valueMC_NUMPY_NUMBA
from ...models.FinModelBlackScholesMC import _valueMC_NUMBA_ONLY
from ...models.FinModelBlackScholesMC import _valueMC_NUMPY_ONLY
from ...models.FinModelBlackScholesMC import _valueMC_NUMBA_PARALLEL

from ...finutils.FinMath import N

###############################################################################

@njit(fastmath=True, cache=True)
def _f(v, args):

    optionTypeValue = int(args[0])
    texp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]
    price = args[6]

    objFn = bsValue(s0, texp, k, r, q, v, optionTypeValue)
    objFn = objFn - price
    return objFn

###############################################################################


def _fvega(v, *args):

    self = args[0]
    texp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]

    fprime = bsVega(s0, texp, k, r, q, v, self._optionType.value)
    return fprime

###############################################################################


class FinEquityVanillaOption():
    ''' Class for managing plain vanilla European calls and puts on equities.
    For American calls and puts see the FinEquityAmericanOption class. '''

    def __init__(self,
                 expiryDate: (FinDate, list),
                 strikePrice: (float, np.ndarray),
                 optionType: (FinOptionTypes, list),
                 numOptions: float = 1.0):
        ''' Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options. '''

        checkArgumentTypes(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        self._expiryDate = expiryDate
        self._strikePrice = strikePrice
        self._optionType = optionType
        self._numOptions = numOptions
        self._texp = None

###############################################################################

    def value(self,
              valueDate: (FinDate, list),
              stockPrice: (np.ndarray, float),
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model: FinModel):
        ''' Equity Vanilla Option valuation using Black-Scholes model. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        self._texp = texp

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            value = bsValue(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        value = value * self._numOptions
        return value

###############################################################################

    def delta(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):
        ''' Calculate the analytical delta of a European vanilla option. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        self._texp = texp

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            delta = bsDelta(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return delta

###############################################################################

    def gamma(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model:FinModel):
        ''' Calculate the analytical gamma of a European vanilla option. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            gamma = bsGamma(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return gamma

###############################################################################

    def vega(self,
             valueDate: FinDate,
             stockPrice: float,
             discountCurve: FinDiscountCurve,
             dividendYield: float,
             model:FinModel):
        ''' Calculate the analytical vega of a European vanilla option. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            vega = bsVega(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return vega

###############################################################################

    def theta(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model:FinModel):
        ''' Calculate the analytical theta of a European vanilla option. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            theta = bsTheta(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return theta

###############################################################################

    def rho(self,
            valueDate: FinDate,
            stockPrice: float,
            discountCurve: FinDiscountCurve,
            dividendYield: float,
            model:FinModel):
        ''' Calculate the analytical rho of a European vanilla option. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)
        df = discountCurve.df(self._expiryDate)
        s0 = stockPrice
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            rho = bsRho(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return rho

###############################################################################

    def impliedVolatility(self,
                          valueDate: FinDate,
                          stockPrice: (float, list, np.ndarray),
                          discountCurve: FinDiscountCurve,
                          dividendYield: float,
                          price):
        ''' Calculate the implied volatility of a European vanilla option. '''

        texp = (self._expiryDate - valueDate) / gDaysInYear

        if texp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999

        if price < 1e-10:
            print("Option value is effectively zero.")
            return -999.0

        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/texp
        q = dividendYield
        k = self._strikePrice
        s0 = stockPrice

        if np.abs(k-s0)/ (k+s0) < 0.05:
            sigma0 = price / 0.4 / stockPrice / np.sqrt(texp)
        else:
            sigma0 = 0.20

        # NEED TO MAP THE OPTION TO AN OTM option!!!

        optionTypeValue = self._optionType.value

        argsv = np.array([optionTypeValue, texp, s0, r, q, k, price])

#        sigma = newton(_f, x0=sigma0, args=argsv, tol=1e-6, maxiter=100)

        sigma = bisection(_f, 0.0, 10.0, args=argsv, xtol=1e-6, maxIter=100)

        return sigma

###############################################################################

    def valueMC_NUMPY_ONLY(self,
                           valueDate: FinDate,
                           stockPrice: float,
                           discountCurve: FinDiscountCurve,
                           dividendYield: float,
                           model:FinModel,
                           numPaths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NUMPY_ONLY(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC_NUMBA_ONLY(self,
                           valueDate: FinDate,
                           stockPrice: float,
                           discountCurve: FinDiscountCurve,
                           dividendYield: float,
                           model:FinModel,
                           numPaths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NUMBA_ONLY(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed, 
                           useSobol)

        return v

###############################################################################

    def valueMC_NUMBA_PARALLEL(self,
                           valueDate: FinDate,
                           stockPrice: float,
                           discountCurve: FinDiscountCurve,
                           dividendYield: float,
                           model:FinModel,
                           numPaths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NUMBA_PARALLEL(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed, 
                           useSobol)

#        _valueMC_NUMBA_ONLY.parallel_diagnostics(level=4)

        return v

###############################################################################

    def valueMC_NUMPY_NUMBA(self,
                      valueDate: FinDate,
                      stockPrice: float,
                      discountCurve: FinDiscountCurve,
                      dividendYield: float,
                      model:FinModel,
                      numPaths: int = 10000,
                      seed: int = 4242,
                      useSobol: int = 0):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NUMPY_NUMBA(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC_NONUMBA_NONUMPY(self,
                      valueDate: FinDate,
                      stockPrice: float,
                      discountCurve: FinDiscountCurve,
                      dividendYield: float,
                      model:FinModel,
                      numPaths: int = 10000,
                      seed: int = 4242,
                      useSobol: int = 0):

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NONUMBA_NONUMPY(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendYield: float,
                model:FinModel,
                numPaths: int = 10000,
                seed: int = 4242,
                useSobol: int = 0):
        ''' Value European style call or put option using Monte Carlo. This is
        mainly for educational purposes. Sobol numbers can be used. '''

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t
        vol = model._volatility

        v = _valueMC_NUMBA_ONLY(stockPrice, 
                           t, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           dividendYield, 
                           vol, 
                           numPaths, 
                           seed, 
                           useSobol)

        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUMBER", self._numOptions, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
