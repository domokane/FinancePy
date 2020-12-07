##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import vectorize

from scipy import optimize


from ...finutils.FinDate import FinDate
from ...finutils.FinMath import nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...models.FinModel import FinModel
from ...models.FinModelBlackScholes import FinModelBlackScholes

from ...models.FinModelBlackScholesMC import _valueMC_NONUMBA_NONUMPY
from ...models.FinModelBlackScholesMC import _valueMC_NUMPY_NUMBA
from ...models.FinModelBlackScholesMC import _valueMC_NUMBA_ONLY
from ...models.FinModelBlackScholesMC import _valueMC_NUMPY_ONLY
from ...models.FinModelBlackScholesMC import _valueMC_NUMBA_PARALLEL

from ...finutils.FinMath import N


###############################################################################


def _f(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    discountCurve = args[3]
    dividendYield = args[4]
    price = args[5]

    model = FinModelBlackScholes(volatility)

    objFn = self.value(valueDate,
                       stockPrice,
                       discountCurve,
                       dividendYield,
                       model)

    objFn = objFn - price

    return objFn

###############################################################################


def _fvega(volatility, *args):

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    discountCurve = args[3]
    dividendYield = args[4]

    model = FinModelBlackScholes(volatility)

    fprime = self.vega(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
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
        valid_types = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]

        if isinstance(optionType, list):
            op_type_filter = lambda x: 1 if x in valid_types else 0
            if not all(map(op_type_filter, optionType)):
                raise FinError("Unknown Option Type" + str(optionType))
        elif optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        self._expiryDate = expiryDate
        self._strikePrice = strikePrice
        self._optionType = optionType
        self._numOptions = numOptions

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

        if type(model) == FinModelBlackScholes:

            v_opt = model.value(s0, texp, k, r, q, self._optionType)

        else:
            raise FinError("Unknown Model Type")

        v = v_opt * self._numOptions
        return v

###############################################################################

    def delta(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):
        ''' Calculate the analytical delta of a European vanilla option. '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        df = discountCurve._df(t)
        r = -np.log(df)/t
        q = dividendYield

        if type(model) == FinModelBlackScholes:

            volatility = model._volatility

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(stockPrice / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = r - q
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                delta = np.exp(-q * t) * N(d1)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                delta = -np.exp(-q * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

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
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        df = discountCurve._df(t)
        interestRate = -np.log(df)/t

        if isinstance(model, FinModelBlackScholes):

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(stockPrice / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = interestRate - dividendYield
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            gamma = np.exp(-dividendYield * t) * nprime(d1) / stockPrice / den

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
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        df = discountCurve._df(t)
        interestRate = -np.log(df)/t

        if isinstance(model, FinModelBlackScholes):

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(stockPrice / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = interestRate - dividendYield
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            vega = stockPrice * sqrtT * np.exp(-dividendYield * t) * nprime(d1)
        else:
            raise FinError("Unknown Model type")

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
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        df = discountCurve._df(t)
        interestRate = -np.log(df)/t

        if isinstance(model, FinModelBlackScholes):

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(stockPrice / self._strikePrice)
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = interestRate - dividendYield
            v2 = volatility * volatility
            d1 = (lnS0k + (mu + v2 / 2.0) * t) / den
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = - stockPrice * np.exp(-dividendYield * t) * \
                    nprime(d1) * volatility / 2.0 / sqrtT
                v = v - interestRate * self._strikePrice * \
                    df * N(d2)
                v = v + dividendYield * stockPrice * \
                    np.exp(-dividendYield * t) * N(d1)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = - stockPrice * np.exp(-dividendYield * t) * \
                    nprime(d1) * volatility / 2.0 / sqrtT
                v = v + interestRate * self._strikePrice * \
                    df * N(-d2)
                v = v - dividendYield * stockPrice * \
                    np.exp(-dividendYield * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v

###############################################################################

    def rho(self,
            valueDate: FinDate,
            stockPrice: float,
            discountCurve: FinDiscountCurve,
            dividendYield: float,
            model:FinModel):
        ''' Calculate the analytical rho of a European vanilla option. '''

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        df = discountCurve._df(t)
        interestRate = -np.log(df)/t

        if isinstance(model, FinModelBlackScholes):

            volatility = model._volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            lnS0k = np.log(stockPrice / self._strikePrice)
            K = self._strikePrice
            sqrtT = np.sqrt(t)
            den = volatility * sqrtT
            mu = interestRate - dividendYield
            v2 = volatility * volatility
            d2 = (lnS0k + (mu - v2 / 2.0) * t) / den

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = K * t * df * N(d2)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = -K * t * df * N(-d2)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v
###############################################################################

    def impliedVolatility(self,
                          valueDate: FinDate,
                          stockPrice: (float, list, np.ndarray),
                          discountCurve: FinDiscountCurve,
                          dividendYield: float,
                          price):
        ''' Calculate the implied volatility of a European vanilla option. '''

        argtuple = (self, valueDate, stockPrice,
                    discountCurve, dividendYield, price)

        sigma = optimize.newton(_f, x0=0.2, fprime=_fvega, args=argtuple,
                                tol=1e-5, maxiter=50, fprime2=None)
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
