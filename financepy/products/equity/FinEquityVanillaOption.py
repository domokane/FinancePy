##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit

# from scipy import optimize
from ...utils.FinSolvers1D import newton_secant, bisection, newton

from ...utils.date import Date
from ...utils.global_variables import gDaysInYear
from ...utils.FinError import FinError
from ...utils.FinGlobalTypes import FinOptionTypes
from ...utils.helper_functions import check_argument_types, labelToString
from ...market.curves.discount_curve import DiscountCurve

from ...models.FinModel import FinModel
from ...models.black_scholes import FinModelBlackScholes
from ...models.black_scholes_analytic import bsValue
from ...models.black_scholes_analytic import bsDelta
from ...models.black_scholes_analytic import bsVega
from ...models.black_scholes_analytic import bsGamma
from ...models.black_scholes_analytic import bsRho
from ...models.black_scholes_analytic import bsTheta
from ...models.black_scholes_analytic import bsImpliedVolatility
from ...models.black_scholes_analytic import bsIntrinsic


from ...models.black_scholes_mc import _valueMC_NONUMBA_NONUMPY
from ...models.black_scholes_mc import _valueMC_NUMPY_NUMBA
from ...models.black_scholes_mc import _valueMC_NUMBA_ONLY
from ...models.black_scholes_mc import _valueMC_NUMPY_ONLY
from ...models.black_scholes_mc import _valueMC_NUMBA_PARALLEL

from ...utils.fin_math import N

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
    """ Class for managing plain vanilla European calls and puts on equities.
    For American calls and puts see the FinEquityAmericanOption class. """

    def __init__(self,
                 expiry_date: (Date, list),
                 strikePrice: (float, np.ndarray),
                 optionType: (FinOptionTypes, list),
                 numOptions: float = 1.0):
        """ Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options. """

        check_argument_types(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        self._expiry_date = expiry_date
        self._strikePrice = strikePrice
        self._optionType = optionType
        self._numOptions = numOptions
        self._texp = None

###############################################################################

    def intrinsic(self,
                  valuation_date: (Date, list),
                  stock_price: (np.ndarray, float),
                  discount_curve: DiscountCurve,
                  dividendCurve: DiscountCurve):
        """ Equity Vanilla Option valuation using Black-Scholes model. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        self._texp = texp

        s0 = stock_price
        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        intrinsicValue = bsIntrinsic(s0, texp, k, r, q,
                                     self._optionType.value)

        intrinsicValue = intrinsicValue * self._numOptions
        return intrinsicValue

###############################################################################

    def value(self,
              valuation_date: (Date, list),
              stock_price: (np.ndarray, float),
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model: FinModel):
        """ Equity Vanilla Option valuation using Black-Scholes model. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        self._texp = texp

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price

        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

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
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model):
        """ Calculate the analytical delta of a European vanilla option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        self._texp = texp

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            delta = bsDelta(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return delta

###############################################################################

    def gamma(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model:FinModel):
        """ Calculate the analytical gamma of a European vanilla option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price

        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            gamma = bsGamma(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return gamma

###############################################################################

    def vega(self,
             valuation_date: Date,
             stock_price: float,
             discount_curve: DiscountCurve,
             dividendCurve: DiscountCurve,
             model:FinModel):
        """ Calculate the analytical vega of a European vanilla option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            vega = bsVega(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return vega

###############################################################################

    def theta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model:FinModel):
        """ Calculate the analytical theta of a European vanilla option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            theta = bsTheta(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return theta

###############################################################################

    def rho(self,
            valuation_date: Date,
            stock_price: float,
            discount_curve: DiscountCurve,
            dividendCurve: DiscountCurve,
            model:FinModel):
        """ Calculate the analytical rho of a European vanilla option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        texp = np.maximum(texp, 1e-10)

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            rho = bsRho(s0, texp, k, r, q, v, self._optionType.value)

        else:
            raise FinError("Unknown Model Type")

        return rho

###############################################################################

    def impliedVolatility(self,
                          valuation_date: Date,
                          stock_price: (float, list, np.ndarray),
                          discount_curve: DiscountCurve,
                          dividendCurve: DiscountCurve,
                          price):
        """ Calculate the Black-Scholes implied volatility of a European 
        vanilla option. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        if texp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strikePrice
        s0 = stock_price

        sigma = bsImpliedVolatility(s0, texp, k, r, q, price, 
                                    self._optionType.value)
        
        return sigma

###############################################################################

    def valueMC_NUMPY_ONLY(self,
                           valuation_date: Date,
                           stock_price: float,
                           discount_curve: DiscountCurve,
                           dividendCurve: DiscountCurve,
                           model:FinModel,
                           num_paths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NUMPY_ONLY(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC_NUMBA_ONLY(self,
                           valuation_date: Date,
                           stock_price: float,
                           discount_curve: DiscountCurve,
                           dividendCurve: DiscountCurve,
                           model:FinModel,
                           num_paths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NUMBA_ONLY(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

        return v

###############################################################################

    def valueMC_NUMBA_PARALLEL(self,
                               valuation_date: Date,
                               stock_price: float,
                               discount_curve: DiscountCurve,
                               dividendCurve: DiscountCurve,
                               model:FinModel,
                               num_paths: int = 10000,
                               seed: int = 4242,
                               useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NUMBA_PARALLEL(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

#        _valueMC_NUMBA_ONLY.parallel_diagnostics(level=4)

        return v

###############################################################################

    def valueMC_NUMPY_NUMBA(self,
                            valuation_date: Date,
                            stock_price: float,
                            discount_curve: DiscountCurve,
                            dividendCurve: DiscountCurve,
                            model:FinModel,
                            num_paths: int = 10000,
                            seed: int = 4242,
                            useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NUMPY_NUMBA(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC_NONUMBA_NONUMPY(self,
                                valuation_date: Date,
                                stock_price: float,
                                discount_curve: DiscountCurve,
                                dividendCurve: DiscountCurve,
                                model:FinModel,
                                num_paths: int = 10000,
                                seed: int = 4242,
                                useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NONUMBA_NONUMPY(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def valueMC(self,
                valuation_date: Date,
                stock_price: float,
                discount_curve: DiscountCurve,
                dividendCurve: DiscountCurve,
                model:FinModel,
                num_paths: int = 10000,
                seed: int = 4242,
                useSobol: int = 0):
        """ Value European style call or put option using Monte Carlo. This is
        mainly for educational purposes. Sobol numbers can be used. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividendCurve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _valueMC_NUMBA_ONLY(stock_price,
                           texp, 
                           self._strikePrice,
                           self._optionType.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUMBER", self._numOptions, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
