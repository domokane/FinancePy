##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit

# from scipy import optimize
from ...utils.solver_1d import newton_secant, bisection, newton

from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
from ...utils.global_types import FinOptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.discount.curve import DiscountCurve

from ...models.FinModel import FinModel
from ...models.black_scholes import BlackScholes
from ...models.black_scholes_analytic import bs_value
from ...models.black_scholes_analytic import bs_delta
from ...models.black_scholes_analytic import bs_vega
from ...models.black_scholes_analytic import bs_gamma
from ...models.black_scholes_analytic import bs_rho
from ...models.black_scholes_analytic import bs_theta
from ...models.black_scholes_analytic import bsImpliedVolatility
from ...models.black_scholes_analytic import bsIntrinsic


from ...models.black_scholes_mc import _value_mc_NONUMBA_NONUMPY
from ...models.black_scholes_mc import _value_mc_NUMPY_NUMBA
from ...models.black_scholes_mc import _value_mc_NUMBA_ONLY
from ...models.black_scholes_mc import _value_mc_NUMPY_ONLY
from ...models.black_scholes_mc import _value_mc_NUMBA_PARALLEL

from ...utils.math import N

###############################################################################

@njit(fastmath=True, cache=True)
def _f(v, args):

    option_type_value = int(args[0])
    texp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]
    price = args[6]

    obj_fn = bs_value(s0, texp, k, r, q, v, option_type_value)
    obj_fn = obj_fn - price
    return obj_fn

###############################################################################


def _fvega(v, *args):

    self = args[0]
    texp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]

    fprime = bs_vega(s0, texp, k, r, q, v, self._option_type.value)
    return fprime

###############################################################################


class EquityVanillaOption():
    """ Class for managing plain vanilla European calls and puts on equities.
    For American calls and puts see the EquityAmericanOption class. """

    def __init__(self,
                 expiry_date: (Date, list),
                 strike_price: (float, np.ndarray),
                 option_type: (FinOptionTypes, list),
                 numOptions: float = 1.0):
        """ Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options. """

        check_argument_types(self.__init__, locals())

        if option_type != FinOptionTypes.EUROPEAN_CALL and \
           option_type != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(option_type))

        self._expiry_date = expiry_date
        self._strike_price = strike_price
        self._option_type = option_type
        self._num_options = numOptions
        self._texp = None

###############################################################################

    def intrinsic(self,
                  valuation_date: (Date, list),
                  stock_price: (np.ndarray, float),
                  discount_curve: DiscountCurve,
                  dividend_curve: DiscountCurve):
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        intrinsicValue = bsIntrinsic(s0, texp, k, r, q,
                                     self._option_type.value)

        intrinsicValue = intrinsicValue * self._num_options
        return intrinsicValue

###############################################################################

    def value(self,
              valuation_date: (Date, list),
              stock_price: (np.ndarray, float),
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            value = bs_value(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        value = value * self._num_options
        return value

###############################################################################

    def delta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            delta = bs_delta(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        return delta

###############################################################################

    def gamma(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            gamma = bs_gamma(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        return gamma

###############################################################################

    def vega(self,
             valuation_date: Date,
             stock_price: float,
             discount_curve: DiscountCurve,
             dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            vega = bs_vega(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        return vega

###############################################################################

    def theta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            theta = bs_theta(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        return theta

###############################################################################

    def rho(self,
            valuation_date: Date,
            stock_price: float,
            discount_curve: DiscountCurve,
            dividend_curve: DiscountCurve,
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

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price

        if isinstance(model, BlackScholes):

            v = model._volatility
            rho = bs_rho(s0, texp, k, r, q, v, self._option_type.value)

        else:
            raise FinError("Unknown Model Type")

        return rho

###############################################################################

    def implied_volatility(self,
                          valuation_date: Date,
                          stock_price: (float, list, np.ndarray),
                          discount_curve: DiscountCurve,
                          dividend_curve: DiscountCurve,
                          price):
        """ Calculate the Black-Scholes implied volatility of a European 
        vanilla option. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        if texp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        k = self._strike_price
        s0 = stock_price

        sigma = bsImpliedVolatility(s0, texp, k, r, q, price, 
                                    self._option_type.value)
        
        return sigma

###############################################################################

    def value_mc_NUMPY_ONLY(self,
                           valuation_date: Date,
                           stock_price: float,
                           discount_curve: DiscountCurve,
                           dividend_curve: DiscountCurve,
                           model:FinModel,
                           num_paths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NUMPY_ONLY(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def value_mc_NUMBA_ONLY(self,
                           valuation_date: Date,
                           stock_price: float,
                           discount_curve: DiscountCurve,
                           dividend_curve: DiscountCurve,
                           model:FinModel,
                           num_paths: int = 10000,
                           seed: int = 4242,
                           useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NUMBA_ONLY(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

        return v

###############################################################################

    def value_mc_NUMBA_PARALLEL(self,
                               valuation_date: Date,
                               stock_price: float,
                               discount_curve: DiscountCurve,
                               dividend_curve: DiscountCurve,
                               model:FinModel,
                               num_paths: int = 10000,
                               seed: int = 4242,
                               useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NUMBA_PARALLEL(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

#        _value_mc_NUMBA_ONLY.parallel_diagnostics(level=4)

        return v

###############################################################################

    def value_mc_NUMPY_NUMBA(self,
                            valuation_date: Date,
                            stock_price: float,
                            discount_curve: DiscountCurve,
                            dividend_curve: DiscountCurve,
                            model:FinModel,
                            num_paths: int = 10000,
                            seed: int = 4242,
                            useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NUMPY_NUMBA(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def value_mc_NONUMBA_NONUMPY(self,
                                valuation_date: Date,
                                stock_price: float,
                                discount_curve: DiscountCurve,
                                dividend_curve: DiscountCurve,
                                model:FinModel,
                                num_paths: int = 10000,
                                seed: int = 4242,
                                useSobol: int = 0):

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NONUMBA_NONUMPY(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed,
                           useSobol)

        return v

###############################################################################

    def value_mc(self,
                valuation_date: Date,
                stock_price: float,
                discount_curve: DiscountCurve,
                dividend_curve: DiscountCurve,
                model:FinModel,
                num_paths: int = 10000,
                seed: int = 4242,
                useSobol: int = 0):
        """ Value European style call or put option using Monte Carlo. This is
        mainly for educational purposes. Sobol numbers can be used. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        df = discount_curve.df(self._expiry_date)
        r = -np.log(df)/texp

        dq = dividend_curve.df(self._expiry_date)
        q = -np.log(dq)/texp

        vol = model._volatility

        v = _value_mc_NUMBA_ONLY(stock_price,
                           texp, 
                           self._strike_price,
                           self._option_type.value,
                           r, 
                           q, 
                           vol, 
                           num_paths,
                           seed, 
                           useSobol)

        return v

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("NUMBER", self._num_options, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
