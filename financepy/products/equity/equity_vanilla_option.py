##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit

# from scipy import optimize
from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve

from ...models.model import Model
from ...models.black_scholes import BlackScholes
from ...models.black_scholes_analytic import bs_value
from ...models.black_scholes_analytic import bs_delta
from ...models.black_scholes_analytic import bs_vega
from ...models.black_scholes_analytic import bs_gamma
from ...models.black_scholes_analytic import bs_rho
from ...models.black_scholes_analytic import bs_vanna
from ...models.black_scholes_analytic import bs_theta
from ...models.black_scholes_analytic import bs_implied_volatility
from ...models.black_scholes_analytic import bs_intrinsic

from ...models.black_scholes_mc import _value_mc_nonumba_nonumpy
from ...models.black_scholes_mc import _value_mc_numpy_numba
from ...models.black_scholes_mc import _value_mc_numba_only
from ...models.black_scholes_mc import _value_mc_numpy_only
from ...models.black_scholes_mc import _value_mc_numba_parallel

###############################################################################


@njit(fastmath=True, cache=True)
def _f(v, args):

    option_type_value = int(args[0])
    t_exp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]
    price = args[6]

    obj_fn = bs_value(s0, t_exp, k, r, q, v, option_type_value)
    obj_fn = obj_fn - price
    return obj_fn


###############################################################################


def _fvega(v, *args):

    self = args[0]
    t_exp = args[1]
    s0 = args[2]
    r = args[3]
    q = args[4]
    k = args[5]

    fprime = bs_vega(s0, t_exp, k, r, q, v, self.option_type.value)
    return fprime


###############################################################################


class EquityVanillaOption:
    """Class for managing plain vanilla European calls and puts on equities.
    For American calls and puts see the EquityAmericanOption class."""

    def __init__(
        self,
        expiry_dt: (Date, list),
        strike_price: (float, np.ndarray),
        option_type: (OptionTypes, list),
        num_options: float = 1.0,
    ):
        """Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options."""

        check_argument_types(self.__init__, locals())

        if isinstance(option_type, OptionTypes):
            option_type_value = option_type.value
        elif isinstance(option_type, list):
            option_type_value = []
            for opt in option_type:
                option_type_value.append(opt.value)
            option_type_value = np.array(option_type_value)

        self.option_type_value = option_type_value

        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.option_type = option_type
        self.num_options = num_options
        self.t_exp = None

    ###########################################################################

    def intrinsic(
        self,
        value_dt: (Date, list),
        stock_price: (np.ndarray, float),
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
    ):
        """Equity Vanilla Option valuation using Black-Scholes model."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if isinstance(self.expiry_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        elif isinstance(self.expiry_dt, list):
            t_exp = []
            for exp_dt in self.expiry_dt:
                t = (exp_dt - value_dt) / g_days_in_year
            t_exp.append(t)
            t_exp = np.array(t_exp)
        else:
            t_exp = value_dt

        self.t_exp = t_exp

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        intrinsic_value = bs_intrinsic(
            s0, t_exp, k, r, q, self.option_type_value
        )

        intrinsic_value = intrinsic_value * self.num_options
        return intrinsic_value

    ###########################################################################

    def value(
        self,
        value_dt: (Date, list),
        stock_price: (np.ndarray, float),
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Equity Vanilla Option valuation using Black-Scholes model."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if isinstance(self.expiry_dt, list):
            if any(value_dt > self.expiry_dt):
                raise FinError("Valuation date after expiry dates.")
        elif value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        if isinstance(self.expiry_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        elif isinstance(self.expiry_dt, list):
            t_exp = []
            for exp_dt in self.expiry_dt:
                t = (exp_dt - value_dt) / g_days_in_year
            t_exp.append(t)
            t_exp = np.array(t_exp)
        else:
            t_exp = value_dt

        self.t_exp = t_exp

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price

        t_exp = np.maximum(t_exp, 1e-10)

        # Extract the discount. Adjust if tvalue date is not same as curve date
        # I decided to put an error message - may reconsider
        df_expiry = discount_curve.df(self.expiry_dt)
        # df_value = discount_curve.df(value_dt)
        # df = df_expiry / df_value
        r = -np.log(df_expiry) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):

            v = model.volatility
            value = bs_value(s0, t_exp, k, r, q, v, self.option_type_value)

        else:
            raise FinError("Unknown Model Type")

        value = value * self.num_options
        return value

    ###########################################################################

    def delta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculate the analytical delta of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        self.t_exp = t_exp

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):

            v = model.volatility
            delta = bs_delta(s0, t_exp, k, r, q, v, self.option_type_value)

        else:
            raise FinError("Unknown Model Type")

        return delta

    ###########################################################################

    def gamma(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate the analytical gamma of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price

        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):

            v = model.volatility
            gamma = bs_gamma(s0, t_exp, k, r, q, v, self.option_type_value)

        else:
            raise FinError("Unknown Model Type")

        return gamma

    ###########################################################################

    def vega(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate the analytical vega of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):

            v = model.volatility
            vega = bs_vega(s0, t_exp, k, r, q, v, self.option_type_value)

        else:
            raise FinError("Unknown Model Type")

        return vega

    ###########################################################################

    def theta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate the analytical theta of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):
            v = model.volatility
            theta = bs_theta(s0, t_exp, k, r, q, v, self.option_type_value)
        else:
            raise FinError("Unknown Model Type")

        return theta

    ###########################################################################

    def rho(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate the analytical rho of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):
            v = model.volatility
            rho = bs_rho(s0, t_exp, k, r, q, v, self.option_type_value)
        else:
            raise FinError("Unknown Model Type")

        return rho

    ###########################################################################

    def vanna(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate the analytical vanna of a European vanilla option."""

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        s0 = stock_price
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price

        if isinstance(model, BlackScholes):
            v = model.volatility
            vanna = bs_vanna(s0, t_exp, k, r, q, v, self.option_type_value)
        else:
            raise FinError("Unknown Model Type")

        return vanna

    ###########################################################################

    def implied_volatility(
        self,
        value_dt: Date,
        stock_price: (float, list, np.ndarray),
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        price,
    ):
        """Calculate the Black-Scholes implied volatility of a European
        vanilla option."""

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        if t_exp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        k = self.strike_price
        s0 = stock_price

        sigma = bs_implied_volatility(
            s0, t_exp, k, r, q, price, self.option_type_value
        )

        return sigma

    ###########################################################################

    def value_mc_numpy_only(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_numpy_only(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type.value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        return v

    ###########################################################################

    def value_mc_numba_only(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_numba_only(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type_value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        return v

    ###########################################################################

    def value_mc_numba_parallel(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_numba_parallel(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type_value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        #        _value_mc_NUMBA_ONLY.parallel_diagnostics(level=4)

        return v

    ###########################################################################

    def value_mc_numpy_numba(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_numpy_numba(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type_value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        return v

    ###########################################################################

    def value_mc_nonumba_nonumpy(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_nonumba_nonumpy(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type.value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
        use_sobol: int = 0,
    ):
        """Value European style call or put option using Monte Carlo. This is
        mainly for educational purposes. Sobol numbers can be used."""

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t_exp

        vol = model.volatility

        v = _value_mc_numba_only(
            stock_price,
            t_exp,
            self.strike_price,
            self.option_type_value,
            r,
            q,
            vol,
            num_paths,
            seed,
            use_sobol,
        )

        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE VALUE", self.option_type)
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
