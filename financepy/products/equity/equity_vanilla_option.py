##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import List, Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError

from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_shapes

from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import option_years

from ...models.model import Model
from ...models.black_scholes import BlackScholes
from ...models.black_scholes_analytic import european_value
from ...models.black_scholes_analytic import delta
from ...models.black_scholes_analytic import vega
from ...models.black_scholes_analytic import gamma
from ...models.black_scholes_analytic import rho
from ...models.black_scholes_analytic import vanna
from ...models.black_scholes_analytic import theta
from ...models.black_scholes_analytic import implied_volatility
from ...models.black_scholes_analytic import intrinsic

from ...models.black_scholes_mc import value_mc_nonumba_nonumpy
from ...models.black_scholes_mc import value_mc_numpy_numba
from ...models.black_scholes_mc import value_mc_numba_only
from ...models.black_scholes_mc import value_mc_numpy_only
from ...models.black_scholes_mc import value_mc_numba_parallel

########################################################################################


class EquityVanillaOption:
    """Class for managing plain vanilla European calls and puts on equities.
    For American calls and puts see the EquityAmericanOption class."""

    def __init__(
        self,
        expiry_dt: Date | List,
        strike_price: Union[float, np.ndarray],
        opt_type: OptionTypes,
        num_options: float = 1.0,
    ):
        """Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options."""

        check_argument_types(self.__init__, locals())

        self.opt_type = opt_type
        self.opt_type_value = opt_type.value
        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.num_options = num_options

    ###########################################################################

    def intrinsic(
        self,
        value_dt: Date,
        stock_price: Union[np.ndarray, float],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
    ):
        """Calculate the discounted intrinsic value of the option."""

        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price
        s0 = stock_price

        check_shapes(s0, k, t_exp)

        intrinsic_value = intrinsic(s0, t_exp, k, r, q, self.opt_type_value)
        intrinsic_value = intrinsic_value * self.num_options
        return intrinsic_value

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float | np.ndarray,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Equity Vanilla Option valuation using Black-Scholes model."""

        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price
        s0 = stock_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):

            v = model.volatility
            value = european_value(s0, t_exp, k, r, q, v, self.opt_type_value)

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):

            v = model.volatility
            d = delta(s0, t_exp, k, r, q, v, self.opt_type_value)

        else:
            raise FinError("Unknown Model Type")

        return d

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):

            v = model.volatility
            g = gamma(s0, t_exp, k, r, q, v, self.opt_type_value)

        else:
            raise FinError("Unknown Model Type")

        return g

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):

            v = model.volatility
            veg = vega(s0, t_exp, k, r, q, v, self.opt_type_value)

        else:
            raise FinError("Unknown Model Type")

        return veg

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):
            v = model.volatility
            thet = theta(s0, t_exp, k, r, q, v, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")

        return thet

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):
            v = model.volatility
            rh = rho(s0, t_exp, k, r, q, v, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")

        return rh

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):
            v = model.volatility
            van = vanna(s0, t_exp, k, r, q, v, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")

        return van

    ###########################################################################

    def implied_volatility(
        self,
        value_dt: Date,
        stock_price: Union[float, list, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        price,
    ):
        """Calculate the Black-Scholes implied volatility of a European
        vanilla option."""

        t_exp = option_years(value_dt, self.expiry_dt)

        if t_exp < 1.0 / 366.0:
            raise FinError("Expiry time is too close to zero.")

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_price

        if np.ndim(k) != 0:
            raise FinError("Strike price must be scalar for implied volatility.")

        sigma = implied_volatility(s0, t_exp, k, r, q, price, self.opt_type_value)

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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_numpy_only(
            s0,
            t_exp,
            self.strike_price,
            self.opt_type_value,
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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_numba_only(
            s0,
            t_exp,
            self.strike_price,
            self.opt_type_value,
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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_numba_parallel(
            s0,
            t_exp,
            self.strike_price,
            self.opt_type_value,
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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_numpy_numba(
            s0,
            t_exp,
            self.strike_price,
            self.opt_type_value,
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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_nonumba_nonumpy(
            s0,
            t_exp,
            self.strike_price,
            r,
            q,
            vol,
            self.opt_type_value,
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

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        s0 = stock_price

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        vol = model.volatility

        v = value_mc_numba_only(
            s0,
            t_exp,
            self.strike_price,
            r,
            q,
            vol,
            self.opt_type_value,
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
        s += label_to_string("OPTION TYPE VALUE", self.opt_type)
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
