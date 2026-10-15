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
from ...utils.mc_result import MCResult

from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import option_years

from ...models.model import Model
from ...models.black_scholes import BlackScholes
from ...models.bs_power_option import power_option_value
from ...models.bs_power_option import power_option_delta

from .equity_option import EquityOption

########################################################################################


class EquityPowerOption(EquityOption):
    """Class for managing European power option calls and puts on equities."""

    def __init__(
        self,
        expiry_dt: Date | List,
        strike_price: Union[float, np.ndarray],
        opt_type: OptionTypes,
        power: float = 1.0,
        num_options: float = 1.0,
    ) -> None:
        """Create the Equity Vanilla option object by specifying the expiry
        date, the option strike, the option type and the number of options."""

        check_argument_types(self.__init__, locals())

        if not np.isfinite(power) or power <= 0.0:
            raise FinError("Power must be strictly positive.")

        self.opt_type = opt_type
        self.opt_type_value = opt_type.value
        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.power = power
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
        n = self.power

        check_shapes(s0, k, t_exp)

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            phi = 1.0
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            phi = -1.0
        else:
            raise FinError("Unknown option type")

        sT = s0 * np.exp((r - q) * t_exp)

        intrinsic_value = (
            np.exp(-r * t_exp)
            * np.maximum(phi * (sT**n - k), 0.0)
        )

        return intrinsic_value * self.num_options

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
        n = self.power

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):

            v = model.volatility
            value = power_option_value(s0, t_exp, k, r, q, v, self.opt_type_value, n)

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
        model: Model,
    ):
        """Analytical spot delta for a European power call or put."""
        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_shapes(stock_price, self.strike_price, t_exp)

        if not isinstance(model, BlackScholes):
            raise FinError("Unknown Model Type")

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        return self.num_options * power_option_delta(
            stock_price, t_exp, self.strike_price,
            r, q, model.volatility, self.opt_type_value, self.power,
        )

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
    ) -> MCResult:
        """Value the power option using Monte Carlo simulation.

        The terminal stock price is simulated under the risk-neutral measure:

            S(T) = S(0) * exp(
                (r - q - 0.5 * sigma**2) * T
                + sigma * sqrt(T) * Z
            )

        and the payoff is:

            Call: max(S(T)**n - K, 0)
            Put:  max(K - S(T)**n, 0)

        Parameters
        ----------
        value_dt : Date
            Valuation date.
        stock_price : float
            Current stock price.
        discount_curve : DiscountCurve
            Risk-free discount curve.
        dividend_curve : DiscountCurve
            Dividend discount curve.
        model : Model
            Black-Scholes model.
        num_paths : int
            Number of Monte Carlo paths.
        seed : int
            Random number seed.

        Returns
        -------
        MCResult
            Monte Carlo value and standard error.
        """

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if not isinstance(model, BlackScholes):
            raise FinError("Unknown Model Type")

        if self.power <= 0:
            raise FinError("Power must be positive")

        if num_paths <= 0:
            raise FinError("Number of paths must be positive")

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        s0 = stock_price
        k = self.strike_price
        n = self.power
        v = model.volatility

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            phi = 1.0
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            phi = -1.0
        else:
            raise FinError("Unknown option type")

        rng = np.random.default_rng(seed)

        num_pairs = (num_paths + 1) // 2
        z = rng.standard_normal(num_pairs)

        z = np.concatenate((z, -z))
        z = z[:num_paths]

        mu = r - q - 0.5 * v * v
        s_t = s0 * np.exp(mu * t_exp + v * np.sqrt(t_exp) * z)

        payoff = np.maximum(phi * (s_t**n - k), 0.0)

        discounted_payoff = np.exp(-r * t_exp) * payoff
        discounted_payoff *= self.num_options

        value = np.mean(discounted_payoff)

        std_error = (
            np.std(discounted_payoff, ddof=1)
            / np.sqrt(num_paths)
        )

        return MCResult(value, std_error)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE VALUE", self.opt_type)
        s += label_to_string("POWER", self.power, "")
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
