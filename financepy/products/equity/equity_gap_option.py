########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import BarrierTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...models.process_simulator import GBMNumericalSchemeTypes
from ...models.process_simulator import ProcessTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...models.bs_barrier_option import barrier_option_value
from ...models.bs_barrier_option_mc import value_barrier_option_mc
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_strike_price
from ...utils.helpers import option_years
from ...models.model import Model
from ...utils.global_types import OptionTypes

# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

########################################################################################


class EquityGapOption(EquityOption):
    """Class to hold details of an Equity Gap Option. It also
    calculates the option price using Black Scholes."""

    def __init__(
        self,
        expiry_dt: Date,
        trigger_price: float,
        strike_price: float,
        opt_type: OptionTypes,
        num_options: float = 1.0,
    ) -> None:
        """Create the EquityBarrierOption by specifying the expiry date,
        strike price, OPTION_TYPE, barrier level, the number of observations
        per year and the notional."""

        check_argument_types(self.__init__, locals())

        check_strike_price(strike_price)

        self.expiry_dt = expiry_dt
        self.trigger_price = float(trigger_price)
        self.strike_price = float(strike_price)
        self.opt_type = opt_type
        self.opt_type_value = opt_type.value

        self.num_options = num_options

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Value a European equity gap option using Black-Scholes.

        The payoff at expiry is:

            Call: (S_T - K2) * I(S_T > K1)
            Put:  (K2 - S_T) * I(S_T < K1)

        where:
            K1 = trigger_price
            K2 = strike_price
        """

        from scipy.stats import norm

        t_exp = option_years(value_dt, self.expiry_dt)

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        if t_exp <= 0.0:
            raise FinError("Expiry date must be after valuation date")

        s = np.asarray(stock_price, dtype=float)
        k1 = self.trigger_price
        k2 = self.strike_price
        sigma = model.volatility

        df = discount_curve.df(self.expiry_dt)
        dq = dividend_curve.df(self.expiry_dt)

        # exp((r-q)T) = dq / df
        sqrt_t = np.sqrt(t_exp)

        d1 = (
            np.log(s / k1)
            + np.log(dq / df)
            + 0.5 * sigma * sigma * t_exp
        ) / (sigma * sqrt_t)

        d2 = d1 - sigma * sqrt_t

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            value = (
                s * dq * norm.cdf(d1)
                - k2 * df * norm.cdf(d2)
            )

        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            value = (
                k2 * df * norm.cdf(-d2)
                - s * dq * norm.cdf(-d1)
            )

        else:
            raise FinError("Unknown option type")

        value *= self.num_options

        if np.ndim(stock_price) == 0:
            return float(value)

        return value

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 42,
    ):
        """Value a European equity gap option using Monte Carlo.

        The payoff at expiry is:

            Call: (S_T - K2) * I(S_T > K1)
            Put:  (K2 - S_T) * I(S_T < K1)

        where:
            K1 = trigger_price
            K2 = strike_price
        """

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date")

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if t_exp <= 0.0:
            raise FinError("Expiry date must be after valuation date")

        if num_paths <= 0:
            raise FinError("Number of paths must be positive")

        s0 = float(stock_price)
        k1 = self.trigger_price
        k2 = self.strike_price
        sigma = model.volatility

        df = discount_curve.df(self.expiry_dt)
        dq = dividend_curve.df(self.expiry_dt)

        # Continuously compounded rates implied by the discount factors
        r = -np.log(df) / t_exp
        q = -np.log(dq) / t_exp

        # Use antithetic variates for variance reduction.
        rng = np.random.default_rng(seed)

        num_half_paths = (num_paths + 1) // 2
        z = rng.standard_normal(num_half_paths)

        z = np.concatenate((z, -z))[:num_paths]

        # Exact terminal GBM simulation under the risk-neutral measure
        st = s0 * np.exp(
            (r - q - 0.5 * sigma * sigma) * t_exp
            + sigma * np.sqrt(t_exp) * z
        )

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            payoff = np.where(
                st > k1,
                st - k2,
                0.0,
            )

        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            payoff = np.where(
                st < k1,
                k2 - st,
                0.0,
            )

        else:
            raise FinError("Unknown option type")

        value = df * np.mean(payoff)

        return value * self.num_options

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("TRIGGER PRICE", self.trigger_price)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("NUM OPTIONS", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
