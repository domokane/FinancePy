########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import BarrierTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.global_types import GBMNumericalSchemeTypes
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...models.equity_barrier_option_bs import value_equity_barrier_option_bs
from ...models.barrier_option_mc import value_barrier_option_mc
from ...models.process_simulator import ProcessTypes
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_strike_price
from ...utils.helpers import option_years

# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

########################################################################################


class EquityBarrierOption(EquityOption):
    """Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum BarrierTypes."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        barrier_type: BarrierTypes,
        barrier_level: float,
        num_obs_per_year: Union[int, float] = 252,
        notional: float = 1.0,
    ) -> None:
        """Create the EquityBarrierOption by specifying the expiry date,
        strike price, OPTION_TYPE, barrier level, the number of observations
        per year and the notional."""

        check_argument_types(self.__init__, locals())

        check_strike_price(strike_price)

        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.barrier_level = float(barrier_level)
        self.num_obs_per_year = int(num_obs_per_year)

        if barrier_type not in BarrierTypes:
            raise FinError("OPTION_TYPE " + str(barrier_type) + " unknown.")

        self.barrier_type = barrier_type
        self.notional = notional

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        """

        t_exp = option_years(value_dt, self.expiry_dt)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        values = []

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        values = value_equity_barrier_option_bs(
            t_exp,
            self.strike_price,
            self.barrier_level,
            stock_price,
            r,
            q,
            model.volatility,
            self.barrier_type.value,
            self.num_obs_per_year,
        )

        values = values * self.notional

        if isinstance(stock_price, float):
            return values
        else:
            return np.array(values)

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float | np.ndarray,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_obs_per_year: int=252,
        num_paths: int = 10000,
        seed: int = 42,
    ):
        """Value the Equity Barrier Option using Monte Carlo."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        check_stock_price(stock_price)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        mu = r - q

        scheme = GBMNumericalSchemeTypes.ANTITHETIC

        model_params = (stock_price, mu, model.volatility, scheme)

        process_type = ProcessTypes.GBM_PROCESS

        value = value_barrier_option_mc(
            t_exp,
            self.strike_price,
            self.barrier_type,
            self.barrier_level,
            stock_price,
            r,
            process_type,
            model_params,
            num_obs_per_year,
            num_paths,
            seed,
        )

        value = value * self.notional
        return value

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("BARRIER TYPE", self.barrier_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_level)
        s += label_to_string("NUM OBSERVATIONS", self.num_obs_per_year)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
