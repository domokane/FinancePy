########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from typing import Union

import numpy as np

from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...models.equity_barrier_option_bs import value_equity_barrier_option_bs
from ...models.equity_barrier_option_mc import value_equity_barrier_option_mc
from ...models.process_simulator import FinGBMNumericalScheme
from ...models.process_simulator import ProcessTypes

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import BarrierTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.global_vars import G_DAYS_IN_YEAR


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
        opt_type: BarrierTypes,
        barrier_level: float,
        num_obs_per_year: Union[int, float] = 252,
        notional: float = 1.0,
    ):
        """Create the EquityBarrierOption by specifying the expiry date,
        strike price, option type, barrier level, the number of observations
        per year and the notional."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.barrier_level = float(barrier_level)
        self.num_obs_per_year = int(num_obs_per_year)

        if opt_type not in BarrierTypes:
            raise FinError("Option Type " + str(opt_type) + " unknown.")

        self.opt_type = opt_type
        self.notional = notional

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        """

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        if isinstance(stock_price, int):
            stock_price = float(stock_price)

        if isinstance(stock_price, float):
            stock_prices = [stock_price]
        else:
            stock_prices = stock_price

        values = []

        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEAR

        if t_exp < 0:
            raise FinError("Option expires before value date.")

        values = value_equity_barrier_option_bs(
            t_exp,
            self.strike_price,
            self.barrier_level,
            stock_prices,
            discount_curve.cc_rate(self.expiry_dt),
            dividend_curve.cc_rate(self.expiry_dt),
            model.volatility,
            self.opt_type.value,
            self.num_obs_per_year,
        )

        values = values * self.notional

        if isinstance(stock_price, float):
            return values[0]
        else:
            return np.array(values)

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_obs_per_year=252,
        num_paths: int = 10000,
        seed: int = 42,
    ):
        """This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        """

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if isinstance(stock_price, int):
            stock_price = float(stock_price)

        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEAR

        if t_exp < 0:
            raise FinError("Option expires before value date.")

        risk_free_rate = discount_curve.cc_rate(self.expiry_dt)
        dividend_rate = dividend_curve.cc_rate(self.expiry_dt)
        drift = risk_free_rate - dividend_rate

        scheme = FinGBMNumericalScheme.NORMAL_SCHEME

        model_params = (stock_price, drift, model.volatility, scheme)

        process_type = ProcessTypes.GBM_PROCESS

        value = value_equity_barrier_option_mc(
            t_exp,
            self.strike_price,
            self.opt_type,
            self.barrier_level,
            self.notional,
            stock_price,
            risk_free_rate,
            process_type,
            model_params,
            num_obs_per_year,
            num_paths,
            seed,
        )

        value = value * self.notional

        #        if isinstance(stock_price, float):
        #            return values[0]
        #        else:
        #            return np.array(values)

        return value

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_level)
        s += label_to_string("NUM OBSERVATIONS", self.num_obs_per_year)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
