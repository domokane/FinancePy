###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

from ...models.equity_barrier_models import value_barrier
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...models.process_simulator import FinProcessSimulator

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import EquityBarrierTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.global_vars import g_days_in_year


# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

###############################################################################


class EquityBarrierOption(EquityOption):
    """Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum EquityBarrierTypes."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        opt_type: EquityBarrierTypes,
        barrier_level: float,
        num_obs_per_year: (int, float) = 252,
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

        if opt_type not in EquityBarrierTypes:
            raise FinError("Option Type " + str(opt_type) + " unknown.")

        self.opt_type = opt_type
        self.notional = notional

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: (float, np.ndarray),
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

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        if t_exp < 0:
            raise FinError("Option expires before value date.")

        values = value_barrier(
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

    ###########################################################################

    def value_mc(
        self,
        t: float,
        k,
        opt_type: int,
        b,
        notional,
        s: float,
        r: float,
        process_type,
        model_params,
        num_ann_obs: int = 252,
        num_paths: int = 10000,
        seed: int = 4242,
    ):
        """A Monte-Carlo based valuation of the barrier option which simulates
        the evolution of the stock price of at a specified number of annual
        observation times until expiry to examine if the barrier has been
        crossed and the corresponding value of the final payoff, if any. It
        assumes a GBM model for the stock price."""
        t = max(t, 1e-6)
        num_time_steps = int(t * num_ann_obs)
        process = FinProcessSimulator()

        #######################################################################

        if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and s <= b:
            return 0.0
        elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and s >= b:
            return 0.0
        elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and s <= b:
            return 0.0
        elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and s >= b:
            return 0.0

        #######################################################################

        simple_call = False
        simple_put = False

        if opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and s <= b:
            simple_call = True
        elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value and s >= b:
            simple_call = True
        elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value and s >= b:
            simple_put = True
        elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and s <= b:
            simple_put = True

        if simple_put or simple_call:
            s_all = process.get_process(
                process_type, t, model_params, 1, num_paths, seed
            )

        if simple_call:
            c = (np.maximum(s_all[:, -1] - k, 0.0)).mean()
            c = c * np.exp(-r * t)
            return c

        if simple_put:
            p = (np.maximum(k - s_all[:, -1], 0.0)).mean()
            p = p * np.exp(-r * t)
            return p

        # Get full set of paths
        s_all = process.get_process(
            process_type, t, model_params, num_time_steps, num_paths, seed
        )

        (num_paths, num_time_steps) = s_all.shape

        if (
            opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value
            or opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value
            or opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value
            or opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value
        ):

            barrier_crossed_from_above = [False] * num_paths

            for p in range(0, num_paths):
                barrier_crossed_from_above[p] = np.any(s_all[p] <= b)

        if (
            opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value
            or opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value
            or opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value
            or opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value
        ):

            barrier_crossed_from_below = [False] * num_paths
            for p in range(0, num_paths):
                barrier_crossed_from_below[p] = np.any(s_all[p] >= b)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * (
                ones - barrier_crossed_from_above
            )
        elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
            payoff = (
                np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_above
            )
        elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
            payoff = (
                np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_below
            )
        elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * (
                ones - barrier_crossed_from_below
            )
        elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
            payoff = (
                np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_below
            )
        elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * (
                ones - barrier_crossed_from_below
            )
        elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * (
                ones - barrier_crossed_from_above
            )
        elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
            payoff = (
                np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_above
            )
        else:
            raise FinError("Unknown barrier option type." + str(opt_type))

        v = payoff.mean() * np.exp(-r * t)

        return v * notional

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


###############################################################################
