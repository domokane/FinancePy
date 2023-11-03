###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

from ...models.equity_barrier_models import value_bs
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...models.process_simulator import FinProcessSimulator

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import EquityBarrierTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.global_vars import gDaysInYear


# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

###############################################################################


class EquityBarrierOption(EquityOption):
    """ Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum EquityBarrierTypes. """

    def __init__(self,
                 expiry_date: Date,
                 strike_price: float,
                 option_type: EquityBarrierTypes,
                 barrier_level: float,
                 num_obs_per_year: (int, float) = 252,
                 notional: float = 1.0):
        """ Create the EquityBarrierOption by specifying the expiry date,
        strike price, option type, barrier level, the number of observations
        per year and the notional. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strike_price = float(strike_price)
        self._barrier_level = float(barrier_level)
        self._num_obs_per_year = int(num_obs_per_year)

        if option_type not in EquityBarrierTypes:
            raise FinError("Option Type " + str(option_type) + " unknown.")

        self._option_type = option_type
        self._notional = notional

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: (float, np.ndarray),
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        """

        if isinstance(valuation_date, Date) is False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option value date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option value date")

        if isinstance(stock_price, int):
            stock_price = float(stock_price)

        if isinstance(stock_price, float):
            stock_prices = [stock_price]
        else:
            stock_prices = stock_price

        values = []

        t_exp = (self._expiry_date - valuation_date) / gDaysInYear

        if t_exp < 0:
            raise FinError("Option expires before value date.")

        values = value_bs(t_exp,
                          self._strike_price,
                          self._barrier_level,
                          stock_prices,
                          discount_curve.cc_rate(self._expiry_date),
                          dividend_curve.cc_rate(self._expiry_date),
                          model._volatility,
                          self._option_type.value,
                          self._num_obs_per_year)

        values = values * self._notional

        if isinstance(stock_price, float):
            return values[0]
        else:
            return np.array(values)

###############################################################################

    def value_mc(self,
                 t: float,
                 k,
                 option_type: int,
                 b,
                 notional,
                 s: float,
                 r: float,
                 process_type,
                 model_params,
                 num_ann_obs: int = 252,
                 num_paths: int = 10000,
                 seed: int = 4242):
        """ A Monte-Carlo based valuation of the barrier option which simulates
        the evolution of the stock price of at a specified number of annual
        observation times until expiry to examine if the barrier has been
        crossed and the corresponding value of the final payoff, if any. It
        assumes a GBM model for the stock price. """
        t = max(t, 1e-6)
        num_time_steps = int(t * num_ann_obs)
        option_type = option_type
        process = FinProcessSimulator()

        #######################################################################

        if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and s <= b:
            return 0.0
        elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and s >= b:
            return 0.0
        elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and s <= b:
            return 0.0
        elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and s >= b:
            return 0.0

        #######################################################################

        simple_call = False
        simple_put = False

        if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and s <= b:
            simple_call = True
        elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value and s >= b:
            simple_call = True
        elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value and s >= b:
            simple_put = True
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and s <= b:
            simple_put = True

        if simple_put or simple_call:
            Sall = process.get_process(
                process_type, t, model_params, 1, num_paths, seed)

        if simple_call:
            c = (np.maximum(Sall[:, -1] - k, 0.0)).mean()
            c = c * np.exp(-r * t)
            return c

        if simple_put:
            p = (np.maximum(k - Sall[:, -1], 0.0)).mean()
            p = p * np.exp(-r * t)
            return p

        # Get full set of paths
        Sall = process.get_process(process_type, t, model_params, num_time_steps,
                                   num_paths, seed)

        (num_paths, num_time_steps) = Sall.shape

        if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value or \
                option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value or \
                option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value or \
                option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:

            barrier_crossed_from_above = [False] * num_paths

            for p in range(0, num_paths):
                barrier_crossed_from_above[p] = np.any(Sall[p] <= b)

        if option_type == EquityBarrierTypes.UP_AND_IN_CALL.value or \
                option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value or \
                option_type == EquityBarrierTypes.UP_AND_IN_PUT.value or \
                option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:

            barrier_crossed_from_below = [False] * num_paths
            for p in range(0, num_paths):
                barrier_crossed_from_below[p] = np.any(Sall[p] >= b)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
            payoff = np.maximum(Sall[:, -1] - k, 0.0) * \
                     (ones - barrier_crossed_from_above)
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
            payoff = np.maximum(Sall[:, -1] - k, 0.0) * \
                barrier_crossed_from_above
        elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
            payoff = np.maximum(Sall[:, -1] - k, 0.0) * \
                barrier_crossed_from_below
        elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
            payoff = np.maximum(Sall[:, -1] - k, 0.0) * \
                     (ones - barrier_crossed_from_below)
        elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
            payoff = np.maximum(k - Sall[:, -1], 0.0) * \
                barrier_crossed_from_below
        elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
            payoff = np.maximum(k - Sall[:, -1], 0.0) * \
                     (ones - barrier_crossed_from_below)
        elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
            payoff = np.maximum(k - Sall[:, -1], 0.0) * \
                     (ones - barrier_crossed_from_above)
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
            payoff = np.maximum(k - Sall[:, -1], 0.0) * \
                barrier_crossed_from_above
        else:
            raise FinError("Unknown barrier option type." +
                           str(option_type))

        v = payoff.mean() * np.exp(- r * t)

        return v * notional

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("BARRIER LEVEL", self._barrier_level)
        s += label_to_string("NUM OBSERVATIONS", self._num_obs_per_year)
        s += label_to_string("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
