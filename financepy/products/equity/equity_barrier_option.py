###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

from ...market.curves.discount_curve import DiscountCurve
from ...models.equity_barrier_models import value_one
from ...products.equity.equity_option import EquityOption
from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import EquityBarrierTypes
from ...utils.helpers import label_to_string, check_argument_types


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
                 num_observations_per_year: (int, float) = 252,
                 notional: float = 1.0):
        """ Create the EquityBarrierOption by specifying the expiry date,
        strike price, option type, barrier level, the number of observations
        per year and the notional. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strike_price = float(strike_price)
        self._barrier_level = float(barrier_level)
        self._num_observations_per_year = int(num_observations_per_year)

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

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date")

        if isinstance(stock_price, int):
            stock_price = float(stock_price)

        if isinstance(stock_price, float):
            stock_prices = [stock_price]
        else:
            stock_prices = stock_price

        values = []
        for s in stock_prices:
            v = value_one(self._expiry_date, self._strike_price, self._option_type.value, self._barrier_level,
                          self._num_observations_per_year, self._notional, valuation_date, s,
                          discount_curve.cc_rate(self._expiry_date),
                          dividend_curve.cc_rate(self._expiry_date), model)
            values.append(v)

        if isinstance(stock_price, float):
            return values[0]
        else:
            return np.array(values)

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("BARRIER LEVEL", self._barrier_level)
        s += label_to_string("NUM OBSERVATIONS",
                             self._num_observations_per_year)
        s += label_to_string("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
