##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np

from ..utils.global_types import OptionTypes
from ..utils.error import FinError

from ..utils.helpers import check_argument_types

from .model import Model
from .equity_crr_tree import crr_tree_val_avg
from .black_scholes_analytic import baw_value
from .black_scholes_analytic import bs_value

from enum import Enum


class BlackScholesTypes(Enum):
    DEFAULT = 0
    ANALYTICAL = 1
    CRR_TREE = 2
    BARONE_ADESI = 3

###############################################################################


class BlackScholes(Model):

    def __init__(self,
                 volatility: (float, np.ndarray),
                 implementationType: BlackScholesTypes = BlackScholesTypes.DEFAULT,
                 num_steps_per_year: int = 100):

        check_argument_types(self.__init__, locals())

        self._volatility = volatility
        self._implementationType = implementationType
        self._num_steps_per_year = num_steps_per_year

    def value(self,
              spotPrice: float,
              time_to_expiry: float,
              strike_price: float,
              risk_free_rate: float,
              dividendRate: float,
              option_type: OptionTypes):

        if option_type == OptionTypes.EUROPEAN_CALL \
                or option_type == OptionTypes.EUROPEAN_PUT:

            if self._implementationType is BlackScholesTypes.DEFAULT:
                self._implementationType = BlackScholesTypes.ANALYTICAL

            if self._implementationType == BlackScholesTypes.ANALYTICAL:

                v = bs_value(spotPrice, time_to_expiry, strike_price,
                             risk_free_rate, dividendRate, self._volatility,
                             option_type.value)

                return v

            elif self._implementationType == BlackScholesTypes.CRR_TREE:

                v = crr_tree_val_avg(spotPrice, risk_free_rate, dividendRate,
                                     self._volatility, self._num_steps_per_year,
                                     time_to_expiry, option_type.value,
                                     strike_price)['value']

                return v

            else:

                raise FinError("Implementation not available for this product")

        elif option_type == OptionTypes.AMERICAN_CALL \
                or option_type == OptionTypes.AMERICAN_PUT:

            if self._implementationType is BlackScholesTypes.DEFAULT:
                self._implementationType = BlackScholesTypes.CRR_TREE

            if self._implementationType == BlackScholesTypes.BARONE_ADESI:

                if option_type == OptionTypes.AMERICAN_CALL:
                    phi = +1
                elif option_type == OptionTypes.AMERICAN_PUT:
                    phi = -1

                v = baw_value(spotPrice, time_to_expiry, strike_price,
                              risk_free_rate, dividendRate, self._volatility,
                              phi)

                return v

            elif self._implementationType == BlackScholesTypes.CRR_TREE:

                v = crr_tree_val_avg(spotPrice, risk_free_rate, dividendRate,
                                     self._volatility, self._num_steps_per_year,
                                     time_to_expiry, option_type.value,
                                     strike_price)['value']

                return v

            else:

                raise FinError("Implementation not available for this product")

        else:

            raise FinError("Should not be here")

###############################################################################
