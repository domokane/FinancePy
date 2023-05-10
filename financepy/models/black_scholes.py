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
from .equity_lsmc import equity_lsmc, FIT_TYPES

from .black_scholes_analytic import (
    bs_value,
    baw_value,
    bjerksund_stensland_value
)
from .finite_difference import black_scholes_finite_difference
from .finite_difference_PSOR import black_scholes_fd_PSOR


from enum import Enum


class BlackScholesTypes(Enum):
    DEFAULT = 0
    ANALYTICAL = 1
    CRR_TREE = 2
    BARONE_ADESI = 3
    LSMC = 4
    Bjerksund_Stensland = 5
    FINITE_DIFFERENCE = 6
    PSOR = 7

###############################################################################


class BlackScholes(Model):

    ###########################################################################

    def __init__(self,
                 volatility: (float, np.ndarray),
                 implementationType: BlackScholesTypes = BlackScholesTypes.DEFAULT,
                 num_steps_per_year=52,
                 num_paths=10000,
                 seed=42,
                 use_sobol=False,
                 params=None):

        check_argument_types(self.__init__, locals())

        self._volatility = volatility
        self._implementationType = implementationType
        self._num_steps_per_year = num_steps_per_year
        self._num_paths = num_paths
        self._seed = seed
        self._use_sobol = use_sobol
        self._params = params if params else {}

    ###########################################################################

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

                v = bs_value(spotPrice,
                             time_to_expiry,
                             strike_price,
                             risk_free_rate,
                             dividendRate,
                             self._volatility,
                             option_type.value)

                return v

            elif self._implementationType == BlackScholesTypes.CRR_TREE:

                v = crr_tree_val_avg(spotPrice,
                                     risk_free_rate,
                                     dividendRate,
                                     self._volatility,
                                     self._num_steps_per_year,
                                     time_to_expiry,
                                     option_type.value,
                                     strike_price)['value']

                return v

            elif self._implementationType == BlackScholesTypes.FINITE_DIFFERENCE:
                v = black_scholes_finite_difference(spot_price=spotPrice,
                                                    time_to_expiry=time_to_expiry,
                                                    strike_price=strike_price,
                                                    risk_free_rate=risk_free_rate,
                                                    dividend_yield=dividendRate,
                                                    volatility=self._volatility,
                                                    option_type=option_type.value,
                                                    **self._params
                                                    )
                return v
            elif self._implementationType == BlackScholesTypes.PSOR:
                v = black_scholes_fd_PSOR(spot_price=spotPrice,
                                          time_to_expiry=time_to_expiry,
                                          strike_price=strike_price,
                                          risk_free_rate=risk_free_rate,
                                          dividend_yield=dividendRate,
                                          volatility=self._volatility,
                                          option_type=option_type.value,
                                          **self._params
                                          )
                return v
            elif self._implementationType == BlackScholesTypes.LSMC:

                v = equity_lsmc(spot_price=spotPrice,
                                risk_free_rate=risk_free_rate,
                                dividend_yield=dividendRate,
                                sigma=self._volatility,
                                num_paths=self._num_paths,
                                num_steps_per_year=self._num_steps_per_year,
                                time_to_expiry=time_to_expiry,
                                option_type_value=option_type.value,
                                strike_price=strike_price,
                                seed=self._seed,
                                use_sobol=self._use_sobol,
                                **self._params)

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

                v = baw_value(spotPrice,
                              time_to_expiry,
                              strike_price,
                              risk_free_rate,
                              dividendRate,
                              self._volatility,
                              phi)

                return v

            elif self._implementationType == BlackScholesTypes.CRR_TREE:

                v = crr_tree_val_avg(spotPrice,
                                     risk_free_rate,
                                     dividendRate,
                                     self._volatility,
                                     self._num_steps_per_year,
                                     time_to_expiry,
                                     option_type.value,
                                     strike_price)['value']

                return v

            elif self._implementationType == BlackScholesTypes.LSMC:

                poly_degree = 3
                fit_type = FIT_TYPES.HERMITE_E

                v = equity_lsmc(spot_price=spotPrice,
                                risk_free_rate=risk_free_rate,
                                dividend_yield=dividendRate,
                                sigma=self._volatility,
                                num_paths=self._num_paths,
                                num_steps_per_year=self._num_steps_per_year,
                                time_to_expiry=time_to_expiry,
                                option_type_value=option_type.value,
                                strike_price=strike_price,
                                poly_degree = poly_degree,
                                fit_type_value = fit_type.value,
                                seed=self._seed,
                                use_sobol=self._use_sobol,
                                **self._params)

                return v
            elif self._implementationType == BlackScholesTypes.Bjerksund_Stensland:
                v = bjerksund_stensland_value(spotPrice,
                                              time_to_expiry,
                                              strike_price,
                                              risk_free_rate,
                                              dividendRate,
                                              self._volatility,
                                              option_type.value)
                return v

            elif self._implementationType == BlackScholesTypes.FINITE_DIFFERENCE:
                v = black_scholes_finite_difference(spot_price=spotPrice,
                                                    time_to_expiry=time_to_expiry,
                                                    strike_price=strike_price,
                                                    risk_free_rate=risk_free_rate,
                                                    dividend_yield=dividendRate,
                                                    volatility=self._volatility,
                                                    option_type=option_type.value,
                                                    **self._params
                                                    )
                return v

            elif self._implementationType == BlackScholesTypes.PSOR:
                v = black_scholes_fd_PSOR(spot_price=spotPrice,
                                          time_to_expiry=time_to_expiry,
                                          strike_price=strike_price,
                                          risk_free_rate=risk_free_rate,
                                          dividend_yield=dividendRate,
                                          volatility=self._volatility,
                                          option_type=option_type.value,
                                          **self._params
                                          )
                return v

            else:

                raise FinError("Implementation not available for this product")

        else:

            raise FinError("Should not be here")

###############################################################################
