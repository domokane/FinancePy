##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np

from ..utils.global_types import FinOptionTypes
from ..utils.error import FinError

from ..utils.helpers import check_argument_types

from .FinModel import FinModel
from .equity_crr_tree import crrTreeValAvg
from .black_scholes_analytic import bawValue
from .black_scholes_analytic import bs_value

from enum import Enum

class FinModelBlackScholesTypes(Enum):
        DEFAULT = 0
        ANALYTICAL = 1
        CRR_TREE = 2
        BARONE_ADESI = 3

###############################################################################

class BlackScholes(FinModel):
    
    def __init__(self,
                 volatility: (float, np.ndarray), 
                 implementationType: FinModelBlackScholesTypes = FinModelBlackScholesTypes.DEFAULT,
                 num_steps_per_year: int = 100):

        check_argument_types(self.__init__, locals())

        self._volatility = volatility
        self._implementationType = implementationType
        self._num_steps_per_year = num_steps_per_year

    def value(self, 
              spotPrice: float, 
              time_to_expiry: float,
              strike_price: float,
              riskFreeRate: float, 
              dividendRate: float, 
              option_type: FinOptionTypes):

        if option_type == FinOptionTypes.EUROPEAN_CALL \
            or option_type == FinOptionTypes.EUROPEAN_PUT:

            if self._implementationType is FinModelBlackScholesTypes.DEFAULT:
                self._implementationType = FinModelBlackScholesTypes.ANALYTICAL

            if self._implementationType == FinModelBlackScholesTypes.ANALYTICAL:

                v =  bs_value(spotPrice, time_to_expiry, strike_price,
                             riskFreeRate, dividendRate, self._volatility,
                             option_type.value)

                return v

            elif self._implementationType == FinModelBlackScholesTypes.CRR_TREE:
                
                v = crrTreeValAvg(spotPrice, riskFreeRate, dividendRate, 
                                  self._volatility, self._num_steps_per_year,
                                  time_to_expiry, option_type.value,
                                  strike_price)['value']

                return v

            else:
                
                raise FinError("Implementation not available for this product")

        elif option_type == FinOptionTypes.AMERICAN_CALL \
            or option_type == FinOptionTypes.AMERICAN_PUT:

            if self._implementationType is FinModelBlackScholesTypes.DEFAULT:
                self._implementationType = FinModelBlackScholesTypes.CRR_TREE

            if self._implementationType == FinModelBlackScholesTypes.BARONE_ADESI:

                if option_type == FinOptionTypes.AMERICAN_CALL:
                    phi = +1
                elif option_type == FinOptionTypes.AMERICAN_PUT:
                    phi = -1

                v =  bawValue(spotPrice, time_to_expiry, strike_price,
                              riskFreeRate, dividendRate, self._volatility,
                              phi)

                return v

            elif self._implementationType == FinModelBlackScholesTypes.CRR_TREE:

                v = crrTreeValAvg(spotPrice, riskFreeRate, dividendRate, 
                                  self._volatility, self._num_steps_per_year,
                                  time_to_expiry, option_type.value,
                                  strike_price)['value']

                return v

            else:
                
                raise FinError("Implementation not available for this product")

        else:
            
            raise FinError("Should not be here")

###############################################################################

