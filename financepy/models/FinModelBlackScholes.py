##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

from scipy.stats import norm
from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinError import FinError

N = norm.cdf

from .FinModel import FinModel
from .FinModelCRRTree import crrTreeValAvg
from .FinModelBlackScholesAnalytical import bawValue
from .FinModelBlackScholesAnalytical import bsValue


from enum import Enum


class FinModelBlackScholesTypes(Enum):
        ANALYTICAL = 1
        CRR_TREE = 2
        BARONE_ADESI = 3

###############################################################################

class FinModelBlackScholes(FinModel):
    
    def __init__(self,
                 volatility: float, 
                 implementationType: FinModelBlackScholesTypes = FinModelBlackScholesTypes.ANALYTICAL,
                 parametersDict: dict = None):

        self._volatility = volatility
        self._implementationType = implementationType
        self._parametersDict = parametersDict

    def value(self, 
              spotPrice: float, 
              timeToExpiry: float, 
              strikePrice: float, 
              riskFreeRate: float, 
              dividendRate: float, 
              optionType: FinOptionTypes):

        if self._implementationType == FinModelBlackScholesTypes.ANALYTICAL:

            if optionType == FinOptionTypes.EUROPEAN_CALL:
                phi = +1
            elif optionType == FinOptionTypes.EUROPEAN_PUT:
                phi = -1
            else:
                print(optionType)
                raise FinError("Unsupported Option Type")

            v =  bsValue(spotPrice, 
                         timeToExpiry,
                         strikePrice,
                         riskFreeRate,
                         dividendRate,
                         self._volatility,
                         phi)

            return v

        elif self._implementationType == FinModelBlackScholesTypes.BARONE_ADESI:

            if optionType == FinOptionTypes.AMERICAN_CALL:
                phi = +1
            elif optionType == FinOptionTypes.AMERICAN_PUT:
                phi = -1
            else:
                print(optionType)
                raise FinError("Unsupported Option Type")

            v =  bawValue(spotPrice, 
                          timeToExpiry,
                          strikePrice,
                          riskFreeRate,
                          dividendRate,
                          self._volatility,
                          phi)

            return v

        elif self._implementationType == FinModelBlackScholesTypes.CRR_TREE:

            numStepsPerYear = self._parametersDict["numStepsPerYear"]

            v = crrTreeValAvg(spotPrice, riskFreeRate, 
                              dividendRate, 
                              self._volatility, 
                              numStepsPerYear,
                              timeToExpiry, optionType.value, 
                              strikePrice)['value']

            return v

        else:
            
            raise FinError("Unsupported implementation type")

###############################################################################

