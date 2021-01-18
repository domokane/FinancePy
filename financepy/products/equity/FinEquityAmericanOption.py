##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...products.equity.FinEquityOption import FinEquityOption

from ...models.FinModel import FinModel

###############################################################################
# TODO: Implement some analytical approximations
# TODO: Tree with discrete dividends
# TODO: Other dynamics such as SABR
###############################################################################


class FinEquityAmericanOption(FinEquityOption):
    ''' Class for American (and European) style options on simple vanilla
    calls and puts - a tree valuation model is used that can handle both. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikePrice: float,
                 optionType: FinOptionTypes,
                 numOptions: float = 1.0):
        ''' Class for American style options on simple vanilla calls and puts.
        Specify the expiry date, strike price, whether the option is a call or
        put and the number of options. '''

        checkArgumentTypes(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
            optionType != FinOptionTypes.EUROPEAN_PUT and \
            optionType != FinOptionTypes.AMERICAN_CALL and \
                optionType != FinOptionTypes.AMERICAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        self._expiryDate = expiryDate
        self._strikePrice = strikePrice
        self._optionType = optionType
        self._numOptions = numOptions

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: (np.ndarray, float),
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              model: FinModel):
        ''' Valuation of an American option using a CRR tree to take into
        account the value of early exercise. '''

        if type(valueDate) == FinDate:
            texp = (self._expiryDate - valueDate) / gDaysInYear
        else:
            texp = valueDate

        if np.any(stockPrice <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)

        r = discountCurve.ccRate(self._expiryDate)        
        q = dividendCurve.ccRate(self._expiryDate)

        s = stockPrice
        k = self._strikePrice

        v = model.value(s, texp, k, r, q, self._optionType)
                    
        v = v * self._numOptions

        if isinstance(s, float):
            return v
        else:
            return v[0]

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUMBER", self._numOptions, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
