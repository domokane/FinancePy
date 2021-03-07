##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.FinError import FinError
from ...utils.global_types import FinOptionTypes
from ...utils.helpers import check_argument_types, labelToString
from ...market.discount.curve import DiscountCurve
from ...products.equity.FinEquityOption import FinEquityOption

from ...models.FinModel import FinModel

###############################################################################
# TODO: Implement some analytical approximations
# TODO: Tree with discrete dividends
# TODO: Other dynamics such as SABR
###############################################################################


class FinEquityAmericanOption(FinEquityOption):
    """ Class for American (and European) style options on simple vanilla
    calls and puts - a tree valuation model is used that can handle both. """

    def __init__(self,
                 expiry_date: Date,
                 strikePrice: float,
                 optionType: FinOptionTypes,
                 numOptions: float = 1.0):
        """ Class for American style options on simple vanilla calls and puts.
        Specify the expiry date, strike price, whether the option is a call or
        put and the number of options. """

        check_argument_types(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
            optionType != FinOptionTypes.EUROPEAN_PUT and \
            optionType != FinOptionTypes.AMERICAN_CALL and \
                optionType != FinOptionTypes.AMERICAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        self._expiry_date = expiry_date
        self._strikePrice = strikePrice
        self._optionType = optionType
        self._numOptions = numOptions

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: (np.ndarray, float),
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model: FinModel):
        """ Valuation of an American option using a CRR tree to take into
        account the value of early exercise. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, FinModel) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)

        r = discount_curve.ccRate(self._expiry_date)        
        q = dividendCurve.ccRate(self._expiry_date)

        s = stock_price
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
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUMBER", self._numOptions, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
