##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption

from ...models.model import Model

###############################################################################
# TODO: Implement some analytical approximations
# TODO: Tree with discrete dividends
# TODO: Other dynamics such as SABR
###############################################################################


class EquityAmericanOption(EquityOption):
    """ Class for American (and European) style options on simple vanilla
    calls and puts - a tree valuation model is used that can handle both. """

    def __init__(self,
                 expiry_date: Date,
                 strike_price: float,
                 option_type: OptionTypes,
                 num_options: float = 1.0):
        """ Class for American style options on simple vanilla calls and puts.
        Specify the expiry date, strike price, whether the option is a call or
        put and the number of options. """

        check_argument_types(self.__init__, locals())

        if option_type != OptionTypes.EUROPEAN_CALL and \
            option_type != OptionTypes.EUROPEAN_PUT and \
            option_type != OptionTypes.AMERICAN_CALL and \
                option_type != OptionTypes.AMERICAN_PUT:
            raise FinError("Unknown Option Type" + str(option_type))

        self._expiry_date = expiry_date
        self._strike_price = strike_price
        self._option_type = option_type
        self._num_options = num_options

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: (np.ndarray, float),
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model: Model):
        """ Valuation of an American option using a CRR tree to take into
        account the value of early exercise. """

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date")

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, Model) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")

        texp = np.maximum(texp, 1e-10)

        r = discount_curve.cc_rate(self._expiry_date)
        q = dividend_curve.cc_rate(self._expiry_date)

        s = stock_price
        k = self._strike_price

        v = model.value(s, texp, k, r, q, self._option_type)
        v = v * self._num_options

        if isinstance(s, float):
            return v
        else:
            return v[0]

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("NUMBER", self._num_options, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
