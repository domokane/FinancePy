##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union, Optional

import numpy as np

from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve

from ...models.model import Model
from ...models.black import Black, implied_volatility


class EquityIndexOption:
    """ Class for managing plain vanilla European/American
    calls and puts on equity indices."""

    def __init__(self,
                 expiry_date: Union[Date, list],
                 strike_price: Union[float, np.ndarray],
                 option_type: OptionTypes,
                 num_options: Optional[float] = 1.0,
                 ):
        """ Create the Equity Index option object by specifying the expiry
        date, the option strike, the option type and the number of options. """

        check_argument_types(self.__init__, locals())

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT,
                           OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            self._option_type = option_type
            self._option_type_value = option_type.value
        else:
            raise FinError("Unknown Option Type")
        self._expiry_date = expiry_date
        self._strike_price = strike_price
        self._num_options = num_options
        self._texp = None

###############################################################################

    def value(self,
              valuation_date: Union[Date, list],
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model,
              ):
        """ Equity Index Option valuation using Black model. """
        if isinstance(valuation_date, Date) is False:
            raise FinError("Valuation date is not a Date")
        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")
        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option value date")
        if isinstance(self._expiry_date, Date):
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        elif isinstance(self._expiry_date, list):
            texp = []
            for expDate in self._expiry_date:
                t = (expDate - valuation_date) / gDaysInYear
            texp.append(t)
            texp = np.array(texp)
        else:
            raise FinError("Valuation date must be Date or list of Date")
        self._texp = texp
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")
        texp = np.maximum(texp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        k = self._strike_price
        if isinstance(model, Black):
            value = model.value(forward_price, k, texp,
                                df, self._option_type)
        else:
            raise FinError("Unknown Model Type")
        value = value * self._num_options
        return value

###############################################################################

    def delta(self,
              valuation_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model):
        """ Calculate delta of a European/American Index option. """
        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date
        self._texp = texp
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")
        texp = np.maximum(texp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        k = self._strike_price
        if isinstance(model, Black):
            delta = model.delta(forward_price, k, texp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return delta

###############################################################################

    def gamma(self,
              valuation_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model):
        """ Calculate gamma of a European/American Index option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")
        texp = np.maximum(texp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        k = self._strike_price
        if isinstance(model, Black):
            gamma = model.gamma(forward_price, k, texp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return gamma

###############################################################################

    def vega(self,
             valuation_date: Date,
             forward_price: float,
             discount_curve: DiscountCurve,
             model: Model):
        """ Calculate vega of a European/American Index option. """

        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")
        texp = np.maximum(texp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        k = self._strike_price
        if isinstance(model, Black):
            vega = model.vega(forward_price, k, texp,
                              df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return vega

###############################################################################

    def theta(self,
              valuation_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model):
        """ Calculate theta of a European/American Index option. """
        if type(valuation_date) == Date:
            texp = (self._expiry_date - valuation_date) / gDaysInYear
        else:
            texp = valuation_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(texp < 0.0):
            raise FinError("Time to expiry must be positive.")
        texp = np.maximum(texp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        k = self._strike_price
        if isinstance(model, Black):
            theta = model.theta(forward_price, k, texp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return theta

###############################################################################

    def implied_volatility(self,
                           valuation_date: Date,
                           forward_price: Union[float, list, np.ndarray],
                           discount_curve: DiscountCurve,
                           model: Model,
                           price: float,
                           ):
        """ Calculate the Black implied volatility of a European/American
        Index option. """
        texp = (self._expiry_date - valuation_date) / gDaysInYear
        if texp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(valuation_date)
        r = -np.log(df)/texp
        if isinstance(model, Black):
            sigma = implied_volatility(
                forward_price, texp, r,
                self._strike_price, price,
                self._option_type)
        else:
            raise FinError("Unknown Model Type")
        return sigma

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE VALUE", self._option_type)
        s += label_to_string("NUMBER", self._num_options, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
