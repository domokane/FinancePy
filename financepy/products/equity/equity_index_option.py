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
        self._t_exp = None

###############################################################################

    def value(self,
              value_date: Union[Date, list],
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model,
              ):
        """ Equity Index Option valuation using Black model. """
        if isinstance(value_date, Date) is False:
            raise FinError("Valuation date is not a Date")
        if value_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")
        if discount_curve._value_date != value_date:
            raise FinError(
                "Discount Curve valuation date not same as option value date")
        if isinstance(self._expiry_date, Date):
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        elif isinstance(self._expiry_date, list):
            t_exp = []
            for expDate in self._expiry_date:
                t = (expDate - value_date) / gDaysInYear
            t_exp.append(t)
            t_exp = np.array(t_exp)
        else:
            raise FinError("Valuation date must be Date or list of Date")
        self._t_exp = t_exp
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        k = self._strike_price
        if isinstance(model, Black):
            value = model.value(forward_price, k, t_exp,
                                df, self._option_type)
        else:
            raise FinError("Unknown Model Type")
        value = value * self._num_options
        return value

###############################################################################

    def delta(self,
              value_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model):
        """ Calculate delta of a European/American Index option. """

        if isinstance(value_date, Date):
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        else:
            t_exp = value_date
        self._t_exp = t_exp
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        k = self._strike_price
        if isinstance(model, Black):
            delta = model.delta(forward_price, k, t_exp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return delta

###############################################################################

    def gamma(self,
              value_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model):
        """ Calculate gamma of a European/American Index option. """

        if isinstance(value_date, Date):
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        else:
            t_exp = value_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        k = self._strike_price
        if isinstance(model, Black):
            gamma = model.gamma(forward_price, k, t_exp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return gamma

###############################################################################

    def vega(self,
             value_date: Date,
             forward_price: float,
             discount_curve: DiscountCurve,
             model: Model):
        """ Calculate vega of a European/American Index option. """

        if isinstance(value_date, Date):
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        else:
            t_exp = value_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        k = self._strike_price
        if isinstance(model, Black):
            vega = model.vega(forward_price, k, t_exp,
                              df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return vega

###############################################################################

    def theta(self,
              value_date: Date,
              forward_price: float,
              discount_curve: DiscountCurve,
              model: Model):
        """ Calculate theta of a European/American Index option. """

        if isinstance(value_date, Date):
            t_exp = (self._expiry_date - value_date) / gDaysInYear
        else:
            t_exp = value_date
        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        k = self._strike_price
        if isinstance(model, Black):
            theta = model.theta(forward_price, k, t_exp,
                                df, self._option_type_value)
        else:
            raise FinError("Unknown Model Type")
        return theta

###############################################################################

    def implied_volatility(self,
                           value_date: Date,
                           forward_price: Union[float, list, np.ndarray],
                           discount_curve: DiscountCurve,
                           model: Model,
                           price: float,
                           ):
        """ Calculate the Black implied volatility of a European/American
        Index option. """
        t_exp = (self._expiry_date - value_date) / gDaysInYear
        if t_exp < 1.0 / 365.0:
            print("Expiry time is too close to zero.")
            return -999
        df = discount_curve.df(self._expiry_date) / \
            discount_curve.df(value_date)
        r = -np.log(df)/t_exp
        if isinstance(model, Black):
            sigma = implied_volatility(
                forward_price, t_exp, r,
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
