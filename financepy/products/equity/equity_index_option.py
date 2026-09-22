##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union, Optional

import numpy as np

from ...utils.date import Date
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve

from ...models.model import Model
from ...models.black import Black, implied_volatility

from ...utils.check_values import check_curve_dt
from ...utils.helpers import option_years


class EquityIndexOption:
    """Class for managing plain vanilla European/American
    calls and puts on equity indices."""

    def __init__(
        self,
        expiry_dt: Union[Date, list],
        strike_price: Union[float, np.ndarray],
        opt_type: OptionTypes,
        num_options: Optional[float] = 1.0,
    ):
        """Create the Equity Index option object by specifying the expiry
        date, the option strike, the option type and the number of options."""

        check_argument_types(self.__init__, locals())

        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):
            self.opt_type = opt_type
            self.opt_type_value = opt_type.value
        else:
            raise FinError("Unknown Option Type")
        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.num_options = num_options
        self.t_exp = None

    ###########################################################################

    def value(
        self,
        value_dt: Union[Date, list],
        forward_price: float,
        discount_curve: DiscountCurve,
        model: Model,
    ):
        """Equity Index Option valuation using Black model."""

        check_curve_dt(value_dt, discount_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")

        df = discount_curve.df(self.expiry_dt) / discount_curve.df(value_dt)

        k = self.strike_price

        if isinstance(model, Black):
            value = model.value(forward_price, k, t_exp, df, self.opt_type)
        else:
            raise FinError("Unknown Model Type")

        value = value * self.num_options
        return value

    ###########################################################################

    def delta(
        self,
        value_dt: Date,
        forward_price: float,
        discount_curve: DiscountCurve,
        model,
    ):
        """Calculate delta of a European/American Index option."""

        check_curve_dt(value_dt, discount_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")

        df = discount_curve.df(self.expiry_dt) / discount_curve.df(value_dt)
        k = self.strike_price

        if isinstance(model, Black):
            delta = model.delta(forward_price, k, t_exp, df, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")

        return delta

    ###########################################################################

    def gamma(
        self,
        value_dt: Date,
        forward_price: float,
        discount_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate gamma of a European/American Index option."""
        check_curve_dt(value_dt, discount_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt) / discount_curve.df(value_dt)
        k = self.strike_price
        if isinstance(model, Black):
            gamma = model.gamma(forward_price, k, t_exp, df, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")
        return gamma

    ###########################################################################

    def vega(
        self,
        value_dt: Date,
        forward_price: float,
        discount_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate vega of a European/American Index option."""
        check_curve_dt(value_dt, discount_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        if np.any(forward_price <= 0.0):
            raise FinError("Forward price must be greater than zero.")
        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")
        t_exp = np.maximum(t_exp, 1e-10)
        df = discount_curve.df(self.expiry_dt) / discount_curve.df(value_dt)
        k = self.strike_price
        if isinstance(model, Black):
            vega = model.vega(forward_price, k, t_exp, df, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")
        return vega

    ###########################################################################

    def theta(
        self,
        value_dt: Date,
        forward_price: float,
        discount_curve: DiscountCurve,
        model: Model,
    ):
        """Calculate theta of a European/American Index option."""
        check_curve_dt(value_dt, discount_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        df = discount_curve.df(self.expiry_dt) / discount_curve.df(value_dt)
        k = self.strike_price
        if isinstance(model, Black):
            theta = model.theta(forward_price, k, t_exp, df, self.opt_type_value)
        else:
            raise FinError("Unknown Model Type")
        return theta

    ###########################################################################

    def implied_volatility(
        self,
        value_dt: Date,
        forward_price: Union[float, list, np.ndarray],
        discount_curve: DiscountCurve,
        model: Model,
        price: float,
    ):
        check_curve_dt(value_dt, discount_curve)

        """Calculate Black implied volatility of a European/American Index option."""
        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        if t_exp < 1.0 / 366.0:
            print("Expiry time is too close to zero.")
            return -999

        r = discount_curve.fwd_zero_rate_cc(value_dt, self.expiry_dt)

        if isinstance(model, Black):
            sigma = implied_volatility(
                forward_price,
                t_exp,
                r,
                self.strike_price,
                price,
                self.opt_type,
            )
        else:
            raise FinError("Unknown Model Type")
        return sigma

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE VALUE", self.opt_type)
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
