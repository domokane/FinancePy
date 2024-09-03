##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption

# from ...models.black_scholes_analytic import baw_value
from ...models.model import Model

###############################################################################
# TODO: Implement some analytical approximations
# TODO: Tree with discrete dividends
# TODO: Other dynamics such as SABR
###############################################################################


class EquityAmericanOption(EquityOption):
    """Class for American (and European) style options on simple vanilla
    calls and puts - a tree valuation model is used that can handle both."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        option_type: OptionTypes,
        num_options: float = 1.0,
    ):
        """Class for American style options on simple vanilla calls and puts.
        Specify the expiry date, strike price, whether the option is a call or
        put and the number of options."""

        check_argument_types(self.__init__, locals())

        if (
            option_type != OptionTypes.EUROPEAN_CALL
            and option_type != OptionTypes.EUROPEAN_PUT
            and option_type != OptionTypes.AMERICAN_CALL
            and option_type != OptionTypes.AMERICAN_PUT
        ):
            raise FinError("Unknown Option Type" + str(option_type))

        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.option_type = option_type
        self.num_options = num_options

    ###############################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: (np.ndarray, float),
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Valuation of an American option using a CRR tree to take into
        account the value of early exercise."""

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        if isinstance(value_dt, Date):
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_exp = value_dt

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be greater than zero.")

        if isinstance(model, Model) is False:
            raise FinError("Model is not inherited off type FinModel.")

        if np.any(t_exp < 0.0):
            raise FinError("Time to expiry must be positive.")

        t_exp = np.maximum(t_exp, 1e-10)

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        s = stock_price
        k = self.strike_price

        v = model.value(s, t_exp, k, r, q, self.option_type)
        v = v * self.num_options

        if isinstance(s, float):
            return v
        else:
            return v[0]

    ###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.option_type)
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###############################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
