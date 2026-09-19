##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...products.equity.equity_option import EquityOption
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.helpers import option_years

# from ...models.black_scholes_analytic import baw_value
from ...models.model import Model

########################################################################################
# TODO: Implement some analytical approximations
# TODO: Tree with discrete dividends
# TODO: Other dynamics such as SABR
########################################################################################


class EquityAmericanOption(EquityOption):
    """Class for American (and European) style options on simple vanilla
    calls and puts - a tree valuation model is used that can handle both."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        opt_type: OptionTypes,
        num_options: float = 1.0,
    ):
        """Class for American style options on simple vanilla calls and puts.
        Specify the expiry date, strike price, whether the option is a call or
        put and the number of options."""

        check_argument_types(self.__init__, locals())

        if (
            opt_type != OptionTypes.EUROPEAN_CALL
            and opt_type != OptionTypes.EUROPEAN_PUT
            and opt_type != OptionTypes.AMERICAN_CALL
            and opt_type != OptionTypes.AMERICAN_PUT
        ):
            raise FinError("Unknown OPTION_TYPE" + str(opt_type))

        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.opt_type = opt_type
        self.num_options = num_options

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[np.ndarray, float],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Valuation of an American option using a CRR tree to take into
        account the value of early exercise."""

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_exp = np.maximum(t_exp, 1e-10)

        s = stock_price
        k = self.strike_price

        v = model.value(s, t_exp, k, r, q, self.opt_type)
        v = v * self.num_options

        if isinstance(s, float):
            return v
        else:
            return v[0]

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("NUMBER", self.num_options, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
