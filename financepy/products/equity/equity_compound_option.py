##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from financepy.models import equity_compound_option_bs

from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.global_vars import G_DAYS_IN_YEARS, G_SMALL

from ...products.equity.equity_option import EquityOption
from ...market.curves.discount_curve_flat import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...models.equity_compound_option_bs import (
    equity_compound_option_bs,
    equity_compound_option_value_tree,
)

########################################################################################
# TODO: Vectorise pricer
# TODO: Monte Carlo pricer
########################################################################################


class EquityCompoundOption(EquityOption):
    """A EquityCompoundOption is an option which allows the holder
    to either buy or sell another underlying option on a first expiry date that
    itself expires on a second expiry date. Both strikes are set at trade
    initiation."""

    def __init__(
        self,
        c_expiry_dt: Date,  # Compound Option expiry date
        c_opt_type: OptionTypes,  # Compound option type
        c_strike_price: float,  # Compound option strike
        u_expiry_dt: Date,  # Underlying option expiry date
        u_opt_type: OptionTypes,  # Underlying option type
        u_strike_price: float,
    ):  # Underlying option strike price
        """Create the EquityCompoundOption by passing in the first and
        second expiry dates as well as the corresponding strike prices and
        option types."""

        check_argument_types(self.__init__, locals())

        if c_expiry_dt > u_expiry_dt:
            raise FinError("Compound expiry date must precede underlying expiry date")

        if c_opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.AMERICAN_CALL,
            OptionTypes.EUROPEAN_PUT,
            OptionTypes.AMERICAN_PUT,
        ):
            raise FinError("Compound option must be European or American call or put.")

        if u_opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.AMERICAN_CALL,
            OptionTypes.EUROPEAN_PUT,
            OptionTypes.AMERICAN_PUT,
        ):
            raise FinError(
                "Underlying option must be European or American call or put."
            )

        self.c_expiry_dt = c_expiry_dt
        self.c_strike_price = float(c_strike_price)
        self.c_opt_type = c_opt_type

        self.u_expiry_dt = u_expiry_dt
        self.u_strike_price = float(u_strike_price)
        self.u_opt_type = u_opt_type

    ####################################################################################

    def _preprocess_inputs(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Validate inputs and compute (tc, tu, kc, ku, ru, qu, vol)."""

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date")

        if value_dt > self.c_expiry_dt:
            raise FinError("Valuation date after compound expiry date.")

        if value_dt > self.u_expiry_dt:
            raise FinError("Valuation date after underlying expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        tc = (self.c_expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tu = (self.u_expiry_dt - value_dt) / G_DAYS_IN_YEARS
        kc = self.c_strike_price
        ku = self.u_strike_price

        df_u = discount_curve.df(self.u_expiry_dt)
        ru = -np.log(df_u) / tu

        dq_u = dividend_curve.df(self.u_expiry_dt)
        qu = -np.log(dq_u) / tu

        vol = np.maximum(model.volatility, G_SMALL)

        return tc, tu, kc, ku, ru, qu, vol

    ####################################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_steps: int = 200,
    ):
        """Value the compound option using an analytical approach if it is
        entirely European style. Otherwise use a Tree approach to handle the
        early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
        Rubinstein (1991). See also Haug page 132."""

        tc, tu, kc, ku, ru, qu, vol = self._preprocess_inputs(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        v = equity_compound_option_bs(
            self.c_opt_type,
            self.u_opt_type,
            tc,
            tu,
            kc,
            ku,
            stock_price,
            ru,
            qu,
            vol,
            num_steps,
        )

        return v

    ####################################################################################

    def value_tree(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_steps: int = 200,
    ):
        """Value the compound option using a Tree approach to handle the
        early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
        Rubinstein (1991). See also Haug page 132."""

        tc, tu, kc, ku, ru, qu, vol = self._preprocess_inputs(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        v = equity_compound_option_value_tree(
            self.c_opt_type,
            self.u_opt_type,
            tc,
            tu,
            kc,
            ku,
            stock_price,
            ru,
            qu,
            vol,
            num_steps,
        )

        return v

    ####################################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("CPD EXPIRY DATE", self.c_expiry_dt)
        s += label_to_string("CPD STRIKE PRICE", self.c_strike_price)
        s += label_to_string("CPD OPTION TYPE", self.c_opt_type)
        s += label_to_string("UND EXPIRY DATE", self.u_expiry_dt)
        s += label_to_string("UND STRIKE PRICE", self.u_strike_price)
        s += label_to_string("UND OPTION TYPE", self.u_opt_type)
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
