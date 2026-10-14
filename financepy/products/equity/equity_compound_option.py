##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.global_vars import G_SMALL

from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...models.bs_equity_compound_option import equity_compound_option_value_tree
from ...models.bs_equity_compound_option import equity_compound_option_bs
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_volatility
from ...utils.helpers import option_years
from ...models.model import Model

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
        c_opt_type: OptionTypes,  # Compound OPTION_TYPE
        c_strike_price: float,  # Compound option strike
        u_expiry_dt: Date,  # Underlying option expiry date
        u_opt_type: OptionTypes,  # Underlying OPTION_TYPE
        u_strike_price: float,
    ) -> None:  # Underlying option strike price
        """Create the EquityCompoundOption by passing in the first and
        second expiry dates as well as the corresponding strike prices and
        OPTION_TYPEs."""

        check_argument_types(self.__init__, locals())

        if c_expiry_dt >= u_expiry_dt:
            raise FinError("Compound expiry date must precede underlying expiry date")

        if c_opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            raise FinError("Compound option must be European call or put.")

        if u_opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            raise FinError("Underlying option must be European call or put.")

        if not np.isfinite(c_strike_price) or c_strike_price < 0.0:
            raise FinError("Invalid compound strike price")

        if not np.isfinite(u_strike_price) or u_strike_price <= 0.0:
            raise FinError("Invalid underlying strike price")

        self.c_expiry_dt = c_expiry_dt
        self.c_strike_price = float(c_strike_price)
        self.c_opt_type = c_opt_type

        self.u_expiry_dt = u_expiry_dt
        self.u_strike_price = float(u_strike_price)
        self.u_opt_type = u_opt_type

    ####################################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_steps: int = 200,
    ):
        """Value the compound option using an analytical approach if it is
        entirely European style. Otherwise use a Tree approach to handle the
        early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
        Rubinstein (1991). See also Haug page 132."""

        if value_dt >= self.c_expiry_dt:
            raise FinError(
                "Valuation date must precede compound option expiry"
            )

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_volatility(model.volatility)

        if not np.isfinite(stock_price) or stock_price <= 0.0:
            raise FinError("Stock price must be positive")

        rc = discount_curve.zero_rate_cc(self.c_expiry_dt)
        qc = dividend_curve.zero_rate_cc(self.c_expiry_dt)

        ru = discount_curve.zero_rate_cc(self.u_expiry_dt)
        qu = dividend_curve.zero_rate_cc(self.u_expiry_dt)

        tc = option_years(value_dt, self.c_expiry_dt)
        tu = option_years(value_dt, self.u_expiry_dt)

        kc = self.c_strike_price
        ku = self.u_strike_price

        vol = np.maximum(model.volatility, G_SMALL)

        v = equity_compound_option_bs(
            self.c_opt_type.value,
            self.u_opt_type.value,
            tc,
            tu,
            kc,
            ku,
            stock_price,
            rc,
            ru,
            qc,
            qu,
            vol
        )

        return v

    ####################################################################################

    def value_tree(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_steps: int = 200,
    ):
        """Value the compound option using a Tree approach to handle the
        early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
        Rubinstein (1991). See also Haug page 132."""

        if value_dt >= self.c_expiry_dt:
            raise FinError(
                "Valuation date must precede compound option expiry"
            )

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_volatility(model.volatility)

        if not np.isfinite(stock_price) or stock_price <= 0.0:
            raise FinError("Stock price must be positive")

        if isinstance(num_steps, (bool, np.bool_)) or not isinstance(
            num_steps, (int, np.integer)
        ) or num_steps < 2:
            raise FinError("num_steps must be an integer >= 2")

        ru = discount_curve.zero_rate_cc(self.u_expiry_dt)
        qu = dividend_curve.zero_rate_cc(self.u_expiry_dt)

        tc = option_years(value_dt, self.c_expiry_dt)
        tu = option_years(value_dt, self.u_expiry_dt)

        kc = self.c_strike_price
        ku = self.u_strike_price

        vol = np.maximum(model.volatility, G_SMALL)

        v = equity_compound_option_value_tree(
            self.c_opt_type.value,
            self.u_opt_type.value,
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
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("CPD EXPIRY DATE", self.c_expiry_dt)
        s += label_to_string("CPD STRIKE PRICE", self.c_strike_price)
        s += label_to_string("CPD OPTION_TYPE", self.c_opt_type)
        s += label_to_string("UND EXPIRY DATE", self.u_expiry_dt)
        s += label_to_string("UND STRIKE PRICE", self.u_strike_price)
        s += label_to_string("UND OPTION_TYPE", self.u_opt_type)
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
