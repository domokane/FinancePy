##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from copy import deepcopy
from typing import Union

from financepy.models.model import Model

from ...utils.global_vars import G_DAYS_IN_YEAR
from ...market.curves.discount_curve import DiscountCurve
from ...utils.date import Date

########################################################################################


########################################################################################


class EquityOption:
    """Parent class for equity options requiring perturbatory risk."""

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: Union[DiscountCurve, float],
        dividend_curve: Union[DiscountCurve, float],
        model: Model,
        **kwargs,
    ):
        raise NotImplementedError("value() must be implemented by subclass")

    ###########################################################################

    def delta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option delta by perturbation of stock price."""

        bump = stock_price * 1e-3

        v_dn = self.value(
            value_dt,
            stock_price - bump,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        v_up = self.value(
            value_dt,
            stock_price + bump,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        return (v_up - v_dn) / (2.0 * bump)

    ###########################################################################

    def gamma(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option gamma by perturbation of stock price."""

        bump = 1e-3 * stock_price

        v = self.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        v_bumped_dn = self.value(
            value_dt,
            stock_price - bump,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        v_bumped_up = self.value(
            value_dt,
            stock_price + bump,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        return (v_bumped_up - 2.0 * v + v_bumped_dn) / bump**2

    ###########################################################################

    def vega(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option vega by perturbing volatility by 1%."""

        bump = 0.01

        v = self.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        bumped_model = deepcopy(model)
        bumped_model.volatility += bump

        v_bumped = self.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            bumped_model,
            **kwargs,
        )

        return v_bumped - v

    ###########################################################################

    def vanna(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option vanna by perturbing delta with respect to vol."""
        bump = 0.01

        delta = self.delta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        bumped_model = deepcopy(model)
        bumped_model.volatility += bump

        delta_bumped = self.delta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            bumped_model,
            **kwargs,
        )

        return (delta_bumped - delta) / bump

    ###########################################################################

    def theta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option theta by moving valuation date one calendar day."""

        v = self.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        next_dt = value_dt.add_days(1)

        old_discount_anchor_dt = discount_curve.anchor_dt
        old_dividend_anchor_dt = dividend_curve.anchor_dt

        discount_curve.anchor_dt = next_dt
        dividend_curve.anchor_dt = next_dt

        time_bump = (next_dt - value_dt) / G_DAYS_IN_YEAR

        try:
            v_bumped = self.value(
                next_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                **kwargs,
            )
        finally:
            # Always restore curves, even if valuation raises.
            discount_curve.anchor_dt = old_discount_anchor_dt
            dividend_curve.anchor_dt = old_dividend_anchor_dt

        return (v_bumped - v) / time_bump

    ###########################################################################

    def rho(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        **kwargs,
    ):
        """Calculate option rho by perturbing the interest-rate curve."""

        bump = 0.0001

        v = self.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            **kwargs,
        )

        v_bumped = self.value(
            value_dt,
            stock_price,
            discount_curve.bump_parallel(bump),
            dividend_curve,
            model,
            **kwargs,
        )

        return (v_bumped - v) / bump

    ###########################################################################
