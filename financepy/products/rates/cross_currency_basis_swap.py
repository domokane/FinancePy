"""Fixed-notional cross-currency floating/floating basis swap.

The FX quote is DOM units per one FOR unit (e.g. EURUSD = USD per EUR).
Each floating leg uses a discount and an index curve in its own currency.
PV is expressed in domestic currency by converting the foreign leg at
valuation-date spot FX. This is a two-curve domestic/foreign valuation,
not a collateral-currency or stochastic-FX basis model.
"""
from __future__ import annotations

import numpy as np

from ...market.curves.discount_curve import DiscountCurve
from ...utils.calendar import CalendarTypes, BusDayAdjustTypes, DateGenRuleTypes
from ...utils.check_values import check_curve_dt, check_currency_pair
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.error import FinError
from ...utils.frequency import FrequencyTypes
from ...utils.global_types import SwapTypes

from .swap_float_leg import SwapFloatLeg


class CrossCurrencyBasisSwap:
    """Two floating-rate legs in different currencies with optional par exchange.

    For a domestic PAY leg, the foreign leg RECEIVEs, and vice versa.
    With exchange_principal=True, initial par is exchanged opposite to the
    coupon-leg direction and terminal par in the coupon-leg direction.
    No intermediate FX reset or mark-to-market notional resets are modeled.
    """

    def __init__(
        self,
        effective_dt: Date,
        term_dt_or_tenor: Date | str,
        currency_pair: str,
        domestic_leg_type: SwapTypes,
        domestic_notional: float,
        foreign_notional: float,
        domestic_spread: float = 0.0,
        foreign_spread: float = 0.0,
        domestic_frequency: FrequencyTypes = FrequencyTypes.QUARTERLY,
        foreign_frequency: FrequencyTypes = FrequencyTypes.QUARTERLY,
        domestic_day_count: DayCountTypes = DayCountTypes.ACT_360,
        foreign_day_count: DayCountTypes = DayCountTypes.ACT_360,
        exchange_principal: bool = True,
        calendar: CalendarTypes = CalendarTypes.WEEKEND,
        business_day_adjust: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        date_gen_rule: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        check_currency_pair(currency_pair)
        if domestic_leg_type not in (SwapTypes.PAY, SwapTypes.RECEIVE):
            raise FinError("Domestic leg must be PAY or RECEIVE.")
        if (not np.isfinite([domestic_notional, foreign_notional,
                             domestic_spread, foreign_spread]).all()
                or domestic_notional <= 0 or foreign_notional <= 0):
            raise FinError("Notionals must be positive; spreads finite.")

        self.currency_pair = currency_pair.upper()
        self.domestic_ccy = self.currency_pair[3:]
        self.foreign_ccy = self.currency_pair[:3]
        self.domestic_notional = domestic_notional
        self.foreign_notional = foreign_notional
        self.effective_dt = effective_dt
        self.exchange_principal = exchange_principal
        self.domestic_leg_type = domestic_leg_type
        self.foreign_leg_type = (
            SwapTypes.RECEIVE if domestic_leg_type == SwapTypes.PAY
            else SwapTypes.PAY
        )
        principal = 1.0 if exchange_principal else 0.0

        self.domestic_leg = SwapFloatLeg(
            effective_dt, term_dt_or_tenor, domestic_leg_type,
            domestic_spread, domestic_frequency, domestic_day_count,
            domestic_notional, principal, 0, calendar,
            business_day_adjust, date_gen_rule,
        )
        self.foreign_leg = SwapFloatLeg(
            effective_dt, term_dt_or_tenor, self.foreign_leg_type,
            foreign_spread, foreign_frequency, foreign_day_count,
            foreign_notional, principal, 0, calendar,
            business_day_adjust, date_gen_rule,
        )
        self.maturity_dt = self.domestic_leg.maturity_dt

    def value(
        self,
        value_dt: Date,
        spot_fx: float,
        domestic_discount: DiscountCurve,
        foreign_discount: DiscountCurve,
        domestic_index: DiscountCurve | None = None,
        foreign_index: DiscountCurve | None = None,
        domestic_first_fixing: float | None = None,
        foreign_first_fixing: float | None = None,
    ) -> float:
        """Present value in domestic currency at the supplied spot FX."""
        if not np.isfinite(spot_fx) or spot_fx <= 0.0:
            raise FinError("spot_fx must be positive and finite.")
        check_curve_dt(value_dt, domestic_discount, foreign_discount)
        domestic_index = (domestic_discount if domestic_index is None
                          else domestic_index)
        foreign_index = (foreign_discount if foreign_index is None
                         else foreign_index)
        check_curve_dt(value_dt, domestic_index, foreign_index)

        pv_dom = self.domestic_leg.value(
            value_dt, domestic_discount, domestic_index, domestic_first_fixing
        )
        pv_for = self.foreign_leg.value(
            value_dt, foreign_discount, foreign_index, foreign_first_fixing
        )
        value = pv_dom + spot_fx * pv_for

        if self.exchange_principal and value_dt <= self.effective_dt:
            # Initial principal cashflows are opposite to the signed
            # coupon/terminal principal cashflows, if not settled already.
            dom_sign = -1.0 if self.domestic_leg_type == SwapTypes.PAY else 1.0
            for_sign = -dom_sign
            df_dom = domestic_discount.df(self.effective_dt) / domestic_discount.df(value_dt)
            df_for = foreign_discount.df(self.effective_dt) / foreign_discount.df(value_dt)
            value -= (
                dom_sign * self.domestic_notional * df_dom
                + spot_fx * for_sign * self.foreign_notional * df_for
            )
        return value

    def par_foreign_spread(
        self,
        value_dt: Date,
        spot_fx: float,
        domestic_discount: DiscountCurve,
        foreign_discount: DiscountCurve,
        domestic_index: DiscountCurve | None = None,
        foreign_index: DiscountCurve | None = None,
        domestic_first_fixing: float | None = None,
        foreign_first_fixing: float | None = None,
    ) -> float:
        """Foreign floating spread that zeroes domestic-currency PV."""
        pv = self.value(
            value_dt, spot_fx, domestic_discount, foreign_discount,
            domestic_index, foreign_index, domestic_first_fixing, foreign_first_fixing
        )
        sign = 1.0 if self.foreign_leg_type == SwapTypes.RECEIVE else -1.0
        df_value = foreign_discount.df(value_dt)
        annuity = sum(
            alpha * foreign_discount.df(payment_dt) / df_value
            for alpha, payment_dt in zip(
                self.foreign_leg.year_fracs, self.foreign_leg.payment_dts
            ) if payment_dt > value_dt
        )
        derivative = sign * spot_fx * self.foreign_notional * annuity
        if derivative == 0.0:
            raise FinError("No outstanding foreign coupon accruals.")
        return self.foreign_leg.spread - pv / derivative

