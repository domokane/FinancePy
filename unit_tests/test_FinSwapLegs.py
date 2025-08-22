# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.rates.swap_float_leg import SwapFloatLeg
from financepy.products.rates.swap_fixed_leg import SwapFixedLeg
from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_types import SwapTypes
from financepy.utils.math import ONE_MILLION

########################################################################################


def test__fin_fixed_ibor_swap_leg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    coupon = -0.44970 / 100.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    leg_pay_rec_type = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swap_fixed_leg = SwapFixedLeg(
        effective_dt,
        maturity_dt,
        leg_pay_rec_type,
        coupon,
        freq_type,
        dc_type,
        notional,
        principal,
        payment_lag,
        cal_type,
        bd_type,
        dg_type,
    )

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    v = swap_fixed_leg.value(effective_dt, libor_curve)
    assert round(v, 4) == 194018.0116

########################################################################################


def test__fin_fixed_ois_swap_leg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    coupon = -0.515039 / 100.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    leg_pay_rec_type = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swap_fixed_leg = SwapFixedLeg(
        effective_dt,
        maturity_dt,
        leg_pay_rec_type,
        coupon,
        freq_type,
        dc_type,
        notional,
        principal,
        payment_lag,
        cal_type,
        bd_type,
        dg_type,
    )

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    v = swap_fixed_leg.value(effective_dt, libor_curve)
    assert round(v, 4) == 225367.1730

########################################################################################


def test__fin_float_ibor_leg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    leg_pay_rec_type = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swap_float_leg = SwapFloatLeg(
        effective_dt,
        maturity_dt,
        leg_pay_rec_type,
        spread,
        freq_type,
        dc_type,
        notional,
        principal,
        payment_lag,
        cal_type,
        bd_type,
        dg_type,
    )

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    first_fixing = 0.03

    v = swap_float_leg.value(effective_dt, libor_curve, libor_curve, first_fixing)
    assert round(v, 4) == -2009695.8385

########################################################################################


def test__fin_float_ois_leg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    leg_pay_rec_type = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swap_float_leg = SwapFloatLeg(
        effective_dt,
        maturity_dt,
        leg_pay_rec_type,
        spread,
        freq_type,
        dc_type,
        notional,
        principal,
        payment_lag,
        cal_type,
        bd_type,
        dg_type,
    )

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    first_fixing = 0.03

    v = swap_float_leg.value(effective_dt, libor_curve, libor_curve, first_fixing)
    assert round(v, 4) == -2038364.5665
