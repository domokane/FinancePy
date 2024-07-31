###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.math import ONE_MILLION
from financepy.utils.global_types import SwapTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date
from financepy.products.rates.swap_fixed_leg import SwapFixedLeg
from financepy.products.rates.swap_float_leg import SwapFloatLeg
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFixedIborSwapLeg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    coupon = -0.44970/100.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFixedLeg = SwapFixedLeg(effective_dt,
                                maturity_dt,
                                legPayRecType,
                                coupon,
                                freq_type,
                                dc_type,
                                notional,
                                principal,
                                payment_lag,
                                cal_type,
                                bd_type,
                                dg_type)

###############################################################################


def test_FinFixedOISSwapLeg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    coupon = -0.515039/100.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFixedLeg = SwapFixedLeg(effective_dt,
                                maturity_dt,
                                legPayRecType,
                                coupon,
                                freq_type,
                                dc_type,
                                notional,
                                principal,
                                payment_lag,
                                cal_type,
                                bd_type,
                                dg_type)

###############################################################################


def test_FinFloatIborLeg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFloatLeg = SwapFloatLeg(effective_dt,
                                maturity_dt,
                                legPayRecType,
                                spread,
                                freq_type,
                                dc_type,
                                notional,
                                principal,
                                payment_lag,
                                cal_type,
                                bd_type,
                                dg_type)

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    firstFixing = 0.03

    v = swapFloatLeg.value(effective_dt,
                           libor_curve,
                           libor_curve,
                           firstFixing)


###############################################################################

def test_FinFloatOISLeg():

    effective_dt = Date(28, 10, 2020)
    maturity_dt = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFloatLeg = SwapFloatLeg(effective_dt,
                                maturity_dt,
                                legPayRecType,
                                spread,
                                freq_type,
                                dc_type,
                                notional,
                                principal,
                                payment_lag,
                                cal_type,
                                bd_type,
                                dg_type)

    libor_curve = DiscountCurveFlat(effective_dt, 0.05)

    firstFixing = 0.03

    v = swapFloatLeg.value(effective_dt, libor_curve, libor_curve,
                           firstFixing)

###############################################################################


def swapFixedLegMonthEnds():

    # Written in response to github issue that has been solved

    fixedleg_1 = SwapFixedLeg(effective_dt=Date(30, 8, 2021),
                              end_dt='2Y',
                              leg_type=SwapTypes.PAY,
                              freq_type=FrequencyTypes.SEMI_ANNUAL,
                              dc_type=DayCountTypes.THIRTY_E_360,
                              cal_type=CalendarTypes.UNITED_STATES,
                              coupon=0.0,
                              end_of_month=False)

    fixedleg_2 = SwapFixedLeg(effective_dt=Date(30, 8, 2021),
                              end_dt='3Y',
                              leg_type=SwapTypes.PAY,
                              freq_type=FrequencyTypes.SEMI_ANNUAL,
                              dc_type=DayCountTypes.THIRTY_E_360,
                              cal_type=CalendarTypes.UNITED_STATES,
                              coupon=0.0,
                              end_of_month=False)

    fixedleg_1.generate_payments()
    fixedleg_2.generate_payments()

    print("leg_1")
    fixedleg_1.print_payments()
    print("leg_2")
    fixedleg_2.print_payments()

###############################################################################


def test_swapFloatLeg():

    effective_dt = Date(1, 9, 2021)

    fixedleg_2 = SwapFixedLeg(effective_dt,
                              end_dt='3Y',
                              leg_type=SwapTypes.PAY,
                              freq_type=FrequencyTypes.SEMI_ANNUAL,
                              dc_type=DayCountTypes.THIRTY_E_360,
                              cal_type=CalendarTypes.UNITED_STATES,
                              coupon=0)

    floatleg_2 = SwapFloatLeg(effective_dt,
                              end_dt='3Y',
                              leg_type=SwapTypes.PAY,
                              freq_type=FrequencyTypes.SEMI_ANNUAL,
                              dc_type=DayCountTypes.THIRTY_E_360,
                              cal_type=CalendarTypes.UNITED_STATES,
                              spread=0)

    fixedleg_2.generate_payments()
    floatleg_2.generate_payment_dts()

    discount_curve = DiscountCurveFlat(
        effective_dt, 0.05, dc_type=DayCountTypes.THIRTY_E_360)
    index_curve = DiscountCurveFlat(
        effective_dt, 0.05, dc_type=DayCountTypes.ACT_ACT_ISDA)

    floatleg_2.value(effective_dt, discount_curve, index_curve)
    # print("leg_2")
    # fixedleg_2.print_payments()
    # print("fleg_2")
    # floatleg_2.print_payments()

###############################################################################


test_swapFloatLeg()
# swapFixedLegMonthEnds()

# Ibor Swap
test_FinFixedIborSwapLeg()
test_FinFloatIborLeg()

# OIS Swap
test_FinFixedOISSwapLeg()
test_FinFloatOISLeg()

test_cases.compareTestCases()
