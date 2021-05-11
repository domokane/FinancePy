###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.curve_flat import DiscountCurveFlat
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
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFixedIborSwapLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)

    coupon = -0.44970/100.0
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFixedLeg = SwapFixedLeg(effective_date,
                                maturity_date,
                                legPayRecType,
                                coupon,
                                freq_type,
                                day_count_type,
                                notional,
                                principal,
                                payment_lag,
                                calendar_type,
                                bus_day_adjust_type,
                                date_gen_rule_type)

###############################################################################


def test_FinFixedOISSwapLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)

    coupon = -0.515039/100.0
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFixedLeg = SwapFixedLeg(effective_date,
                                maturity_date,
                                legPayRecType,
                                coupon,
                                freq_type,
                                day_count_type,
                                notional,
                                principal,
                                payment_lag,
                                calendar_type,
                                bus_day_adjust_type,
                                date_gen_rule_type)

###############################################################################


def test_FinFloatIborLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFloatLeg = SwapFloatLeg(effective_date,
                                maturity_date,
                                legPayRecType,
                                spread,
                                freq_type,
                                day_count_type,
                                notional,
                                principal,
                                payment_lag,
                                calendar_type,
                                bus_day_adjust_type,
                                date_gen_rule_type)

    libor_curve = DiscountCurveFlat(effective_date, 0.05)

    firstFixing = 0.03

    v = swapFloatLeg.value(effective_date, libor_curve, libor_curve,
                           firstFixing)


###############################################################################

def test_FinFloatOISLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)

    spread = 0.0
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = SwapTypes.PAY
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFloatLeg = SwapFloatLeg(effective_date,
                                maturity_date,
                                legPayRecType,
                                spread,
                                freq_type,
                                day_count_type,
                                notional,
                                principal,
                                payment_lag,
                                calendar_type,
                                bus_day_adjust_type,
                                date_gen_rule_type)

    libor_curve = DiscountCurveFlat(effective_date, 0.05)

    firstFixing = 0.03

    v = swapFloatLeg.value(effective_date, libor_curve, libor_curve,
                           firstFixing)

###############################################################################


# Ibor Swap
test_FinFixedIborSwapLeg()
test_FinFloatIborLeg()

# OIS Swap
test_FinFixedOISSwapLeg()
test_FinFloatOISLeg()

testCases.compareTestCases()
