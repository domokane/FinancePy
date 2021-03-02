###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Math import ONE_MILLION
from financepy.utils.FinGlobalTypes import FinSwapTypes
from financepy.utils.Calendar import FinBusDayAdjustTypes
from financepy.utils.Calendar import FinDateGenRuleTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.Amount import FinAmount
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.Calendar import FinCalendarTypes
from financepy.utils.Date import Date
from financepy.products.rates.FinFixedLeg import FinFixedLeg
from financepy.products.rates.FinFloatLeg import FinFloatLeg
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFixedIborSwapLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)
    
    coupon = -0.44970/100.0
    freq_type = FinFrequencyTypes.ANNUAL
    day_count_type = FinDayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = FinSwapTypes.PAY
    calendar_type = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFixedLeg = FinFixedLeg(effective_date,
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
    freq_type = FinFrequencyTypes.ANNUAL
    day_count_type = FinDayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = FinSwapTypes.PAY
    calendar_type = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFixedLeg = FinFixedLeg(effective_date,
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
    freq_type = FinFrequencyTypes.ANNUAL
    day_count_type = FinDayCountTypes.THIRTY_360_BOND
    notional = 10.0 * ONE_MILLION
    legPayRecType = FinSwapTypes.PAY
    calendar_type = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    payment_lag = 0
    principal = 0.0

    swapFloatLeg = FinFloatLeg(effective_date,
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

    libor_curve = FinDiscountCurveFlat(effective_date, 0.05)

    firstFixing = 0.03

    v = swapFloatLeg.value(effective_date, libor_curve, libor_curve,
                           firstFixing)


###############################################################################

def test_FinFloatOISLeg():

    effective_date = Date(28, 10, 2020)
    maturity_date = Date(28, 10, 2025)
    
    spread = 0.0
    freq_type = FinFrequencyTypes.ANNUAL
    day_count_type = FinDayCountTypes.ACT_360
    notional = 10.0 * ONE_MILLION
    legPayRecType = FinSwapTypes.PAY
    calendar_type = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    payment_lag = 1
    principal = 0.0

    swapFloatLeg = FinFloatLeg(effective_date,
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

    libor_curve = FinDiscountCurveFlat(effective_date, 0.05)

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
