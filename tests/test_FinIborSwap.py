###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from helpers import buildIborSingleCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.utils.math import ONE_MILLION
import numpy as np


def test_LiborSwap():
    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve
    start_date = Date(27, 12, 2017)
    end_date = Date(27, 12, 2067)

    fixed_coupon = 0.015
    fixedFreqType = FrequencyTypes.ANNUAL
    fixed_day_count_type = DayCountTypes.THIRTY_E_360

    float_spread = 0.0
    floatFreqType = FrequencyTypes.SEMI_ANNUAL
    float_day_count_type = DayCountTypes.ACT_360
    firstFixing = -0.00268

    swapCalendarType = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    fixed_leg_type = SwapTypes.RECEIVE

    notional = 10.0 * ONE_MILLION

    swap = IborSwap(start_date,
                    end_date,
                    fixed_leg_type,
                    fixed_coupon,
                    fixedFreqType,
                    fixed_day_count_type,
                    notional,
                    float_spread,
                    floatFreqType,
                    float_day_count_type,
                    swapCalendarType,
                    bus_day_adjust_type,
                    date_gen_rule_type)

    """ Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. """

    valuation_date = Date(30, 11, 2018)
    settlement_date = valuation_date.add_days(2)
    libor_curve = buildIborSingleCurve(valuation_date)
    v = swap.value(settlement_date, libor_curve, libor_curve, firstFixing)

    assert round(v, 4) == 318901.6015


def test_dp_example():
    #  http://www.derivativepricing.com/blogpage.asp?id=8

    start_date = Date(14, 11, 2011)
    end_date = Date(14, 11, 2016)
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL
    swapCalendarType = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    fixed_day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_leg_type = SwapTypes.PAY
    fixed_coupon = 0.0124
    notional = ONE_MILLION

    swap = IborSwap(start_date,
                    end_date,
                    fixed_leg_type,
                    fixed_coupon=fixed_coupon,
                    fixed_freq_type=fixedFreqType,
                    fixed_day_count_type=fixed_day_count_type,
                    float_freq_type=FrequencyTypes.SEMI_ANNUAL,
                    float_day_count_type=DayCountTypes.ACT_360,
                    notional=notional,
                    calendar_type=swapCalendarType,
                    bus_day_adjust_type=bus_day_adjust_type,
                    date_gen_rule_type=date_gen_rule_type)

    dts = [Date(14, 11, 2011), Date(14, 5, 2012), Date(14, 11, 2012),
           Date(14, 5, 2013), Date(14, 11, 2013), Date(14, 5, 2014),
           Date(14, 11, 2014), Date(14, 5, 2015), Date(16, 11, 2015),
           Date(16, 5, 2016), Date(14, 11, 2016)]

    dfs = [0.9999843, 0.9966889, 0.9942107, 0.9911884, 0.9880738, 0.9836490,
           0.9786276, 0.9710461, 0.9621778, 0.9514315, 0.9394919]

    valuation_date = start_date

    curve = DiscountCurve(valuation_date, dts, np.array(dfs),
                          InterpTypes.FLAT_FWD_RATES)

    v = swap.value(valuation_date, curve, curve)

    # This is essentially zero
    assert round(v * 1000, 4) == 785300.0566
