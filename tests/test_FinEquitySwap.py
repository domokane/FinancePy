###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_swap import EquitySwap
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_types import SwapTypes, ReturnTypes
from financepy.utils.math import ONE_MILLION

def test_equity_swap_at_inception():

    effective_date = Date(13, 2, 2018)
    maturity_date = effective_date.add_months(24)

    leg_type = SwapTypes.RECEIVE
    freq_type = FrequencyTypes.SEMI_ANNUAL
    day_count_type = DayCountTypes.THIRTY_360_BOND
    return_type = ReturnTypes.TOTAL_RETURN
    notional = ONE_MILLION
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    rate_spread = 0

    stock_price = 130.0
    stock_qty = notional/stock_price
    discountRate = 0.05
    dividendRate = 0.00

    discount_curve = DiscountCurveFlat(effective_date, discountRate)
    dividend_curve = DiscountCurveFlat(effective_date, dividendRate)

    index_curve = discount_curve

    equity_swap = EquitySwap(effective_date,
                                maturity_date,
                                leg_type,
                                freq_type,
                                day_count_type,
                                stock_price,
                                stock_qty,
                                payment_lag,
                                return_type,
                                freq_type,
                                day_count_type,
                                rate_spread,
                                payment_lag,
                                calendar_type,
                                bus_day_adjust_type,
                                date_gen_rule_type)
    
    value = equity_swap.value(effective_date, 
                                discount_curve,
                                index_curve,
                                dividend_curve)
    
    assert round(value, 5) == 0.00000
