###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_types import SwapTypes, ReturnTypes
from financepy.utils.math import ONE_MILLION

from financepy.products.equity.equity_swap import EquitySwap
from financepy.products.equity.equity_swap_leg import EquitySwapLeg

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

testCases = FinTestCases(__file__, globalTestCaseMode)

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

def test_equity_swap_not_in_inception():

    ## According to http://www-2.rotman.utoronto.ca/~hull/technicalnotes/TechnicalNote19.pdf
    ## We can engineer a price to which the equity and float leg balance each other. This is
    ## relatively easy for a single period swap.

    effective_date = Date(13, 2, 2018)
    valuation_date = effective_date.add_months(6)
    maturity_date = effective_date.add_months(12)

    leg_type = SwapTypes.RECEIVE
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.THIRTY_360_BOND
    return_type = ReturnTypes.TOTAL_RETURN
    notional = ONE_MILLION
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    payment_lag = 0
    rate_spread = 0

    stock_strike = 125.0 
    stock_qty = notional/stock_strike
    discountRate = 0.05
    dividendRate = 0.00

    discount_curve = DiscountCurveFlat(valuation_date, discountRate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividendRate)
    
    index_curve = discount_curve

    ## Rate determined at last reset date, from that date to maturity
    index_curve_first = DiscountCurveFlat(effective_date, discountRate)
    index_alpha_first = DayCount(index_curve_first._day_count_type).year_frac(effective_date, maturity_date)[0]
    firstFixing = ((index_curve_first.df(effective_date) / index_curve_first.df(maturity_date))  - 1.0 ) / index_alpha_first

    ## Rate between valuation date to maturity
    index_curve_period = DiscountCurveFlat(valuation_date, discountRate)
    index_alpha_period = DayCount(index_curve_period._day_count_type).year_frac(valuation_date, maturity_date)[0]
    periodFixing = ((index_curve_period.df(valuation_date) / index_curve_period.df(maturity_date))  - 1.0 ) / index_alpha_period

    ## This is the price at which abs_value(equity leg) == abs_value(float leg)
    stock_price = stock_strike * (1 + firstFixing * index_alpha_first) / (1 + periodFixing * index_alpha_period)

    equity_swap = EquitySwap(effective_date,
                             maturity_date,
                             leg_type,
                             freq_type,
                             day_count_type,
                             stock_strike,
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
    
    value = equity_swap.value(valuation_date, 
                              discount_curve,
                              index_curve,
                              dividend_curve,
                              stock_price,
                              firstFixing)
    

###############################################################################

def test_equity_swap_with_dividends():

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
    dividendRate = 0.02
    indexRate = 0.03

    discount_curve = DiscountCurveFlat(effective_date, discountRate)
    dividend_curve = DiscountCurveFlat(effective_date, dividendRate)
    index_curve    = DiscountCurveFlat(effective_date, indexRate)

    equity_swap_leg = SwapEquityLeg(effective_date,
                                    maturity_date,
                                    leg_type,
                                    freq_type,
                                    day_count_type,
                                    stock_price,
                                    stock_qty,
                                    payment_lag,
                                    return_type,
                                    calendar_type,
                                    bus_day_adjust_type,
                                    date_gen_rule_type)
    
    value_higher_disc = equity_swap_leg.value(effective_date, 
                                              discount_curve,
                                              discount_curve,)
    
    value_with_divs = equity_swap_leg.value(effective_date, 
                                            discount_curve,
                                            index_curve,
                                            dividend_curve,
                                            stock_price,)
    
