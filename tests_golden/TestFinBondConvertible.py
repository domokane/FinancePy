###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

# TODO
from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.products.bonds.bond_convertible import BondConvertible
import time
import numpy as np



testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_BondConvertible():

    settlement_date = Date(31, 12, 2003)
    start_convert_date = Date(31, 12, 2003)
    maturity_date = Date(15, 3, 2022)
    conversion_ratio = 38.4615  # adjust for face
    coupon = 0.0575
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrualBasis = DayCountTypes.ACT_365F
    face = 1000.0

    call_price = 1100
    call_dates = [Date(20, 3, 2007),
                  Date(15, 3, 2012),
                  Date(15, 3, 2017)]
    call_prices = np.array([call_price, call_price, call_price])

    putPrice = 90
    put_dates = [Date(20, 3, 2007),
                 Date(15, 3, 2012),
                 Date(15, 3, 2017)]
    put_prices = np.array([putPrice, putPrice, putPrice])

    bond = BondConvertible(maturity_date,
                           coupon,
                           freq_type,
                           start_convert_date,
                           conversion_ratio,
                           call_dates,
                           call_prices,
                           put_dates,
                           put_prices,
                           accrualBasis,
                           face)
#    print(bond)

    dividend_dates = [Date(20, 3, 2007),
                      Date(15, 3, 2008),
                      Date(15, 3, 2009),
                      Date(15, 3, 2010),
                      Date(15, 3, 2011),
                      Date(15, 3, 2012),
                      Date(15, 3, 2013),
                      Date(15, 3, 2014),
                      Date(15, 3, 2015),
                      Date(15, 3, 2016),
                      Date(15, 3, 2017),
                      Date(15, 3, 2018),
                      Date(15, 3, 2019),
                      Date(15, 3, 2020),
                      Date(15, 3, 2021),
                      Date(15, 3, 2022)]

    dividend_yields = [0.00] * 16
    stock_price = 28.5
    stock_volatility = 0.370
    rate = 0.04
    discount_curve = DiscountCurveFlat(settlement_date,
                                       rate,
                                       FrequencyTypes.CONTINUOUS)
    credit_spread = 0.00
    recovery_rate = 0.40
    num_steps_per_year = 20

    testCases.header("LABEL")
    testCases.print("NO CALLS OR PUTS")

    testCases.header("TIME", "NUMSTEPS", "PRICE")

    for num_steps_per_year in [5, 10, 20, 80]:
        start = time.time()
        res = bond.value(settlement_date,
                         stock_price,
                         stock_volatility,
                         dividend_dates,
                         dividend_yields,
                         discount_curve,
                         credit_spread,
                         recovery_rate,
                         num_steps_per_year)

        end = time.time()
        period = end - start
        testCases.print(period, num_steps_per_year, res)

    dividend_yields = [0.02] * 16
    testCases.header("LABEL")
    testCases.print("DIVIDENDS")

    testCases.header("TIME", "NUMSTEPS", "PRICE")
    for num_steps_per_year in [5, 20, 80]:
        start = time.time()
        res = bond.value(settlement_date,
                         stock_price,
                         stock_volatility,
                         dividend_dates,
                         dividend_yields,
                         discount_curve,
                         credit_spread,
                         recovery_rate,
                         num_steps_per_year)
        end = time.time()
        period = end - start
        testCases.print(period, num_steps_per_year, res)

###############################################################################


test_BondConvertible()
testCases.compareTestCases()
