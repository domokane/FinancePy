###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# TODO
import time
import numpy as np

import sys
sys.path.append("..")

from financepy.products.bonds.FinBondConvertible import FinBondConvertible
from financepy.utils.Date import Date
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBondConvertible():

    settlement_date = Date(31, 12, 2003)
    startConvertDate = Date(31, 12, 2003)
    maturity_date = Date(15, 3, 2022)
    conversionRatio = 38.4615  # adjust for face
    coupon = 0.0575
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrualBasis = FinDayCountTypes.ACT_365F
    face = 1000.0

    callPrice = 1100
    callDates = [Date(20, 3, 2007),
                 Date(15, 3, 2012),
                 Date(15, 3, 2017)]
    callPrices = np.array([callPrice, callPrice, callPrice])

    putPrice = 90
    putDates = [Date(20, 3, 2007),
                Date(15, 3, 2012),
                Date(15, 3, 2017)]
    putPrices = np.array([putPrice, putPrice, putPrice])

    bond = FinBondConvertible(maturity_date,
                              coupon,
                              freq_type,
                              startConvertDate,
                              conversionRatio,
                              callDates,
                              callPrices,
                              putDates,
                              putPrices,
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

    dividendYields = [0.00] * 16
    stockPrice = 28.5
    stockVolatility = 0.370
    rate = 0.04
    discount_curve = FinDiscountCurveFlat(settlement_date,
                                         rate,
                                         FinFrequencyTypes.CONTINUOUS)
    creditSpread = 0.00
    recovery_rate = 0.40
    num_steps_per_year = 20

    testCases.header("LABEL")
    testCases.print("NO CALLS OR PUTS")

    testCases.header("TIME", "NUMSTEPS", "PRICE")

    for num_steps_per_year in [5, 10, 20, 80]:
        start = time.time()
        res = bond.value(settlement_date,
                         stockPrice,
                         stockVolatility,
                         dividend_dates,
                         dividendYields,
                         discount_curve,
                         creditSpread,
                         recovery_rate,
                         num_steps_per_year)

        end = time.time()
        period = end - start
        testCases.print(period, num_steps_per_year, res)

    dividendYields = [0.02] * 16
    testCases.header("LABEL")
    testCases.print("DIVIDENDS")

    testCases.header("TIME", "NUMSTEPS", "PRICE")
    for num_steps_per_year in [5, 20, 80]:
        start = time.time()
        res = bond.value(settlement_date,
                         stockPrice,
                         stockVolatility,
                         dividend_dates,
                         dividendYields,
                         discount_curve,
                         creditSpread,
                         recovery_rate,
                         num_steps_per_year)
        end = time.time()
        period = end - start
        testCases.print(period, num_steps_per_year, res)

###############################################################################


test_FinBondConvertible()
testCases.compareTestCases()
