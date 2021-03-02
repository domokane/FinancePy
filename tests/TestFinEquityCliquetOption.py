###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityCliquetOption import FinEquityCliquetOption
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.Date import Date
from financepy.utils.FinGlobalTypes import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityCliquetOption():

    start_date = Date(1, 1, 2014)
    finalExpiryDate = Date(1, 1, 2017)
    freq_type = FinFrequencyTypes.QUARTERLY
    optionType = FinOptionTypes.EUROPEAN_CALL

    cliquetOption = FinEquityCliquetOption(start_date,
                                           finalExpiryDate,
                                           optionType,
                                           freq_type)

    valuation_date = Date(1, 1, 2015)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    model = FinModelBlackScholes(volatility)
    discount_curve = FinDiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = FinDiscountCurveFlat(valuation_date, dividendYield)

    v = cliquetOption.value(valuation_date,
                            stockPrice,
                            discount_curve,
                            dividendCurve,
                            model)

    testCases.header("LABEL", "VALUE")
    testCases.print("FINANCEPY", v)

###############################################################################


test_FinEquityCliquetOption()
testCases.compareTestCases()
