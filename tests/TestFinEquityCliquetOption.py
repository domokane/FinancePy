###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityCliquetOption():

    startDate = FinDate(1, 1, 2014)
    finalExpiryDate = FinDate(1, 1, 2017)
    freqType = FinFrequencyTypes.QUARTERLY
    option_type = FinOptionTypes.EUROPEAN_CALL

    cliquetOption = EquityCliquetOption(startDate,
                                           finalExpiryDate,
                                           option_type,
                                           freqType)

    valueDate = FinDate(1, 1, 2015)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    v = cliquetOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividend_curve,
                            model)

    testCases.header("LABEL", "VALUE")
    testCases.print("FINANCEPY", v)

###############################################################################


test_EquityCliquetOption()
testCases.compareTestCases()
