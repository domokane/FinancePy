###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityCliquetOption import FinEquityCliquetOption
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityCliquetOptionHaug():
    ''' Following example in Haug Page 130 '''

    startDate = FinDate(1, 1, 2015)
    finalExpiryDate = FinDate(1, 1, 2017)
    freqType = FinFrequencyTypes.QUARTERLY
    optionType = FinOptionTypes.EUROPEAN_CALL

    cliquetOption = FinEquityCliquetOption(startDate,
                                           finalExpiryDate,
                                           optionType,
                                           freqType)

    valueDate = FinDate(1, 8, 2016)
    stockPrice = 50.0
    volatility = 0.35
    interestRate = 0.10
    dividendYield = 0.05
    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    v = cliquetOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividendYield,
                            model)

    testCases.header("LABEL", "VALUE")
    testCases.print("FINANCEPY", v)

###############################################################################


test_FinEquityCliquetOptionHaug()
testCases.compareTestCases()
