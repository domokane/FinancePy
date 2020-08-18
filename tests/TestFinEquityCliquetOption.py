# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.equity.FinEquityCliquetOption import FinEquityCliquetOption
from financepy.products.equity.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinOptionTypes import FinOptionTypes

import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityCliquetOptionHaug():
    ''' Following example in Haug Page 130 '''

    startDate = FinDate(1, 1, 2015)
    finalExpiryDate = FinDate(1, 1, 2017)
    frequencyType = FinFrequencyTypes.QUARTERLY
    optionType = FinOptionTypes.EUROPEAN_CALL

    cliquetOption = FinEquityCliquetOption(startDate,
                                           finalExpiryDate,
                                           optionType,
                                           frequencyType)

    valueDate = FinDate(1, 8, 2016)
    stockPrice = 50.0
    volatility = 0.35
    interestRate = 0.10
    dividendYield = 0.05
    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    v = cliquetOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividendYield,
                            model)

    cliquetOption.printFlows()
    print("FINANCEPY", v)

###############################################################################


test_FinEquityCliquetOptionHaug()
# test_FinEquityChooserOptionMatlab()
# test_FinEquityChooserOptionDerivicom()
testCases.compareTestCases()
