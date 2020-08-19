# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.equity.FinEquityChooserOption import FinEquityChooserOption
from financepy.products.equity.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityChooserOptionHaug():
    ''' Following example in Haug Page 130 '''

    valueDate = FinDate(1, 1, 2015)
    chooseDate = FinDate(2, 4, 2015)
    callExpiryDate = FinDate(1, 7, 2015)
    putExpiryDate = FinDate(2, 8, 2015)
    callStrike = 55.0
    putStrike = 48.0
    stockPrice = 50.0
    volatility = 0.35
    interestRate = 0.10
    dividendYield = 0.05

    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividendYield,
                            model)

    v_mc = chooserOption.valueMC(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividendYield,
                                 model, 20000)

    v_haug = 6.0508
    print("FINANCEPY", v, "HAUG", v_haug, "MC", v_mc)

##########################################################################


def test_FinEquityChooserOptionMatlab():
    '''https://fr.mathworks.com/help/fininst/chooserbybls.html '''

    valueDate = FinDate(1, 6, 2007)
    chooseDate = FinDate(31, 8, 2007)
    callExpiryDate = FinDate(2, 12, 2007)
    putExpiryDate = FinDate(2, 12, 2007)
    callStrike = 60.0
    putStrike = 60.0
    stockPrice = 50.0
    volatility = 0.20
    interestRate = 0.10
    dividendYield = 0.05

    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividendYield,
                            model)

    v_mc = chooserOption.valueMC(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividendYield,
                                 model, 20000)

    v_matlab = 8.9308
    print("FINANCEPY", v, "MATLAB", v_matlab, "MC", v_mc)

##########################################################################


def test_FinEquityChooserOptionDerivicom():
    '''http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm '''

    valueDate = FinDate(1, 1, 2007)
    chooseDate = FinDate(1, 2, 2007)
    callExpiryDate = FinDate(1, 4, 2007)
    putExpiryDate = FinDate(1, 5, 2007)
    callStrike = 40.0
    putStrike = 35.0
    stockPrice = 38.0
    volatility = 0.20
    interestRate = 0.08
    dividendYield = 0.0625

    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividendYield,
                            model)

    v_mc = chooserOption.valueMC(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividendYield,
                                 model, 20000)

    v_derivicom = 1.0989
    print("FINANCEPY", v, "DERIVICOM", v_derivicom, "MC", v_mc)

##########################################################################


test_FinEquityChooserOptionHaug()
test_FinEquityChooserOptionMatlab()
test_FinEquityChooserOptionDerivicom()
testCases.compareTestCases()
