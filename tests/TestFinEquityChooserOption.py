###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_chooser_option import EquityChooserOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_EquityChooserOptionHaug():
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

    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    chooserOption = EquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividend_curve,
                                 model, 20000)

    v_haug = 6.0508
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "HAUG", v_haug, "MC", v_mc)

##########################################################################


def test_EquityChooserOptionMatlab():
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

    model = BlackScholes(volatility)

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    chooserOption = EquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividend_curve,
                                 model, 20000)

    v_matlab = 8.9308
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "MATLAB", v_matlab, "MC", v_mc)

##########################################################################


def test_EquityChooserOptionDerivicom():
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

    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    chooserOption = EquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valueDate,
                            stockPrice,
                            discountCurve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valueDate,
                                 stockPrice,
                                 discountCurve,
                                 dividend_curve,
                                 model, 20000)

    v_derivicom = 1.0989
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "DERIVICOM", v_derivicom, "MC", v_mc)

##########################################################################


test_EquityChooserOptionHaug()
test_EquityChooserOptionMatlab()
test_EquityChooserOptionDerivicom()
testCases.compareTestCases()
