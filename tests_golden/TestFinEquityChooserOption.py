###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityChooserOption import FinEquityChooserOption
from financepy.models.black_scholes import FinModelBlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinEquityChooserOptionHaug():
    """ Following example in Haug Page 130 """

    valuation_date = Date(1, 1, 2015)
    chooseDate = Date(2, 4, 2015)
    callExpiryDate = Date(1, 7, 2015)
    putExpiryDate = Date(2, 8, 2015)
    callStrike = 55.0
    putStrike = 48.0
    stock_price = 50.0
    volatility = 0.35
    interestRate = 0.10
    dividendYield = 0.05

    model = FinModelBlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividendCurve,
                            model)

    v_mc = chooserOption.valueMC(valuation_date,
                                 stock_price,
                                 discount_curve,
                                 dividendCurve,
                                 model, 20000)

    v_haug = 6.0508
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "HAUG", v_haug, "MC", v_mc)

##########################################################################


def test_FinEquityChooserOptionMatlab():
    """https://fr.mathworks.com/help/fininst/chooserbybls.html """

    valuation_date = Date(1, 6, 2007)
    chooseDate = Date(31, 8, 2007)
    callExpiryDate = Date(2, 12, 2007)
    putExpiryDate = Date(2, 12, 2007)
    callStrike = 60.0
    putStrike = 60.0
    stock_price = 50.0
    volatility = 0.20
    interestRate = 0.10
    dividendYield = 0.05

    model = FinModelBlackScholes(volatility)

    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividendCurve,
                            model)

    v_mc = chooserOption.valueMC(valuation_date,
                                 stock_price,
                                 discount_curve,
                                 dividendCurve,
                                 model, 20000)

    v_matlab = 8.9308
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "MATLAB", v_matlab, "MC", v_mc)

##########################################################################


def test_FinEquityChooserOptionDerivicom():
    """http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm """

    valuation_date = Date(1, 1, 2007)
    chooseDate = Date(1, 2, 2007)
    callExpiryDate = Date(1, 4, 2007)
    putExpiryDate = Date(1, 5, 2007)
    callStrike = 40.0
    putStrike = 35.0
    stock_price = 38.0
    volatility = 0.20
    interestRate = 0.08
    dividendYield = 0.0625

    model = FinModelBlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    chooserOption = FinEquityChooserOption(chooseDate,
                                           callExpiryDate,
                                           putExpiryDate,
                                           callStrike,
                                           putStrike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividendCurve,
                            model)

    v_mc = chooserOption.valueMC(valuation_date,
                                 stock_price,
                                 discount_curve,
                                 dividendCurve,
                                 model, 20000)

    v_derivicom = 1.0989
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "DERIVICOM", v_derivicom, "MC", v_mc)

##########################################################################


test_FinEquityChooserOptionHaug()
test_FinEquityChooserOptionMatlab()
test_FinEquityChooserOptionDerivicom()
testCases.compareTestCases()
