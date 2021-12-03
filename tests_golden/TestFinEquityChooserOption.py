###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_chooser_option import EquityChooserOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_EquityChooserOptionHaug():
    """ Following example in Haug Page 130 """

    valuation_date = Date(1, 1, 2015)
    choose_date = Date(2, 4, 2015)
    call_expiry_date = Date(1, 7, 2015)
    put_expiry_date = Date(2, 8, 2015)
    call_strike = 55.0
    put_strike = 48.0
    stock_price = 50.0
    volatility = 0.35
    interest_rate = 0.10
    dividend_yield = 0.05

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    chooserOption = EquityChooserOption(choose_date,
                                        call_expiry_date,
                                        put_expiry_date,
                                        call_strike,
                                        put_strike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valuation_date,
                                  stock_price,
                                  discount_curve,
                                  dividend_curve,
                                  model, 20000)

    v_haug = 6.0508
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "HAUG", v_haug, "MC", v_mc)

##########################################################################


def test_EquityChooserOptionMatlab():
    """https://fr.mathworks.com/help/fininst/chooserbybls.html """

    valuation_date = Date(1, 6, 2007)
    chooseDate = Date(31, 8, 2007)
    call_expiry_date = Date(2, 12, 2007)
    put_expiry_date = Date(2, 12, 2007)
    call_strike = 60.0
    put_strike = 60.0
    stock_price = 50.0
    volatility = 0.20
    interest_rate = 0.10
    dividend_yield = 0.05

    model = BlackScholes(volatility)

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    chooserOption = EquityChooserOption(chooseDate,
                                        call_expiry_date,
                                        put_expiry_date,
                                        call_strike,
                                        put_strike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valuation_date,
                                  stock_price,
                                  discount_curve,
                                  dividend_curve,
                                  model, 20000)

    v_matlab = 8.9308
    testCases.header("", "", "", "", "", "")
    testCases.print("FINANCEPY", v, "MATLAB", v_matlab, "MC", v_mc)

##########################################################################


def test_EquityChooserOptionDerivicom():
    """http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm """

    valuation_date = Date(1, 1, 2007)
    chooseDate = Date(1, 2, 2007)
    call_expiry_date = Date(1, 4, 2007)
    put_expiry_date = Date(1, 5, 2007)
    call_strike = 40.0
    put_strike = 35.0
    stock_price = 38.0
    volatility = 0.20
    interest_rate = 0.08
    dividend_yield = 0.0625

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    chooserOption = EquityChooserOption(chooseDate,
                                        call_expiry_date,
                                        put_expiry_date,
                                        call_strike,
                                        put_strike)

    v = chooserOption.value(valuation_date,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            model)

    v_mc = chooserOption.value_mc(valuation_date,
                                  stock_price,
                                  discount_curve,
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
