###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.equity.equity_asian_option import AsianOptionValuationMethods
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

testConvergence = False
testTimeEvolution = False
testMCTimings = True

###############################################################################


def testConvergence():

    valuation_date = Date(1, 1, 2014)
    startAveragingDate = Date(1, 6, 2014)
    expiry_date = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 120  # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    asianOption = EquityAsianOption(startAveragingDate,
                                    expiry_date,
                                    K,
                                    OptionTypes.EUROPEAN_CALL,
                                    num_observations)

    testCases.header(
        "K",
        "Geometric",
        "Turnbull_Wakeman",
        "Curran",
        "FastMC",
        "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    num_paths_list = [5000]

    for num_paths in num_paths_list:

        accruedAverage = stock_price * 1.1

        value_mc_fast = asianOption._value_mc_fast(valuation_date,
                                                   stock_price,
                                                   discount_curve,
                                                   dividend_curve,
                                                   model,
                                                   num_paths,
                                                   seed,
                                                   accruedAverage)

        value_mc_CV = asianOption.value_mc(valuation_date,
                                           stock_price,
                                           discount_curve,
                                           dividend_curve,
                                           model,
                                           num_paths,
                                           seed,
                                           accruedAverage)

        valueGeometric = asianOption.value(valuation_date,
                                           stock_price,
                                           discount_curve,
                                           dividend_curve,
                                           model,
                                           AsianOptionValuationMethods.GEOMETRIC,
                                           accruedAverage)

        valueTurnbullWakeman = asianOption.value(valuation_date,
                                                 stock_price,
                                                 discount_curve,
                                                 dividend_curve,
                                                 model,
                                                 AsianOptionValuationMethods.TURNBULL_WAKEMAN,
                                                 accruedAverage)

        valueCurran = asianOption.value(valuation_date,
                                        stock_price,
                                        discount_curve,
                                        dividend_curve,
                                        model,
                                        AsianOptionValuationMethods.CURRAN,
                                        accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(value_mc_fast)
        valuesMC_CV.append(value_mc_CV)

        testCases.print(
            num_paths,
            valueGeometric,
            valueTurnbullWakeman,
            valueCurran,
            value_mc_fast,
            value_mc_CV)

#    import matplotlib.pyplot as plt
#    x = num_paths_list
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesGeometric,label="Geometric")
#    plt.plot(x,valuesTurnbull,label="Turbull_Wakeman")
#    plt.plot(x,valuesCurran,label="Curran")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_CV,label="MC_CV")
#    plt.legend()
#    plt.xlabel("Number of Paths")
#    plt.show()

###############################################################################


def testTimeEvolution():

    startAveragingDate = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 100  # weekly as we have a year
    accruedAverage = None
    K = 100
    seed = 1976

    model = BlackScholes(volatility)

    asianOption = EquityAsianOption(startAveragingDate,
                                    expiry_date,
                                    K,
                                    OptionTypes.EUROPEAN_CALL,
                                    num_observations)

    testCases.header(
        "Date",
        "Geometric",
        "Turnbull_Wakeman",
        "Curran",
        "FastMC",
        "FastMC_CV")

    valuesTurnbull = []
    valuesCurran = []
    valuesGeometric = []
    valuesMC_fast = []
    valuesMC_CV = []

    valuation_dates = []
    valuation_dates.append(Date(1, 4, 2014))
    valuation_dates.append(Date(1, 6, 2014))
    valuation_dates.append(Date(1, 8, 2014))
    valuation_dates.append(Date(1, 2, 2015))
    valuation_dates.append(Date(1, 4, 2015))
    valuation_dates.append(Date(1, 6, 2015))
    valuation_dates.append(Date(1, 8, 2015))

    num_paths = 10000

    for valuation_date in valuation_dates:

        accruedAverage = stock_price * 0.9

        discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
        dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

        value_mc_fast = asianOption._value_mc_fast(valuation_date,
                                                   stock_price,
                                                   discount_curve,
                                                   dividend_curve,
                                                   model,
                                                   num_paths,
                                                   seed,
                                                   accruedAverage)

        value_mc_CV = asianOption.value_mc(valuation_date,
                                           stock_price,
                                           discount_curve,
                                           dividend_curve,
                                           model,
                                           num_paths,
                                           seed,
                                           accruedAverage)

        valueGeometric = asianOption.value(valuation_date,
                                           stock_price,
                                           discount_curve,
                                           dividend_curve,
                                           model,
                                           AsianOptionValuationMethods.GEOMETRIC,
                                           accruedAverage)

        valueTurnbullWakeman = asianOption.value(valuation_date,
                                                 stock_price,
                                                 discount_curve,
                                                 dividend_curve,
                                                 model,
                                                 AsianOptionValuationMethods.TURNBULL_WAKEMAN,
                                                 accruedAverage)

        valueCurran = asianOption.value(valuation_date,
                                        stock_price,
                                        discount_curve,
                                        dividend_curve,
                                        model,
                                        AsianOptionValuationMethods.CURRAN,
                                        accruedAverage)

        valuesGeometric.append(valueGeometric)
        valuesTurnbull.append(valueTurnbullWakeman)
        valuesCurran.append(valueCurran)
        valuesMC_fast.append(value_mc_fast)
        valuesMC_CV.append(value_mc_CV)

        testCases.print(
            str(valuation_date),
            valueGeometric,
            valueTurnbullWakeman,
            valueCurran,
            value_mc_fast,
            value_mc_CV)

#    import matplotlib.pyplot as plt
#    x = [ dt.date() for dt in valuation_dates]
#
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesGeometric,label="Geometric")
#    plt.plot(x,valuesTurnbull,label="Turbull_Wakeman")
#    plt.plot(x,valuesCurran,label="Curran")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_CV,label="MC_CV")
#    plt.legend()
#    plt.xlabel("Valuation Date")
#    plt.show()

##########################################################################


def testMCTimings():

    valuation_date = Date(1, 1, 2014)
    startAveragingDate = Date(1, 6, 2014)
    expiry_date = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 120  # daily as we have a half year
    accruedAverage = None
    K = 100
    seed = 1976

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    asianOption = EquityAsianOption(startAveragingDate,
                                    expiry_date,
                                    K,
                                    OptionTypes.EUROPEAN_CALL,
                                    num_observations)

    testCases.header(
        "NUMPATHS",
        "VALUE",
        "TIME",
        "VALUE_MC",
        "TIME",
        "VALUE_MC_CV",
        "TIME")

    valuesMC = []
    valuesMC_fast = []
    valuesMC_fast_CV = []

    tvaluesMC = []
    tvaluesMC_fast = []
    tvaluesMC_fast_CV = []

    num_paths_list = [5000]

    for num_paths in num_paths_list:

        accruedAverage = stock_price * 1.1

        start = time.time()
        value_mc = asianOption.value_mc(valuation_date,
                                        stock_price,
                                        discount_curve,
                                        dividend_curve,
                                        model,
                                        num_paths,
                                        seed,
                                        accruedAverage)

        end = time.time()
        t_MC = end - start

        start = time.time()
        value_mc_fast = asianOption._value_mc_fast(valuation_date,
                                                   stock_price,
                                                   discount_curve,
                                                   dividend_curve,
                                                   model,
                                                   num_paths,
                                                   seed,
                                                   accruedAverage)

        end = time.time()
        t_MC_fast = end - start

        start = time.time()
        value_mc_fast_CV = asianOption.value_mc(valuation_date,
                                                stock_price,
                                                discount_curve,
                                                dividend_curve,
                                                model,
                                                num_paths,
                                                seed,
                                                accruedAverage)

        end = time.time()
        t_MC_fast_CV = end - start

        valuesMC.append(value_mc)
        valuesMC_fast.append(value_mc_fast)
        valuesMC_fast_CV.append(value_mc_fast_CV)

        tvaluesMC.append(t_MC)
        tvaluesMC_fast.append(t_MC_fast)
        tvaluesMC_fast_CV.append(t_MC_fast_CV)

        testCases.print(
            num_paths,
            value_mc,
            t_MC,
            value_mc_fast,
            t_MC_fast,
            value_mc_fast_CV,
            t_MC_fast_CV)

#    import matplotlib.pyplot as plt
#    x = num_paths_list
#    plt.figure(figsize=(8,6))
#    plt.plot(x,valuesMC,label="Basic MC")
#    plt.plot(x,valuesMC_fast,label="MC_Fast")
#    plt.plot(x,valuesMC_fast_CV,label="MC_Fast CV")
#    plt.legend()
#    plt.xlabel("Number of Paths")
#    plt.show()


testConvergence()
testMCTimings()
testTimeEvolution()
testCases.compareTestCases()
