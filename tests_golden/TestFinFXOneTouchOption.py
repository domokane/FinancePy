###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.global_types import TouchOptionTypes
from financepy.products.fx import FXOneTouchOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFXOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(2, 7, 2016)
    volatility = 0.20
    barrier_level = 1.0  # H
    model = BlackScholes(volatility)

    domesticRate = 0.10
    foreignRate = 0.03

    num_paths = 50000
    num_steps_per_year = 252 * 2

    domCurve = DiscountCurveFlat(valuation_date, domesticRate)
    forCurve = DiscountCurveFlat(valuation_date, foreignRate)

    spot_fx_rate = 1.050
    payment_size = 1.5

    testCases.header("================================= CASH ONLY")

    downTypes = [TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT,
                 TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for downType in downTypes:

        option = FXOneTouchOption(expiry_date,
                                  downType,
                                  barrier_level,
                                  payment_size)

        v = option.value(valuation_date,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

        print(downType, v, v_mc)

    spot_fx_rate = 0.950
    payment_size = 1.5

    upTypes = [TouchOptionTypes.UP_AND_IN_CASH_AT_HIT,
               TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
               TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for upType in upTypes:

        option = FXOneTouchOption(expiry_date,
                                  upType,
                                  barrier_level,
                                  payment_size)

        v = option.value(valuation_date,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

        print(upType, v, v_mc)


###############################################################################

def test_BBGOneTouchOption():

    # 1YR ONETOUCH ON EURUSD 

    valuation_date = Date(3, 12, 2021)
    expiry_date = Date(5, 12, 2022)
    barrier_level = 1.1865  # THIS IS NUMBER OF DOLLARS PER EURO

    spot_fx_rate = 1.1300 # EURUSD
    volatility = 0.06075

    model = BlackScholes(volatility)

    forRate = 0.00593 # EUR
    domRate = -0.00414 # USD

    num_paths = 50000
    num_steps_per_year = 252

    domCurve = DiscountCurveFlat(valuation_date, domRate)
    forCurve = DiscountCurveFlat(valuation_date, forRate)

    payment_size = 1000000 # EUR

    optionType = TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY

    option = FXOneTouchOption(expiry_date,
                              optionType,
                              barrier_level,
                              payment_size)

    v = option.value(valuation_date,
                     spot_fx_rate,
                     domCurve,
                     forCurve,
                     model)

    v_mc = option.value_mc(valuation_date,
                           spot_fx_rate,
                           domCurve,
                           forCurve,
                           model,
                           num_steps_per_year,
                           num_paths)

    d = option.delta(valuation_date,
                     spot_fx_rate,
                     domCurve,
                     forCurve,
                     model)

    g = option.gamma(valuation_date,
                     spot_fx_rate,
                     domCurve,
                     forCurve,
                     model)

    v = option.vega(valuation_date,
                     spot_fx_rate,
                     domCurve,
                     forCurve,
                     model)

    # I SHOULD GET 49.4934% OR 494,934 in EUR
    # VEGA IS 68,777.26
    # GAMMA IS 916,285
    # DELTA IS -9560266
    
    print(optionType)
    print("Value:", v)
    print("Value MC:", v_mc)
    print("Delta: ", d)
    print("Gamma:", g)
    print("Vega:", v)

###############################################################################


test_FinFXOneTouchOption()
test_BBGOneTouchOption()
testCases.compareTestCases()
