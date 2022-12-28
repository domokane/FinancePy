###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")


import time
from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def testEquityAmericanOption():

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    stock_price = 50.0
    interest_rate = 0.06
    dividend_yield = 0.04
    volatility = 0.40
    strike_price = 50.0

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    testCases.banner("================== EUROPEAN PUT =======================")

    put_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_PUT)

    num_steps = 100

    model = BlackScholes(volatility,
                         BlackScholesTypes.CRR_TREE,
                         num_steps)

    value = put_option.value(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    delta = put_option.delta(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    gamma = put_option.gamma(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    theta = put_option.theta(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS",
                     "VALUE DELTA GAMMA THETA", "TIME")

    num_steps_list = [100]

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price,
                               discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", num_steps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price,
                               discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", num_steps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    call_option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.EUROPEAN_CALL)

    value = call_option.value(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    delta = call_option.delta(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    gamma = call_option.gamma(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    theta = call_option.theta(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)
        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_CALL_TREE", num_steps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")
    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    option = EquityAmericanOption(expiry_date, strike_price,
                                  OptionTypes.AMERICAN_CALL)

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_TREE_CALL", num_steps, results, duration)

        ############### NO DO IT USING AMERICAN MONTE CARLO 
        
        num_paths = 50000

        model = BlackScholes(volatility,
                             BlackScholesTypes.LSMC,
                             num_steps, 
                             num_paths)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_LSMC_CALL", num_steps, results, duration)
 
    testCases.banner(
        "================== AMERICAN PUT =======================")
    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    option = EquityAmericanOption(expiry_date, strike_price,
                                  OptionTypes.AMERICAN_PUT)

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_TREE_PUT", num_steps, results, duration)
#        print("PUT Tree Results", results)

        ############### NO DO IT USING AMERICAN MONTE CARLO 

        model = BlackScholes(volatility,
                             BlackScholesTypes.LSMC,
                             num_steps, 
                             num_paths)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", num_steps, results, duration)
#        print("PUT LSMC Results", results)

#    FinTest.TestReport(filename)

###############################################################################

from financepy.models.equity_lsmc import equity_lsmc
from financepy.models.black_scholes_analytic import bs_value
from financepy.models.equity_crr_tree import crr_tree_val

def replicateLSPaper():

    amer_option_call_value = OptionTypes.AMERICAN_CALL.value
    amer_option_put_value = OptionTypes.AMERICAN_PUT.value

    option_type_values = [amer_option_call_value, amer_option_put_value]    
    option_type_values = [amer_option_put_value]    
    stock_prices = [36.0, 38.0, 40.0, 42.0, 44.0]
    volatilities = [0.20, 0.40]
    times_to_expiry = [1.0, 2.0]

    num_paths = 50000 
    num_steps_per_year = 10

    r = 0.06
    q = 0.0
    k = 40.0

    use_sobol = False
    seed = 1912 

    print("   S     v    T    v_tree  v_eur   v_ls")

    for option_type_value in option_type_values:
        for s in stock_prices:
            for v in volatilities:
                for t in times_to_expiry:
    
                    v_ls = equity_lsmc(s,
                                       r,  # continuously compounded
                                       q,  # continuously compounded
                                       v,  # Black scholes volatility
                                       50, # Steps per year
                                       num_paths,
                                       t,
                                       option_type_value,
                                       k,
                                       use_sobol,
                                       seed)
    
                    v_tree = crr_tree_val(s,
                                          r,  # continuously compounded
                                          q,  # continuously compounded
                                          v,  # Black scholes volatility
                                          2000,
                                          t,
                                          option_type_value,
                                          k,
                                          True)[0]
        
                    if option_type_value == OptionTypes.AMERICAN_CALL.value:
                        euro_option_type_value = OptionTypes.EUROPEAN_CALL.value
                    else:
                        euro_option_type_value = OptionTypes.EUROPEAN_PUT.value

                    v_eur = bs_value(s, t, k, r, q, v, euro_option_type_value)

                    print("%5.1f %5.2f %4.1f %7.3f %7.3f %7.3f" % 
                          (s, v, t, v_tree, v_eur, v_ls))
        
###############################################################################

# replicateLSPaper()
testEquityAmericanOption()
testCases.compareTestCases()
