###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes

from financepy.utils.global_types import FinOptionTypes

from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.models.black_scholes import FinModelBlackScholesTypes

from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_american_option import EquityAmericanOption

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

# TODO Complete output of results to log files

def testFinModelBlackScholes():

    valuation_date = Date(8, 5, 2015)
    expiry_date = Date(15, 1, 2016)

    strike_price = 130.0
    stock_price = 127.62
    volatility = 0.20
    interest_rate = 0.001
    dividendYield = 0.0163

    option_type = FinOptionTypes.AMERICAN_CALL
    euOptionType = FinOptionTypes.EUROPEAN_CALL
    
    amOption = EquityAmericanOption(expiry_date, strike_price,
                                       option_type)
    
    ameuOption = EquityAmericanOption(expiry_date, strike_price,
                                         euOptionType)
    
    euOption = EquityVanillaOption(expiry_date, strike_price,
                                      euOptionType)
    
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate,
                                       FrequencyTypes.CONTINUOUS,
                                       DayCountTypes.ACT_365F)

    dividend_curve = DiscountCurveFlat(valuation_date, dividendYield,
                                      FrequencyTypes.CONTINUOUS,
                                      DayCountTypes.ACT_365F)
    
    num_steps_per_year = 400
    
    modelTree = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             num_steps_per_year)
    
    v = amOption.value(valuation_date, stock_price, discount_curve,
                           dividend_curve, modelTree)
#    print(v)

    modelApprox = BlackScholes(volatility,
                               FinModelBlackScholesTypes.BARONE_ADESI)

    v = amOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, modelApprox)

#    print(v)

    v = ameuOption.value(valuation_date, stock_price, discount_curve,
                           dividend_curve, modelTree)

#    print(v)

    v = euOption.value(valuation_date, stock_price, discount_curve,
                         dividend_curve, modelTree)

#    print(v)

    amTreeValue = []
    amBAWValue = []
    euTreeValue = []
    euAnalValue = []
    volatility = 0.20

    # num_steps_per_year = range(5, 200, 1)
    
    # for num_steps in num_steps_per_year:

    #     modelTree = FinModelBlackScholes(volatility,
    #                                      FinModelBlackScholesTypes.CRR_TREE,
    #                                      {'num_steps_per_year':num_steps})

    #     modelAnal = FinModelBlackScholes(volatility, 
    #                                      FinModelBlackScholesTypes.ANALYTICAL)

    #     modelBAW = FinModelBlackScholes(volatility, 
    #                                     FinModelBlackScholesTypes.BARONE_ADESI)


    #     v_am = amOption.value(valuation_date, stock_price, discount_curve,
    #                           dividendYield, modelTree)

    #     v_eu = ameuOption.value(valuation_date, stock_price, discount_curve,
    #                             dividendYield, modelTree)
 
    #     v_bs = euOption.value(valuation_date, stock_price, discount_curve,
    #                           dividendYield, modelAnal)

    #     v_am_baw = amOption.value(valuation_date, stock_price, discount_curve,
    #                               dividendYield, modelBAW)
        
    #     amTreeValue.append(v_am)
    #     euTreeValue.append(v_eu)
    #     euAnalValue.append(v_bs)
    #     amBAWValue.append(v_am_baw)
        
    
    # plt.title("American PUT Option Price Convergence Analysis")
    # plt.plot(num_steps_per_year, amTreeValue, label="American Tree")
    # plt.plot(num_steps_per_year, amBAWValue, label="American BAW")
    # plt.plot(num_steps_per_year, euTreeValue, label="European Tree")
    # plt.plot(num_steps_per_year, euAnalValue, label="European Anal", lw =2)
    # plt.xlabel("Num Steps")
    # plt.ylabel("Value")
    # plt.legend();

###############################################################################


testFinModelBlackScholes()
testCases.compareTestCases()
