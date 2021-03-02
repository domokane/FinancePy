###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Date import Date
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes

from financepy.utils.FinGlobalTypes import FinOptionTypes

from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.models.FinModelBlackScholes import FinModelBlackScholesTypes

from financepy.products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from financepy.products.equity.FinEquityAmericanOption import FinEquityAmericanOption

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

# TODO Complete output of results to log files

def testFinModelBlackScholes():

    valuation_date = Date(8, 5, 2015)
    expiry_date = Date(15, 1, 2016)

    strikePrice = 130.0
    stockPrice = 127.62
    volatility = 0.20
    interestRate = 0.001
    dividendYield = 0.0163

    optionType = FinOptionTypes.AMERICAN_CALL
    euOptionType = FinOptionTypes.EUROPEAN_CALL
    
    amOption = FinEquityAmericanOption(expiry_date, strikePrice,
                                       optionType)
    
    ameuOption = FinEquityAmericanOption(expiry_date, strikePrice, 
                                         euOptionType)
    
    euOption = FinEquityVanillaOption(expiry_date, strikePrice,
                                      euOptionType)
    
    discount_curve = FinDiscountCurveFlat(valuation_date, interestRate,
                                         FinFrequencyTypes.CONTINUOUS, 
                                         FinDayCountTypes.ACT_365F)

    dividendCurve = FinDiscountCurveFlat(valuation_date, dividendYield,
                                         FinFrequencyTypes.CONTINUOUS, 
                                         FinDayCountTypes.ACT_365F)
    
    num_steps_per_year = 400
    
    modelTree = FinModelBlackScholes(volatility, 
                                     FinModelBlackScholesTypes.CRR_TREE, 
                                     num_steps_per_year)
    
    v = amOption.value(valuation_date, stockPrice, discount_curve, 
                           dividendCurve, modelTree)
#    print(v)

    modelApprox = FinModelBlackScholes(volatility, 
                                       FinModelBlackScholesTypes.BARONE_ADESI)

    v = amOption.value(valuation_date, stockPrice, discount_curve, 
                       dividendCurve, modelApprox)

#    print(v)

    v = ameuOption.value(valuation_date, stockPrice, discount_curve, 
                           dividendCurve, modelTree)

#    print(v)

    v = euOption.value(valuation_date, stockPrice, discount_curve, 
                         dividendCurve, modelTree)

#    print(v)

    amTreeValue = []
    amBAWValue = []
    euTreeValue = []
    euAnalValue = []
    volatility = 0.20

    # num_steps_per_year = range(5, 200, 1)
    
    # for numSteps in num_steps_per_year:

    #     modelTree = FinModelBlackScholes(volatility,
    #                                      FinModelBlackScholesTypes.CRR_TREE,
    #                                      {'num_steps_per_year':numSteps})

    #     modelAnal = FinModelBlackScholes(volatility, 
    #                                      FinModelBlackScholesTypes.ANALYTICAL)

    #     modelBAW = FinModelBlackScholes(volatility, 
    #                                     FinModelBlackScholesTypes.BARONE_ADESI)


    #     v_am = amOption.value(valuation_date, stockPrice, discount_curve, 
    #                           dividendYield, modelTree)

    #     v_eu = ameuOption.value(valuation_date, stockPrice, discount_curve, 
    #                             dividendYield, modelTree)
 
    #     v_bs = euOption.value(valuation_date, stockPrice, discount_curve, 
    #                           dividendYield, modelAnal)

    #     v_am_baw = amOption.value(valuation_date, stockPrice, discount_curve, 
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
