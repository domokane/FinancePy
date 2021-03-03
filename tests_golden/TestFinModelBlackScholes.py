###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes

from financepy.finutils.FinGlobalTypes import FinOptionTypes

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

    valueDate = FinDate(8, 5, 2015)
    expiryDate = FinDate(15, 1, 2016)

    strikePrice = 130.0
    stockPrice = 127.62
    volatility = 0.20
    interestRate = 0.001
    dividendYield = 0.0163

    optionType = FinOptionTypes.AMERICAN_CALL
    euOptionType = FinOptionTypes.EUROPEAN_CALL
    
    amOption = FinEquityAmericanOption(expiryDate, strikePrice,
                                       optionType)
    
    ameuOption = FinEquityAmericanOption(expiryDate, strikePrice, 
                                         euOptionType)
    
    euOption = FinEquityVanillaOption(expiryDate, strikePrice,
                                      euOptionType)
    
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate,
                                         FinFrequencyTypes.CONTINUOUS, 
                                         FinDayCountTypes.ACT_365F)

    dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield,
                                         FinFrequencyTypes.CONTINUOUS, 
                                         FinDayCountTypes.ACT_365F)
    
    numStepsPerYear = 400
    
    modelTree = FinModelBlackScholes(volatility, 
                                     FinModelBlackScholesTypes.CRR_TREE, 
                                     numStepsPerYear)
    
    v = amOption.value(valueDate, stockPrice, discountCurve, 
                           dividendCurve, modelTree)
#    print(v)

    modelApprox = FinModelBlackScholes(volatility, 
                                       FinModelBlackScholesTypes.BARONE_ADESI)

    v = amOption.value(valueDate, stockPrice, discountCurve, 
                       dividendCurve, modelApprox)

#    print(v)

    v = ameuOption.value(valueDate, stockPrice, discountCurve, 
                           dividendCurve, modelTree)

#    print(v)

    v = euOption.value(valueDate, stockPrice, discountCurve, 
                         dividendCurve, modelTree)

#    print(v)

    amTreeValue = []
    amBAWValue = []
    euTreeValue = []
    euAnalValue = []
    volatility = 0.20

    # numStepsPerYear = range(5, 200, 1)
    
    # for numSteps in numStepsPerYear:

    #     modelTree = FinModelBlackScholes(volatility,
    #                                      FinModelBlackScholesTypes.CRR_TREE,
    #                                      {'numStepsPerYear':numSteps})

    #     modelAnal = FinModelBlackScholes(volatility, 
    #                                      FinModelBlackScholesTypes.ANALYTICAL)

    #     modelBAW = FinModelBlackScholes(volatility, 
    #                                     FinModelBlackScholesTypes.BARONE_ADESI)


    #     v_am = amOption.value(valueDate, stockPrice, discountCurve, 
    #                           dividendYield, modelTree)

    #     v_eu = ameuOption.value(valueDate, stockPrice, discountCurve, 
    #                             dividendYield, modelTree)
 
    #     v_bs = euOption.value(valueDate, stockPrice, discountCurve, 
    #                           dividendYield, modelAnal)

    #     v_am_baw = amOption.value(valueDate, stockPrice, discountCurve, 
    #                               dividendYield, modelBAW)
        
    #     amTreeValue.append(v_am)
    #     euTreeValue.append(v_eu)
    #     euAnalValue.append(v_bs)
    #     amBAWValue.append(v_am_baw)
        
    
    # plt.title("American PUT Option Price Convergence Analysis")
    # plt.plot(numStepsPerYear, amTreeValue, label="American Tree")
    # plt.plot(numStepsPerYear, amBAWValue, label="American BAW")
    # plt.plot(numStepsPerYear, euTreeValue, label="European Tree")
    # plt.plot(numStepsPerYear, euAnalValue, label="European Anal", lw =2)
    # plt.xlabel("Num Steps")
    # plt.ylabel("Value")
    # plt.legend();

###############################################################################


testFinModelBlackScholes()
testCases.compareTestCases()
