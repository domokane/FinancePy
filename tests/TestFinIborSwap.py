###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
import numpy as np
sys.path.append("..")

from financepy.utils.Math import ONE_MILLION
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborFRA import FinIborFRA
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.utils.Calendar import FinBusDayAdjustTypes
from financepy.utils.Calendar import FinDateGenRuleTypes
from financepy.utils.Calendar import FinCalendarTypes
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.Date import Date
from financepy.utils.FinGlobalTypes import FinSwapTypes
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinInterpolator import FinInterpTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def buildIborSingleCurve(valuation_date):

    settlement_date = valuation_date.addDays(2)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    maturity_date = settlement_date.addMonths(1)
    depo1 = FinIborDeposit(valuation_date, maturity_date, -0.00251, dcType)
    depos.append(depo1)

    # Series of 1M futures
    start_date = settlement_date.nextIMMDate()
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.0023, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00234, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00225, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00226, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00219, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00213, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00186, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00189, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00175, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00143, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00126, dcType)
    fras.append(fra)

    start_date = start_date.addMonths(1)
    end_date = start_date.addMonths(1)
    fra = FinIborFRA(start_date, end_date, -0.00126, dcType)
    fras.append(fra)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    fixedFreq = FinFrequencyTypes.ANNUAL
    dcType = FinDayCountTypes.THIRTY_E_360
    fixed_legType = FinSwapTypes.PAY

    #######################################
    maturity_date = settlement_date.addMonths(24) 
    swapRate = -0.001506    
    swap1 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap1)

    #######################################
    maturity_date = settlement_date.addMonths(36)
    swapRate = -0.000185 
    swap2 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap2)

    #######################################
    maturity_date = settlement_date.addMonths(48)   
    swapRate = 0.001358
    swap3 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap3)

    #######################################
    maturity_date = settlement_date.addMonths(60)   
    swapRate = 0.0027652
    swap4 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap4)

    #######################################
    maturity_date = settlement_date.addMonths(72)
    swapRate = 0.0041539
    swap5 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap5)

    #######################################
    maturity_date = settlement_date.addMonths(84)
    swapRate = 0.0054604
    swap6 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap6)

    #######################################
    maturity_date = settlement_date.addMonths(96)
    swapRate = 0.006674
    swap7 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap7)

    #######################################
    maturity_date = settlement_date.addMonths(108)
    swapRate = 0.007826
    swap8 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap8)

    #######################################
    maturity_date = settlement_date.addMonths(120)
    swapRate = 0.008821
    swap9 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap9)

    #######################################
    maturity_date = settlement_date.addMonths(132)
    swapRate = 0.0097379
    swap10 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap10)

    #######################################
    maturity_date = settlement_date.addMonths(144)
    swapRate = 0.0105406
    swap11 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap11)

    #######################################
    maturity_date = settlement_date.addMonths(180)
    swapRate = 0.0123927
    swap12 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap12)

    #######################################
    maturity_date = settlement_date.addMonths(240)
    swapRate = 0.0139882
    swap13 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap13)

    #######################################
    maturity_date = settlement_date.addMonths(300)
    swapRate = 0.0144972
    swap14 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap14)

    #######################################
    maturity_date = settlement_date.addMonths(360)
    swapRate = 0.0146081
    swap15 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap15)

    #######################################
    maturity_date = settlement_date.addMonths(420)
    swapRate = 0.01461897
    swap16 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap16)

    #######################################
    maturity_date = settlement_date.addMonths(480)
    swapRate = 0.014567455
    swap17 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap17)

    #######################################
    maturity_date = settlement_date.addMonths(540)
    swapRate = 0.0140826
    swap18 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap18)

    #######################################
    maturity_date = settlement_date.addMonths(600)
    swapRate = 0.01436822
    swap19 = FinIborSwap(settlement_date, maturity_date, fixed_legType,
                             swapRate, fixedFreq, dcType)
    swaps.append(swap19)
    
    ########################################
    
    libor_curve = FinIborSingleCurve(valuation_date, depos, fras, swaps)

    testCases.header("LABEL", "DATE", "VALUE")

    """ Check calibration """
    for depo in depos:
        v = depo.value(settlement_date, libor_curve)
        testCases.print("DEPO VALUE:", depo._maturity_date, v)

    for fra in fras:
        v = fra.value(settlement_date, libor_curve)
        testCases.print("FRA VALUE:", fra._maturity_date, v)
    
    for swap in swaps:
        v = swap.value(settlement_date, libor_curve)
        testCases.print("SWAP VALUE:", swap._maturity_date, v)

    return libor_curve

###############################################################################


def test_LiborSwap():

    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve
    start_date = Date(27, 12, 2017)
    end_date = Date(27, 12, 2067)

    fixedCoupon = 0.015
    fixedFreqType = FinFrequencyTypes.ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360

    floatSpread = 0.0
    floatFreqType = FinFrequencyTypes.SEMI_ANNUAL
    floatDayCountType = FinDayCountTypes.ACT_360
    firstFixing = -0.00268

    swapCalendarType = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    fixed_legType = FinSwapTypes.RECEIVE
    
    notional = 10.0 * ONE_MILLION

    swap = FinIborSwap(start_date,
                            end_date,
                            fixed_legType,
                            fixedCoupon,
                            fixedFreqType,
                            fixedDayCountType,
                            notional,
                            floatSpread,
                            floatFreqType,
                            floatDayCountType,
                            swapCalendarType,
                            bus_day_adjust_type,
                            date_gen_rule_type)

    """ Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. """

    valuation_date = Date(30, 11, 2018)
    settlement_date = valuation_date.addDays(2)
    libor_curve = buildIborSingleCurve(valuation_date)
    v = swap.value(settlement_date, libor_curve, libor_curve, firstFixing)

    v_bbg = 388147.0
    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE USING ONE_CURVE", v)
    testCases.print("BLOOMBERG VALUE", v_bbg)
    testCases.print("DIFFERENCE VALUE", v_bbg - v)

###############################################################################


def test_dp_example():

    #  http://www.derivativepricing.com/blogpage.asp?id=8

    start_date = Date(14, 11, 2011)
    end_date = Date(14, 11, 2016)
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL
    swapCalendarType = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixed_legType = FinSwapTypes.PAY
    fixedCoupon = 0.0124
    notional = ONE_MILLION

    swap = FinIborSwap(start_date,
                        end_date,
                        fixed_legType,
                        fixedCoupon=fixedCoupon,
                        fixedFreqType=fixedFreqType,
                        fixedDayCountType=fixedDayCountType,
                        floatFreqType=FinFrequencyTypes.SEMI_ANNUAL,
                        floatDayCountType=FinDayCountTypes.ACT_360,
                        notional=notional,
                        calendar_type=swapCalendarType,
                        bus_day_adjust_type=bus_day_adjust_type,
                        date_gen_rule_type=date_gen_rule_type)

#    swap.printFixedLegFlows()

    dts = [Date(14, 11, 2011), Date(14, 5, 2012), Date(14, 11, 2012),
           Date(14, 5, 2013), Date(14, 11, 2013), Date(14, 5, 2014),
           Date(14, 11, 2014), Date(14, 5, 2015), Date(16, 11, 2015),
           Date(16, 5, 2016), Date(14, 11, 2016)]

    dfs = [0.9999843, 0.9966889, 0.9942107, 0.9911884, 0.9880738, 0.9836490,
           0.9786276, 0.9710461, 0.9621778, 0.9514315, 0.9394919]

    valuation_date = start_date

    curve = FinDiscountCurve(valuation_date, dts, np.array(dfs),
                             FinInterpTypes.FLAT_FWD_RATES)

    v = swap.value(valuation_date, curve, curve)

#    swap.printFixedLegPV()
#    swap.printFloatLegPV()

    # This is essentially zero
    testCases.header("LABEL", "VALUE")
    testCases.print("Swap Value on a Notional of $1M:", v)

###############################################################################

test_LiborSwap()
test_dp_example()
testCases.compareTestCases()
