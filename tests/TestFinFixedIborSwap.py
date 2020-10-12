###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode

import numpy as np

from financepy.finutils.FinMath import ONE_MILLION
from financepy.products.funding.FinLiborCurve import FinLiborCurve
from financepy.products.funding.FinIborSwap import FinIborSwap
from financepy.products.funding.FinFixedIborSwap import FinFixedIborSwap
from financepy.products.funding.FinIborFRA import FinIborFRA
from financepy.products.funding.FinIborDeposit import FinIborDeposit
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinInterpolate import FinInterpTypes
from financepy.finutils.FinGlobalTypes import FinSwapTypes
from financepy.finutils.FinSchedule import FinSchedule

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def buildLiborCurve(valuationDate):

    settlementDate = valuationDate.addDays(2)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    maturityDate = settlementDate.addMonths(1)
    depo1 = FinIborDeposit(settlementDate, maturityDate, -0.00251, dcType)
    depos.append(depo1)

    # Series of 1M futures
    startDate = settlementDate.nextIMMDate()
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.0023, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00234, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00225, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00226, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00219, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00213, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00186, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00189, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00175, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00143, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00126, dcType)
    fras.append(fra)

    startDate = startDate.addMonths(1)
    endDate = startDate.addMonths(1)
    fra = FinIborFRA(startDate, endDate, -0.00126, dcType)
    fras.append(fra)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    fixedFreq = FinFrequencyTypes.ANNUAL
    floatFreq = FinFrequencyTypes.QUARTERLY
    calType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    dcType = FinDayCountTypes.THIRTY_E_360
    swapType = FinSwapTypes.PAYER

    #######################################
    maturityDate = settlementDate.addMonths(24) 
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                               calType, busDayAdjustType, dateGenRuleType).scheduleDates()   
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = -0.001506    
    swap1 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap1)

    #######################################
    maturityDate = settlementDate.addMonths(36)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()   
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = -0.000185 
    swap2 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap2)

    #######################################
    maturityDate = settlementDate.addMonths(48)   
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                calType, busDayAdjustType, dateGenRuleType).scheduleDates()  
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.001358
    swap3 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap3)

    #######################################
    maturityDate = settlementDate.addMonths(60)   
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()   
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0027652
    swap4 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap4)

    #######################################
    maturityDate = settlementDate.addMonths(72)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0041539
    swap5 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap5)

    #######################################
    maturityDate = settlementDate.addMonths(84)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()    
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0054604
    swap6 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap6)

    #######################################
    maturityDate = settlementDate.addMonths(96)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()    
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.006674
    swap7 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap7)

    #######################################
    maturityDate = settlementDate.addMonths(108)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.007826
    swap8 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap8)

    #######################################
    maturityDate = settlementDate.addMonths(120)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()    
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.008821
    swap9 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap9)

    #######################################
    maturityDate = settlementDate.addMonths(132)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0097379
    swap10 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap10)

    #######################################
    maturityDate = settlementDate.addMonths(144)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0105406
    swap11 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap11)

    #######################################
    maturityDate = settlementDate.addMonths(180)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0123927
    swap12 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap12)

    #######################################
    maturityDate = settlementDate.addMonths(240)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0139882
    swap13 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                             fixedLegDates, dcType, swapRate,
                             floatLegDates, dcType)
    swaps.append(swap13)

    #######################################
    maturityDate = settlementDate.addMonths(300)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()

    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0144972
    swap14 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap14)

    #######################################
    maturityDate = settlementDate.addMonths(360)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0146081
    swap15 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap15)

    #######################################
    maturityDate = settlementDate.addMonths(420)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.01461897
    swap16 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap16)

    #######################################
    maturityDate = settlementDate.addMonths(480)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.014567455
    swap17 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap17)

    #######################################
    maturityDate = settlementDate.addMonths(540)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.0140826
    swap18 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap18)

    #######################################
    maturityDate = settlementDate.addMonths(600)
    floatLegDates = FinSchedule(settlementDate, maturityDate, floatFreq, 
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(settlementDate, maturityDate, fixedFreq,
                                   calType, busDayAdjustType, dateGenRuleType).scheduleDates()
    swapRate = 0.01436822
    swap19 = FinFixedIborSwap(settlementDate, maturityDate, swapType,
                              fixedLegDates, dcType, swapRate,
                              floatLegDates, dcType)
    swaps.append(swap19)
    
    liborCurve = FinLiborCurve(settlementDate, depos, fras, swaps)

    testCases.header("LABEL", "DATE", "VALUE")

    ''' Check calibration '''
    for depo in depos:
        v = depo.value(settlementDate, liborCurve)
        testCases.print("DEPO VALUE:", depo._maturityDate, v)

    for fra in fras:
        v = fra.value(settlementDate, liborCurve)
        testCases.print("FRA VALUE:", fra._maturityDate, v)
    
    for swap in swaps:
        v = swap.value(settlementDate, liborCurve, liborCurve, None)
        testCases.print("SWAP VALUE:", swap._maturityDate, v)

    return liborCurve

###############################################################################


def test_LiborSwap():

    # I have tried to reproduce the example from the blog by Ioannis Rigopoulos
    # https://blog.deriscope.com/index.php/en/excel-interest-rate-swap-price-dual-bootstrapping-curve
    startDate = FinDate(27, 12, 2017)
    endDate = FinDate(27, 12, 2067)

    fixedCoupon = 0.015
    fixedFreqType = FinFrequencyTypes.ANNUAL
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360

    floatSpread = 0.0
    floatFreqType = FinFrequencyTypes.SEMI_ANNUAL
    floatDayCountType = FinDayCountTypes.ACT_360
    firstFixing = -0.00268

    swapCalendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    swapType = FinSwapTypes.RECEIVER
    
    notional = 10.0 * ONE_MILLION

    floatLegDates = FinSchedule(startDate, endDate, floatFreqType, swapCalendarType,
                                busDayAdjustType, dateGenRuleType).scheduleDates()
    fixedLegDates = FinSchedule(startDate, endDate, fixedFreqType, swapCalendarType,
                                busDayAdjustType, dateGenRuleType).scheduleDates()
    
    swap = FinFixedIborSwap(startDate,
                            endDate,
                            swapType,
                            fixedLegDates,
                            fixedDayCountType,
                            fixedCoupon,
                            floatLegDates,
                            floatDayCountType,
                            floatSpread,
                            notional)

    ''' Now perform a valuation after the swap has seasoned but with the
    same curve being used for discounting and working out the implied
    future Libor rates. '''

    valuationDate = FinDate(30, 11, 2018)
    settlementDate = valuationDate.addDays(2)
    liborCurve = buildLiborCurve(valuationDate)
    v = swap.value(settlementDate, liborCurve, liborCurve, firstFixing)

    v_bbg = 388147.0
    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE USING ONE_CURVE", v)
    testCases.print("BLOOMBERG VALUE", v_bbg)
    testCases.print("DIFFERENCE VALUE", v_bbg - v)

###############################################################################

test_LiborSwap()
testCases.compareTestCases()
