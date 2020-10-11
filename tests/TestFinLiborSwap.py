###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode

import numpy as np

from financepy.finutils.FinMath import ONE_MILLION
from financepy.products.funding.FinLiborCurve import FinLiborCurve
from financepy.products.funding.FinLiborSwap import FinLiborSwap
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
from financepy.products.funding.FinLiborSwaption import FinSwapTypes

import sys
sys.path.append("..//..")


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

    fixedFreq = FinFrequencyTypes.ANNUAL
    dcType = FinDayCountTypes.THIRTY_E_360

    swapType = FinSwapTypes.PAYER
    
    maturityDate = settlementDate.addMonths(24)
    swap1 = FinLiborSwap(settlementDate, maturityDate, swapType, -
                         0.001506, fixedFreq, dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(settlementDate,maturityDate,  swapType, -
                         0.000185, fixedFreq, dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.001358,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0027652,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0041539,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0054604,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.006674,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.007826,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.008821,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    maturityDate = settlementDate.addMonths(132)
    swap10 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0097379,
        fixedFreq,
        dcType)
    swaps.append(swap10)

    maturityDate = settlementDate.addMonths(144)
    swap11 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0105406,
        fixedFreq,
        dcType)
    swaps.append(swap11)

    maturityDate = settlementDate.addMonths(180)
    swap12 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0123927,
        fixedFreq,
        dcType)
    swaps.append(swap12)

    maturityDate = settlementDate.addMonths(240)
    swap13 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0139882,
        fixedFreq,
        dcType)
    swaps.append(swap13)

    maturityDate = settlementDate.addMonths(300)
    swap14 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0144972,
        fixedFreq,
        dcType)
    swaps.append(swap14)

    maturityDate = settlementDate.addMonths(360)
    swap15 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0146081,
        fixedFreq,
        dcType)
    swaps.append(swap15)

    maturityDate = settlementDate.addMonths(420)
    swap16 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.01461897,
        fixedFreq,
        dcType)
    swaps.append(swap16)

    maturityDate = settlementDate.addMonths(480)
    swap17 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.014567455,
        fixedFreq,
        dcType)
    swaps.append(swap17)

    maturityDate = settlementDate.addMonths(540)
    swap18 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.0140826,
        fixedFreq,
        dcType)
    swaps.append(swap18)

    maturityDate = settlementDate.addMonths(600)
    swap19 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapType,
        0.01436822,
        fixedFreq,
        dcType)
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
    startDate = FinDate(2017, 12, 27)
    endDate = FinDate(2067, 12, 27)

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

    swap = FinLiborSwap(startDate,
                        endDate,
                        swapType,
                        fixedCoupon,
                        fixedFreqType,
                        fixedDayCountType,
                        notional,
                        floatSpread,
                        floatFreqType,
                        floatDayCountType,
                        swapCalendarType,
                        busDayAdjustType,
                        dateGenRuleType)

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


def test_dp_example():

    #  http://www.derivativepricing.com/blogpage.asp?id=8

    startDate = FinDate(14, 11, 2011)
    endDate = FinDate(14, 11, 2016)
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL
    swapCalendarType = FinCalendarTypes.TARGET
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    fixedDayCountType = FinDayCountTypes.THIRTY_E_360_ISDA
    swapType = FinSwapTypes.PAYER
    fixedCoupon = 0.0124
    notional = ONE_MILLION

    swap = FinLiborSwap(startDate,
                        endDate,
                        swapType,
                        fixedCoupon=fixedCoupon,
                        fixedFreqType=fixedFreqType,
                        fixedDayCountType=fixedDayCountType,
                        floatFreqType=FinFrequencyTypes.SEMI_ANNUAL,
                        floatDayCountType=FinDayCountTypes.ACT_360,
                        notional=notional,
                        calendarType=swapCalendarType,
                        busDayAdjustType=busDayAdjustType,
                        dateGenRuleType=dateGenRuleType)

#    swap.printFixedLegFlows()

    dts = [FinDate(14, 11, 2011), FinDate(14, 5, 2012), FinDate(14, 11, 2012),
           FinDate(14, 5, 2013), FinDate(14, 11, 2013), FinDate(14, 5, 2014),
           FinDate(14, 11, 2014), FinDate(14, 5, 2015), FinDate(16, 11, 2015),
           FinDate(16, 5, 2016), FinDate(14, 11, 2016)]

    dfs = [0.9999843, 0.9966889, 0.9942107, 0.9911884, 0.9880738, 0.9836490,
           0.9786276, 0.9710461, 0.9621778, 0.9514315, 0.9394919]

    valuationDate = startDate

    curve = FinDiscountCurve(valuationDate, dts, np.array(dfs),
                             FinInterpTypes.FLAT_FORWARDS)

    v = swap.value(valuationDate, curve, curve)

#    swap.printFixedLegPV()
#    swap.printFloatLegPV()

    # This is essentially zero
    testCases.header("LABEL", "VALUE")
    testCases.print("Swap Value on a Notional of $1M:", v)

###############################################################################


test_dp_example()
test_LiborSwap()
testCases.compareTestCases()
