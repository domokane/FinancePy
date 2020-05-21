# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""
import time

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorType
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloor
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate

from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.market.curves.FinZeroCurve import FinZeroCurve
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloor, FinLiborCapFloorType
from financepy.market.curves.FinInterpolate import interpolate, FinInterpMethods

from financepy.products.libor.FinLiborModelTypes import FinLiborModelBlack
from financepy.products.libor.FinLiborModelTypes import FinLiborModelShiftedBlack
from financepy.products.libor.FinLiborModelTypes import FinLiborModelSABR

import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinLiborDepositsAndSwaps(valuationDate):

    depoBasis = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.020
    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    maturityDate = settlementDate.addMonths(9)
    depo5 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    maturityDate = settlementDate.addMonths(12)
    depo6 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoBasis)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)
    depos.append(depo6)

    fras = []

    swaps = []
    fixedBasis = FinDayCountTypes.ACT_365_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.02
    maturityDate = settlementDate.addMonths(24)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap1)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap2)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap3)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap4)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap5)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap6)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap7)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap8)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap9)

    liborCurve = FinLiborCurve("USD_LIBOR",
                                  settlementDate,
                                  depos,
                                  fras,
                                  swaps)

    return liborCurve

##########################################################################


def test_FinLiborCapFloor():

    import time

    todayDate = FinDate(2019, 6, 20)
    # Need to check this how it aligns with curve
    startDate = FinDate(2019, 6, 24)
    valuationDate = startDate
    maturityDate = startDate.addMonths(60)
    strikeRate = 0.050
    liborCurve = test_FinLiborDepositsAndSwaps(todayDate)

    # The capfloor has begun
    #lastFixing = 0.028

    ##########################################################################
    ############################ BLACK  ######################################
    ##########################################################################

    start = time.time()

    testCases.header("LABEL", "VALUE")
    testCases.banner("==================== BLACK =======================")

    model = FinLiborModelBlack(0.25)

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    
    print(capfloor)
    
    value = capfloor.value(valuationDate, liborCurve, model)
#    capfloor.print()

    testCases.print("CAP Value:", value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    value = capfloor.value(valuationDate, liborCurve, model)
#    capfloor.print()

    testCases.print("FLOOR Value:", value)
    end = time.time()

    testCases.header("TIME")
    testCases.print(end - start)

    ##########################################################################
    #######################  SHIFTED BLACK  ##################################
    ##########################################################################

    testCases.banner("============== SHIFTED BLACK ==================")

    testCases.header("LABEL", "VALUE")

    start = time.time()

    model = FinLiborModelShiftedBlack(0.25, -0.01)

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    value = capfloor.value(valuationDate, liborCurve, model)
#    capfloor.print()

    testCases.print("CAP Value:", value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    value = capfloor.value(valuationDate, liborCurve, model)
#    capfloor.print()

    testCases.print("FLOOR Value:", value)
    end = time.time()

    testCases.header("TIME")
    testCases.print(end - start)

    ##########################################################################
    ################################ SABR  ###################################
    ##########################################################################

    testCases.header("LABEL", "VALUE")

    start = time.time()

    testCases.banner("============== SABR ==================")

    model = FinLiborModelSABR(0.28, 1.0, -0.09, 0.21)

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    value = capfloor.value(valuationDate, liborCurve, model)
#    capfloor.print()

    testCases.print("CAP Value:", value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(
        startDate,
        maturityDate,
        capFloorType,
        strikeRate)
    value = capfloor.value(valuationDate, liborCurve, model)

    ###########################################################################

    if 1 == 1: # TESTING PRINT LEGS

        model = FinLiborModelBlack(0.25)

        maturityDate = startDate.addMonths(12)
        strikeRate = 0.025

        capFloorType = FinLiborCapFloorType.CAP
        capfloor = FinLiborCapFloor(
            startDate,
            maturityDate,
            capFloorType,
            strikeRate)

        value = capfloor.value(valuationDate, liborCurve, model)
        capfloor.print()
        capfloor.printLeg()

    if 1 == 0: # TESTING PRINT LEGS

        capFloorType = FinLiborCapFloorType.FLOOR
        capfloor = FinLiborCapFloor(
            startDate,
            maturityDate,
            capFloorType,
            strikeRate)

        value = capfloor.value(valuationDate, liborCurve, model)

        print(value)
        capfloor.print()
        capfloor.printLeg()

###############################################################################




def test_FinLiborCapFloorQLExample():

    valuationDate = FinDate(14, 6, 2016)

    dates = [FinDate(14, 6, 2016), FinDate(14, 9, 2016),
             FinDate(14, 12, 2016), FinDate(14, 6, 2017),
             FinDate(14, 6, 2019), FinDate(14, 6, 2021),
             FinDate(15, 6, 2026), FinDate(16, 6, 2031),
             FinDate(16, 6, 2036), FinDate(14, 6, 2046)]

    rates = [0.000000, 0.006616, 0.007049, 0.007795,
             0.009599, 0.011203, 0.015068, 0.017583,
             0.018998, 0.020080]

    frequencyType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA

    discountCurve = FinZeroCurve(valuationDate, dates, rates, frequencyType,
                                 dayCountType,
                                 FinInterpMethods.LINEAR_ZERO_RATES)

    startDate = FinDate(14, 6, 2016)
    endDate = FinDate(14, 6, 2026)
    calendarType = FinCalendarTypes.US
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    frequencyType = FinFrequencyTypes.QUARTERLY
    dateGenRuleType = FinDateGenRuleTypes.FORWARD
    lastFixing = 0.0065560
    notional = 1000000
    dayCountType = FinDayCountTypes.ACT_360
    optionType = FinLiborCapFloorType.CAP
    strikeRate = 0.02

    cap = FinLiborCapFloor(startDate, endDate, optionType, strikeRate,
                           lastFixing, frequencyType,  dayCountType, notional,
                           calendarType, busDayAdjustType, dateGenRuleType)

    blackVol = 0.547295
    model = FinLiborModelBlack(blackVol)

    start = time.time()
    numRepeats = 10
    for i in range(0, numRepeats):
        v = cap.value(valuationDate, discountCurve, model)

    end = time.time()
    period = end - start
#    print(v, period/numRepeats)


test_FinLiborCapFloor()
test_FinLiborCapFloorQLExample()
# testCases.compareTestCases()
