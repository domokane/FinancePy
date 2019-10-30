# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""
import sys
sys.path.append("..//..")


from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes

from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve

from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap

from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloor
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorType
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorModelTypes

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinLiborDepositsAndSwaps(valuationDate):
    
    depoBasis = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)
 
    depositRate = 0.020
    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

    maturityDate = settlementDate.addMonths(9)
    depo5 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

    maturityDate = settlementDate.addMonths(12)
    depo6 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoBasis)

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
    swap1 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap1)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap2)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap3)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap4)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap5)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap6)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap7)

    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap8)
        
    swapRate += 0.0025
    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreq, fixedBasis)
    swaps.append(swap9)
    
    liborCurve = FinLiborOneCurve("USD_LIBOR", 
                                  settlementDate, 
                                  depos, 
                                  fras, 
                                  swaps) 
        
    return liborCurve

################################################################################
    
def test_FinLiborCapFloor():       

    import time

    todayDate = FinDate(2019,6,20)
    startDate = FinDate(2019,6,24) # Need to check this how it aligns with curve
    valuationDate = startDate
    maturityDate = startDate.addMonths(60)
    strikeRate = 0.050
    liborCurve = test_FinLiborDepositsAndSwaps(todayDate)

    # The capfloor has begun
    #lastFixing = 0.028
    
    ############################################################################
    ############################ BLACK  ########################################
    ############################################################################

    start = time.time()

    testCases.header("LABEL","VALUE")    
    testCases.banner("==================== BLACK =======================")

    modelType = FinLiborCapFloorModelTypes.BLACK
    modelParams = {'volatility':0.25}

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)
#    capfloor.print()

    testCases.print("CAP Value:",value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)
#    capfloor.print()

    testCases.print("FLOOR Value:",value)
    end = time.time()
    
    testCases.header("TIME")    
    testCases.print(end-start)

    ############################################################################
    #######################  SHIFTED BLACK  ####################################
    ############################################################################

    testCases.banner("============== SHIFTED BLACK ==================")
 
    testCases.header("LABEL","VALUE")    
 
    start = time.time()

    modelType = FinLiborCapFloorModelTypes.SHIFTED_BLACK
    modelParams = {'volatility':0.25, 'shift':-0.01}

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)    
#    capfloor.print()

    testCases.print("CAP Value:",value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)    
#    capfloor.print()

    testCases.print("FLOOR Value:",value)
    end = time.time()

    testCases.header("TIME")    
    testCases.print(end-start)

    ############################################################################
    ################################ SABR  #####################################
    ############################################################################

    testCases.header("LABEL","VALUE")    
 
    start = time.time()

    testCases.banner("============== SABR ==================")

    modelType = FinLiborCapFloorModelTypes.SABR
    modelParams = {'alpha':0.28, 'beta':1.0, 'rho':-0.09, 'nu':0.21}

    capFloorType = FinLiborCapFloorType.CAP
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)    
#    capfloor.print()

    testCases.print("CAP Value:",value)

    capFloorType = FinLiborCapFloorType.FLOOR
    capfloor = FinLiborCapFloor(startDate, maturityDate, capFloorType, strikeRate)
    value = capfloor.value(valuationDate, liborCurve, modelType, modelParams)    
 #   capfloor.print()

    testCases.print("FLOOR Value:",value)
    end = time.time()
    
    testCases.header("TIME")    
    testCases.print(end-start)
    
test_FinLiborCapFloor()
testCases.compareTestCases()
