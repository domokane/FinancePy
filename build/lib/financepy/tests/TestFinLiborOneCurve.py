# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""

import numpy as np

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinLiborDepositsOnly(depoDCCType):
    
    testCases.header("DEPO CONVENTION")
    testCases.print(str(depoDCCType))

    valuationDate = FinDate(2019,9,18)

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.050

    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(9)
    depo5 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(12)
    depo6 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    depos = []
    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)
    depos.append(depo6)

    fras = []
    swaps = []

    liborCurve = FinLiborOneCurve("USD_LIBOR", 
                                  settlementDate, 
                                  depos, 
                                  fras, 
                                  swaps) 

    testCases.header("DATE","DF")

    for deposit in depos:
        df = liborCurve.df(deposit._maturityDate)
        testCases.print(str(deposit._maturityDate),df)

    numSteps = 40
    dt = 10/numSteps
    times = np.linspace(0.0,10.0,numSteps+1)

    testCases.header("DATE","DF","FWD")

    df0 = 1.0
    for t in times[1:]:
        df1 = liborCurve.df(t)
        fwd = (df0/df1-1.0)/dt
        testCases.print(t,df1,fwd)
        df0 = df1
   
################################################################################

def test_FinLiborDepositsAndSwaps():

    valuationDate = FinDate(2019,9,18)

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)
 
    depositRate = 0.050
    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(9)
    depo5 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    maturityDate = settlementDate.addMonths(12)
    depo6 = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)
    depos.append(depo6)

    fras = []
    
    swaps = []
    fixedDCCType = FinDayCountTypes.ACT_365_ISDA
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.05        
    maturityDate = settlementDate.addMonths(24)
    swap1 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap8)
        
    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap9)

    maturityDate = settlementDate.addMonths(132)
    swap10 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap10)

    maturityDate = settlementDate.addMonths(144)
    swap11 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap11)

    maturityDate = settlementDate.addMonths(180)
    swap12 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap12)

    maturityDate = settlementDate.addMonths(240)
    swap13 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap13)

    maturityDate = settlementDate.addMonths(300)
    swap14 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap14)

    maturityDate = settlementDate.addMonths(360)
    swap15 = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap15)
    
    liborCurve = FinLiborOneCurve("USD_LIBOR", 
                                  settlementDate, 
                                  depos, 
                                  fras, 
                                  swaps) 
                                  
    df = liborCurve.df(settlementDate)

    testCases.header("SETTLEMENT DATE","DF")
    testCases.print(str(settlementDate),df)

    testCases.header("DATE","DF")

    for deposit in depos:
        df = liborCurve.df(deposit._maturityDate)
        testCases.print(str(deposit._maturityDate),df)

################################################################################

depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
test_FinLiborDepositsOnly(depoDCCType)

depoDCCType = FinDayCountTypes.ACT_365_ISDA
test_FinLiborDepositsOnly(depoDCCType)

depoDCCType = FinDayCountTypes.ACT_ACT_ISDA
test_FinLiborDepositsOnly(depoDCCType)

################################################################################

test_FinLiborDepositsAndSwaps()
testCases.compareTestCases()
