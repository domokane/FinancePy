# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""

import sys
sys.path.append("..//..")

import numpy as np

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve

from financepy.products.libor.FinLiborFRA import FinLiborFRA
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap

from financepy.finutils.FinGlobalVariables import gDaysInYear

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinLiborDepositsOnly():
    
    # I have used the following useful blog post by Ioannis Rigopoulos for this
    # https://blog.deriscope.com/index.php/en/yield-curve-excel-quantlib-deposit

    valuationDate = FinDate(2018,2,23)

    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depoDCCType = FinDayCountTypes.ACT_360
    notional = 100.0
    calendarType = FinCalendarTypes.TARGET
    depos = []

    # 1 month
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(1)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, 
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 2 months
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(2)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, 
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 6 months
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(6)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, 
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    # 1 year
    depositRate = 0.04
    maturityDate = settlementDate.addMonths(12)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, 
                           depoDCCType, notional, calendarType)
    depos.append(depo)

    fras = []
    swaps = []

    liborCurve = FinLiborOneCurve("USD_LIBOR", 
                                  settlementDate, 
                                  depos, 
                                  fras, 
                                  swaps) 

    liborCurve.print()

    testCases.header("DATE","DF")

    ''' Check calibration '''
    for depo in depos:
        v = depo.value(settlementDate, liborCurve )
        print("DEPO:",depo._maturityDate,v) 

    times = []
    dfs = []

    dt = settlementDate
    for i in range(0,11):
        df = liborCurve.df(dt)
        zero = liborCurve.zeroContinuous(dt)
        fwd = liborCurve.fwdContinuous(dt)
        t = (dt - settlementDate) / gDaysInYear
        print(dt,"%10.7f"%t,"%10.7f"%df,"%10.7f"%zero,"%10.7f"%fwd)
        dt = dt.addMonths(1)

###############################################################################

def test_FinLiborFRAsOnly():
    
    valuationDate = FinDate(2018,2,23)

    spotDays = 0
    settlementDate = valuationDate.addWorkDays(spotDays)

    depoDCCType = FinDayCountTypes.ACT_360
    notional = 100.0

    payFixed = True

    calendarType = FinCalendarTypes.TARGET
    fras = []

    # 1 x 4 FRA
    fraRate = 0.04
    fraSettlementDate = settlementDate.addMonths(1)
    fraMaturityDate = settlementDate.addMonths(4)
    fra = FinLiborFRA(settlementDate, fraMaturityDate, fraRate, payFixed, 
                           depoDCCType, notional, calendarType)
    fras.append(fra)

    # 4 x 7 FRA
    fraRate = 0.04
    fraSettlementDate = settlementDate.addMonths(4)
    fraMaturityDate = settlementDate.addMonths(7)
    fra = FinLiborFRA(settlementDate, fraMaturityDate, fraRate, payFixed,
                           depoDCCType, notional, calendarType)
    fras.append(fra)

    depos = []
    swaps = []

    liborCurve = FinLiborOneCurve("USD_LIBOR", 
                                  settlementDate, 
                                  depos, 
                                  fras, 
                                  swaps) 

    liborCurve.print()

    testCases.header("DATE","DF")

    ''' Check calibration '''
    for fra in fras:
        v = fra.value(settlementDate, liborCurve)
        print("FRA:",fra._maturityDate,v) 

    dt = settlementDate
    for i in range(0,11):
        df = liborCurve.df(dt)
        zero = liborCurve.zeroContinuous(dt)
        fwd = liborCurve.fwdContinuous(dt)
        t = (dt - settlementDate) / gDaysInYear
        print(dt,"%10.7f"%t,"%10.7f"%df,"%10.7f"%zero,"%10.7f"%fwd)
        dt = dt.addMonths(1)

################################################################################

def test_FinLiborDepositsAndSwaps():

    valuationDate = FinDate(2019,9,18)

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)
 
    depositRate = 0.050
    maturityDate = settlementDate.addMonths(1)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(2)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(3)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(6)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(9)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    maturityDate = settlementDate.addMonths(12)
    depo = FinLiborDeposit(settlementDate, maturityDate, depositRate, depoDCCType)
    depos.append(depo)

    fras = []
    
    swaps = []
    fixedDCCType = FinDayCountTypes.ACT_365_ISDA
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.05        
    maturityDate = settlementDate.addMonths(24)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(36)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(48)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(60)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(72)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(84)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(96)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(108)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)
        
    maturityDate = settlementDate.addMonths(120)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(132)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(144)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(180)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(240)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(300)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)

    maturityDate = settlementDate.addMonths(360)
    swap = FinLiborSwap(settlementDate, maturityDate, swapRate, fixedFreqType, fixedDCCType)
    swaps.append(swap)
    
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

    for swap in swaps:
        df = liborCurve.df(deposit._maturityDate)
        testCases.print(str(deposit._maturityDate),df)

################################################################################

test_FinLiborDepositsOnly()
test_FinLiborFRAsOnly()

test_FinLiborDepositsAndSwaps()

testCases.compareTestCases()
