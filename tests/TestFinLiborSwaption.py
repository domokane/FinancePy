# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.libor.FinLiborSwaption import FinLiborSwaptionType
from financepy.products.libor.FinLiborSwaption import FinLiborSwaption
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate

from financepy.products.libor.FinLiborModelTypes import FinLiborModelBlack
from financepy.products.libor.FinLiborModelTypes import FinLiborModelSABR


import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinLiborDepositsAndSwaps(valuationDate):

    depoBasis = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.030
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

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    swaps = []
    fixedBasis = FinDayCountTypes.ACT_365_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL

    swapRate = 0.03
    maturityDate = settlementDate.addMonths(12)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap1)

    swapRate = 0.034
    maturityDate = settlementDate.addMonths(24)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap2)

    swapRate = 0.037
    maturityDate = settlementDate.addMonths(36)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap3)

    swapRate = 0.039
    maturityDate = settlementDate.addMonths(48)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap4)

    swapRate = 0.040
    maturityDate = settlementDate.addMonths(60)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreq,
        fixedBasis)
    swaps.append(swap5)

    liborCurve = FinLiborCurve("USD_LIBOR",
                               settlementDate,
                               depos,
                               fras,
                               swaps)

    return liborCurve

##########################################################################


def test_FinLiborSwaption():

    import time

    valuationDate = FinDate(2011, 1, 1)
    exerciseDate = FinDate(2012, 1, 1)
    swapMaturityDate = FinDate(2017, 1, 1)

    swapFixedCoupon = 0.030
    swapFixedFrequencyType = FinFrequencyTypes.SEMI_ANNUAL
    swapFixedDayCountType = FinDayCountTypes.ACT_365_ISDA

    liborCurve = test_FinLiborDepositsAndSwaps(valuationDate)

    start = time.time()

    swaptionType = FinLiborSwaptionType.PAYER
    swaption = FinLiborSwaption(exerciseDate,
                                swapMaturityDate,
                                swaptionType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType)

    model = FinLiborModelBlack(0.25)
    settlementDate = valuationDate.addWorkDays(2)
    value = swaption.value(settlementDate, liborCurve, model)

#    swaption.print()

    testCases.header("LABEL", "VALUE")
    testCases.print("PAYER Swaption Value", value)

    swaptionType = FinLiborSwaptionType.RECEIVER
    swaption = FinLiborSwaption(exerciseDate,
                                swapMaturityDate,
                                swaptionType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType)

    model = FinLiborModelBlack(0.25)
    settlementDate = valuationDate.addWorkDays(2)
    value = swaption.value(settlementDate, liborCurve, model)

#    swaption.print()

    testCases.print("RECEIVER Swaption Value", value)

    end = time.time()

    testCases.header("TIME")
    testCases.print(end - start)

##########################################################################
##########################################################################

    start = time.time()

    swaptionType = FinLiborSwaptionType.PAYER

    testCases.banner(
        "===================================================================")

    swaption = FinLiborSwaption(exerciseDate,
                                swapMaturityDate,
                                swaptionType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType)

    model = FinLiborModelSABR(0.28, 1.0, -0.09, 0.21)
    value = swaption.value(settlementDate, liborCurve, model)

#    swaption.print()

    testCases.header("LABEL", "VALUE")
    testCases.print("PAYER Swaption Value", value)

    swaptionType = FinLiborSwaptionType.RECEIVER

    swaption = FinLiborSwaption(exerciseDate,
                                swapMaturityDate,
                                swaptionType,
                                swapFixedCoupon,
                                swapFixedFrequencyType,
                                swapFixedDayCountType)

    model = FinLiborModelSABR(0.28, 1.0, -0.09, 0.21)
    value = swaption.value(settlementDate, liborCurve, model)

#    swaption.print()

    testCases.print("RECEIVER Swaption Value", value)

    end = time.time()

    testCases.header("TIME")
    testCases.print(end - start)

##########################################################################


test_FinLiborSwaption()
testCases.compareTestCases()
