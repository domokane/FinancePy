###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.bonds.FinBondFRN import FinBondFRN
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import *
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.finutils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def buildIborCurve(valuationDate):

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    payFixed = FinSwapTypes.PAY

    spotDays = 2
    settlementDate = valuationDate.addWeekDays(spotDays)

    depositRate = 0.050
    maturityDate = settlementDate.addMonths(1)
    depo1 = FinIborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(3)
    depo2 = FinIborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(6)
    depo3 = FinIborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(9)
    depo4 = FinIborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(12)
    depo5 = FinIborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []
    fixedDCCType = FinDayCountTypes.ACT_365F
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL

    swaps = []

    swapRate = 0.05
    maturityDate = settlementDate.addMonths(24)
    swap1 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinIborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    liborCurve = FinIborSingleCurve(valuationDate,
                                    depos,
                                    fras,
                                    swaps)
    
    if 1 == 0:
        import numpy as np
        numSteps = 40
        dt = 10 / numSteps
        times = np.linspace(0.0, 10.0, numSteps + 1)

        df0 = 1.0
        for t in times[1:]:
            df1 = liborCurve.df(t)
            fwd = (df0 / df1 - 1.0) / dt
            print(t, df1, fwd)
            df0 = df1

    return liborCurve

##########################################################################


def test_FinBondFRN():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # I have a day out problem on the accrued interest - should be 71 and not 72 days
    # Other than that agreement on the DM is very good.

    ##########################################################################
    # CITIGROUP FRN SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE")
    issueDate = FinDate(10, 11, 2010)
    maturityDate = FinDate(10, 11, 2021)
    quotedMargin = 0.0025
    freqType = FinFrequencyTypes.QUARTERLY
    accrualType = FinDayCountTypes.THIRTY_E_360
    face = 1000000

    bond = FinBondFRN(issueDate,
                      maturityDate,
                      quotedMargin,
                      freqType,
                      accrualType,
                      face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 96.793
    resetIbor = 0.0143456 - quotedMargin
    currentIbor = 0.0120534
    futureIbors = 0.0130522

    settlementDate = FinDate(21, 7, 2017)

    dm = bond.discountMargin(settlementDate,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    fullPrice = bond.fullPriceFromDM(settlementDate,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accruedDays
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accruedInterest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlementDate,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollarDuration(settlementDate,
                                   resetIbor,
                                   currentIbor,
                                   futureIbors,
                                   dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modifiedDuration = bond.modifiedRateDuration(settlementDate,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyRateDuration(settlementDate,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDM(settlementDate,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Convexity = ", convexity)

    duration = bond.dollarCreditDuration(settlementDate,
                                         resetIbor,
                                         currentIbor,
                                         futureIbors,
                                         dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modifiedDuration = bond.modifiedCreditDuration(settlementDate,
                                                   resetIbor,
                                                   currentIbor,
                                                   futureIbors,
                                                   dm)

    testCases.print("Modified Credit Duration = ", modifiedDuration)

##########################################################################
# EXAMPLE
# https://ebrary.net/14293/economics/actual_floater
##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE II")
    issueDate = FinDate(28, 3, 2000)
    settlementDate = FinDate(28, 3, 2014)
    maturityDate = FinDate(3, 2, 2021)
    quotedMargin = 0.0020
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.THIRTY_E_360_ISDA
    face = 1000000.0

    bond = FinBondFRN(issueDate,
                      maturityDate,
                      quotedMargin,
                      freqType,
                      accrualType,
                      face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 93.08
    resetIbor = 0.00537 - quotedMargin
    currentIbor = 0.027558
    futureIbors = 0.03295

    dm = bond.discountMargin(settlementDate,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)
    
    fullPrice = bond.fullPriceFromDM(settlementDate,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accruedDays
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accruedInterest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlementDate,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollarDuration(settlementDate,
                                       resetIbor,
                                       currentIbor,
                                       futureIbors,
                                       dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modifiedDuration = bond.modifiedRateDuration(settlementDate,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyRateDuration(settlementDate,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDM(settlementDate,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Convexity = ", convexity)

    principal = bond.principal(settlementDate,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Principal = ", principal)

    duration = bond.dollarCreditDuration(settlementDate,
                                         resetIbor,
                                         currentIbor,
                                         futureIbors,
                                         dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modifiedDuration = bond.modifiedCreditDuration(settlementDate,
                                                   resetIbor,
                                                   currentIbor,
                                                   futureIbors,
                                                   dm)

    testCases.print("Modified Credit Duration = ", modifiedDuration)

##########################################################################


test_FinBondFRN()
testCases.compareTestCases()
