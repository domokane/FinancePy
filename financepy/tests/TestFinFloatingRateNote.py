# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:43:40 2019

@author: Dominic
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.products.bonds.FinFloatingRateNote import FinFloatingRateNote
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def buildLiborCurve(valuationDate):

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    spotDays = 2
    settlementDate = valuationDate.addWorkDays(spotDays)

    depositRate = 0.050
    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(3)
    depo2 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(6)
    depo3 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(9)
    depo4 = FinLiborDeposit(
        settlementDate,
        maturityDate,
        depositRate,
        depoDCCType)

    maturityDate = settlementDate.addMonths(12)
    depo5 = FinLiborDeposit(
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
    fixedDCCType = FinDayCountTypes.ACT_365_ISDA
    fixedFreqType = FinFrequencyTypes.SEMI_ANNUAL

    swaps = []

    swapRate = 0.05
    maturityDate = settlementDate.addMonths(24)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(
        settlementDate,
        maturityDate,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    liborCurve = FinLiborCurve("USD_LIBOR",
                                  settlementDate,
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


def test_FinFloatingRateNote():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # I have a day out problem on the accrued interest - should be 71 and not 72 days
    # Other than that agreement on the DM is very good.

    ##########################################################################
    # CITIGROUP FRN SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE")
    settlementDate = FinDate(2017, 7, 21)
    maturityDate = FinDate(2021, 11, 10)
    quotedMargin = 0.0025
    freqType = FinFrequencyTypes.QUARTERLY
    accrualType = FinDayCountTypes.THIRTY_360
    face = 100.0
    redemption = 1.0

    bond = FinFloatingRateNote(maturityDate,
                               quotedMargin,
                               freqType,
                               accrualType,
                               face,
                               redemption)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 96.793
    nextCoupon = 0.0143456
    futureLibors = 0.0130522

    dm = bond.discountMargin(settlementDate,
                             nextCoupon,
                             futureLibors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    fullPrice = bond.fullPriceFromDiscountMargin(settlementDate, nextCoupon,
                                                 futureLibors, dm)
    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond.pcd(settlementDate)
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond.accruedDays(settlementDate)
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond.accruedInterest(settlementDate, nextCoupon)
    testCases.print("Accrued Amount = ", accdAmount)

    duration = bond.dollarDuration(settlementDate,
                                   nextCoupon,
                                   futureLibors,
                                   dm)

    testCases.print("Dollar Duration = ", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate,
                                             nextCoupon,
                                             futureLibors,
                                             dm)

    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate,
                                             nextCoupon,
                                             futureLibors,
                                             dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDiscountMargin(settlementDate,
                                                 nextCoupon,
                                                 futureLibors,
                                                 dm)

    testCases.print("Convexity = ", convexity)

##########################################################################
# EXAMPLE
# https://ebrary.net/14293/economics/actual_floater
##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE II")
    settlementDate = FinDate(2014, 3, 28)
    maturityDate = FinDate(2021, 2, 3)
    quotedMargin = 0.0020
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.THIRTY_360
    face = 100.0
    redemption = 1.0

    bond = FinFloatingRateNote(maturityDate,
                               quotedMargin,
                               freqType,
                               accrualType,
                               face,
                               redemption)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 93.08
    nextCoupon = 0.00537
    futureLibors = 0.027558

    dm = bond.discountMargin(settlementDate,
                             nextCoupon,
                             futureLibors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    fullPrice = bond.fullPriceFromDiscountMargin(settlementDate, nextCoupon,
                                                 futureLibors, dm)
    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond.pcd(settlementDate)
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond.accruedDays(settlementDate)
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond.accruedInterest(settlementDate, nextCoupon)
    testCases.print("Accrued Amount = ", accdAmount)

    duration = bond.dollarDuration(settlementDate,
                                   nextCoupon,
                                   futureLibors,
                                   dm)

    testCases.print("Dollar Duration = ", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate,
                                             nextCoupon,
                                             futureLibors,
                                             dm)

    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate,
                                             nextCoupon,
                                             futureLibors,
                                             dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDiscountMargin(settlementDate,
                                                 nextCoupon,
                                                 futureLibors,
                                                 dm)

    testCases.print("Convexity = ", convexity)


##########################################################################
test_FinFloatingRateNote()
testCases.compareTestCases()
