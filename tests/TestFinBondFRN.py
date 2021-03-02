###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Date import *
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.bonds.BondFRN import FinBondFRN
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def buildIborCurve(valuation_date):

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []

    payFixed = FinSwapTypes.PAY

    spotDays = 2
    settlement_date = valuation_date.addWeekDays(spotDays)

    depositRate = 0.050
    maturity_date = settlement_date.addMonths(1)
    depo1 = FinIborDeposit(
        settlement_date,
        maturity_date,
        depositRate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(3)
    depo2 = FinIborDeposit(
        settlement_date,
        maturity_date,
        depositRate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(6)
    depo3 = FinIborDeposit(
        settlement_date,
        maturity_date,
        depositRate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(9)
    depo4 = FinIborDeposit(
        settlement_date,
        maturity_date,
        depositRate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(12)
    depo5 = FinIborDeposit(
        settlement_date,
        maturity_date,
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
    maturity_date = settlement_date.addMonths(24)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(36)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(48)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(60)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(72)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    maturity_date = settlement_date.addMonths(84)
    swap6 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    maturity_date = settlement_date.addMonths(96)
    swap7 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    maturity_date = settlement_date.addMonths(108)
    swap8 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    maturity_date = settlement_date.addMonths(120)
    swap9 = FinIborSwap(
        settlement_date,
        maturity_date,
        swapRate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    libor_curve = FinIborSingleCurve(valuation_date,
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
            df1 = libor_curve.df(t)
            fwd = (df0 / df1 - 1.0) / dt
            print(t, df1, fwd)
            df0 = df1

    return libor_curve

##########################################################################


def test_FinBondFRN():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # I have a day out problem on the accrued interest - should be 71 and not 72 days
    # Other than that agreement on the DM is very good.

    ##########################################################################
    # CITIGROUP FRN SCREENSHOT
    ##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE")
    issue_date = Date(10, 11, 2010)
    maturity_date = Date(10, 11, 2021)
    quotedMargin = 0.0025
    freq_type = FinFrequencyTypes.QUARTERLY
    accrual_type = FinDayCountTypes.THIRTY_E_360
    face = 1000000

    bond = FinBondFRN(issue_date,
                      maturity_date,
                      quotedMargin,
                      freq_type,
                      accrual_type,
                      face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 96.793
    resetIbor = 0.0143456 - quotedMargin
    currentIbor = 0.0120534
    futureIbors = 0.0130522

    settlement_date = Date(21, 7, 2017)

    dm = bond.discountMargin(settlement_date,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    fullPrice = bond.fullPriceFromDM(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accruedInterest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollarDuration(settlement_date,
                                   resetIbor,
                                   currentIbor,
                                   futureIbors,
                                   dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modifiedDuration = bond.modifiedRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDM(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Convexity = ", convexity)

    duration = bond.dollarCreditDuration(settlement_date,
                                         resetIbor,
                                         currentIbor,
                                         futureIbors,
                                         dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modifiedDuration = bond.modifiedCreditDuration(settlement_date,
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
    issue_date = Date(28, 3, 2000)
    settlement_date = Date(28, 3, 2014)
    maturity_date = Date(3, 2, 2021)
    quotedMargin = 0.0020
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    accrual_type = FinDayCountTypes.THIRTY_E_360_ISDA
    face = 1000000.0

    bond = FinBondFRN(issue_date,
                      maturity_date,
                      quotedMargin,
                      freq_type,
                      accrual_type,
                      face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 93.08
    resetIbor = 0.00537 - quotedMargin
    currentIbor = 0.027558
    futureIbors = 0.03295

    dm = bond.discountMargin(settlement_date,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             cleanPrice)

    testCases.print("Discount Margin (bp) = ", dm * 10000)
    
    fullPrice = bond.fullPriceFromDM(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", fullPrice)

    lastCouponDt = bond._pcd
    testCases.print("Last Coupon Date = ", str(lastCouponDt))

    accddays = bond._accrued_days
    testCases.print("Accrued Days = ", accddays)

    accdAmount = bond._accruedInterest
    testCases.print("Accrued Amount = ", accdAmount)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Dollar Principal = ", principal)

    duration = bond.dollarDuration(settlement_date,
                                       resetIbor,
                                       currentIbor,
                                       futureIbors,
                                       dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modifiedDuration = bond.modifiedRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauleyDuration)

    convexity = bond.convexityFromDM(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Convexity = ", convexity)

    principal = bond.principal(settlement_date,
                               resetIbor,
                               currentIbor,
                               futureIbors,
                               dm)

    testCases.print("Principal = ", principal)

    duration = bond.dollarCreditDuration(settlement_date,
                                         resetIbor,
                                         currentIbor,
                                         futureIbors,
                                         dm)

    testCases.print("Dollar Credit Duration = ", duration)

    modifiedDuration = bond.modifiedCreditDuration(settlement_date,
                                                   resetIbor,
                                                   currentIbor,
                                                   futureIbors,
                                                   dm)

    testCases.print("Modified Credit Duration = ", modifiedDuration)

##########################################################################


test_FinBondFRN()
testCases.compareTestCases()
