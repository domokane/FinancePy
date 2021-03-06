###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.date import *
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.products.bonds.floating_rate_note import FloatingRateNote
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def buildIborCurve(valuation_date):

    depoDCCType = DayCountTypes.THIRTY_E_360_ISDA
    depos = []

    payFixed = FinSwapTypes.PAY

    spotDays = 2
    settlement_date = valuation_date.addWeekDays(spotDays)

    deposit_rate = 0.050
    maturity_date = settlement_date.addMonths(1)
    depo1 = FinIborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(3)
    depo2 = FinIborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(6)
    depo3 = FinIborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(9)
    depo4 = FinIborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    maturity_date = settlement_date.addMonths(12)
    depo5 = FinIborDeposit(
        settlement_date,
        maturity_date,
        deposit_rate,
        depoDCCType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []
    fixedDCCType = DayCountTypes.ACT_365F
    fixedFreqType = FrequencyTypes.SEMI_ANNUAL

    swaps = []

    swap_rate = 0.05
    maturity_date = settlement_date.addMonths(24)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(36)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(48)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(60)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(72)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

    maturity_date = settlement_date.addMonths(84)
    swap6 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

    maturity_date = settlement_date.addMonths(96)
    swap7 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

    maturity_date = settlement_date.addMonths(108)
    swap8 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

    maturity_date = settlement_date.addMonths(120)
    swap9 = FinIborSwap(
        settlement_date,
        maturity_date,
        swap_rate,
        payFixed,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

    libor_curve = IborSingleCurve(valuation_date,
                                  depos,
                                  fras,
                                  swaps)
    
    if 1 == 0:
        import numpy as np
        num_steps = 40
        dt = 10 / num_steps
        times = np.linspace(0.0, 10.0, num_steps + 1)

        df0 = 1.0
        for t in times[1:]:
            df1 = libor_curve.df(t)
            fwd = (df0 / df1 - 1.0) / dt
            print(t, df1, fwd)
            df0 = df1

    return libor_curve

##########################################################################


def test_BondFRN():

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
    freq_type = FrequencyTypes.QUARTERLY
    accrual_type = DayCountTypes.THIRTY_E_360
    face = 1000000

    bond = FloatingRateNote(issue_date,
                            maturity_date,
                            quotedMargin,
                            freq_type,
                            accrual_type,
                            face)

    testCases.header("FIELD", "VALUE")
    clean_price = 96.793
    resetIbor = 0.0143456 - quotedMargin
    currentIbor = 0.0120534
    futureIbors = 0.0130522

    settlement_date = Date(21, 7, 2017)

    dm = bond.discountMargin(settlement_date,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             clean_price)

    testCases.print("Discount Margin (bp) = ", dm * 10000)

    full_price = bond.full_price_from_dm(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", full_price)

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

    duration = bond.dollar_duration(settlement_date,
                                   resetIbor,
                                   currentIbor,
                                   futureIbors,
                                   dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modifiedRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauleyRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauley_duration)

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

    modified_duration = bond.modifiedCreditDuration(settlement_date,
                                                   resetIbor,
                                                   currentIbor,
                                                   futureIbors,
                                                   dm)

    testCases.print("Modified Credit Duration = ", modified_duration)

##########################################################################
# EXAMPLE
# https://ebrary.net/14293/economics/actual_floater
##########################################################################

    testCases.banner("BLOOMBERG CITIGROUP FRN EXAMPLE II")
    issue_date = Date(28, 3, 2000)
    settlement_date = Date(28, 3, 2014)
    maturity_date = Date(3, 2, 2021)
    quotedMargin = 0.0020
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 1000000.0

    bond = FloatingRateNote(issue_date,
                            maturity_date,
                            quotedMargin,
                            freq_type,
                            accrual_type,
                            face)

    testCases.header("FIELD", "VALUE")
    clean_price = 93.08
    resetIbor = 0.00537 - quotedMargin
    currentIbor = 0.027558
    futureIbors = 0.03295

    dm = bond.discountMargin(settlement_date,
                             resetIbor,
                             currentIbor,
                             futureIbors,
                             clean_price)

    testCases.print("Discount Margin (bp) = ", dm * 10000)
    
    full_price = bond.full_price_from_dm(settlement_date,
                                     resetIbor,
                                     currentIbor,
                                     futureIbors,
                                     dm)

    testCases.print("Full Price = ", full_price)

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

    duration = bond.dollar_duration(settlement_date,
                                       resetIbor,
                                       currentIbor,
                                       futureIbors,
                                       dm)

    testCases.print("Dollar Rate Duration = ", duration)

    modified_duration = bond.modifiedRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Modified Rate Duration = ", modified_duration)

    macauley_duration = bond.macauleyRateDuration(settlement_date,
                                                 resetIbor,
                                                 currentIbor,
                                                 futureIbors,
                                                 dm)

    testCases.print("Macauley Duration = ", macauley_duration)

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

    modified_duration = bond.modifiedCreditDuration(settlement_date,
                                                   resetIbor,
                                                   currentIbor,
                                                   futureIbors,
                                                   dm)

    testCases.print("Modified Credit Duration = ", modified_duration)

##########################################################################


test_BondFRN()
testCases.compareTestCases()
