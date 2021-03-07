##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import os
import datetime as dt

import sys
sys.path.append("..")

from financepy.finutils.FinGlobalTypes import FinSwapTypes
from financepy.products.bonds.FinBond import FinYTMCalcType
from financepy.products.bonds.FinBond import FinBond
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.finutils.FinMath import ONE_MILLION
from financepy.finutils.FinDate import FinDate, fromDatetime
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)



##########################################################################

def buildIborCurve(valuationDate):

    depoDCCType = FinDayCountTypes.THIRTY_E_360_ISDA
    depos = []
    depositRate = 0.050

    depo0 = FinIborDeposit(
        valuationDate,
        "1D",
        depositRate,
        depoDCCType)

    spotDays = 2
    settlementDate = valuationDate.addWeekDays(spotDays)

    maturityDate = settlementDate.addMonths(1)
    depo1 = FinIborDeposit(settlementDate,
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

    depos.append(depo0)
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
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    
#    print(swap1._fixedLeg._paymentDates)

    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap2)

#    print(swap2._fixedLeg._paymentDates)

    
    maturityDate = settlementDate.addMonths(48)
    swap3 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap3)

#    print(swap3._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap4)

#    print(swap4._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap5)

#    print(swap5._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap6)

#    print(swap6._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap7)

#    print(swap7._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap8)

#    print(swap8._fixedLeg._paymentDates)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        swapRate,
        fixedFreqType,
        fixedDCCType)
    swaps.append(swap9)

#    print(swap9._fixedLeg._paymentDates)

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


def test_FinBond():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freqType = FinFrequencyTypes.SEMI_ANNUAL
    settlementDate = FinDate(19, 9, 2012)
    face = ONE_MILLION

    for accrualType in FinDayCountTypes:

        testCases.header("MATURITY", "COUPON", "CLEAN_PRICE", "ACCD_DAYS",
                         "ACCRUED", "YTM")

        for _, bond in bondDataFrame.iterrows():

            dateString = bond['maturity']
            matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
            maturityDt = fromDatetime(matDatetime)
            issueDt = FinDate(maturityDt._d, maturityDt._m, 2000)

            coupon = bond['coupon']/100.0
            cleanPrice = bond['mid']
            bond = FinBond(issueDt, maturityDt,
                           coupon, freqType, accrualType, 100)

            ytm = bond.yieldToMaturity(settlementDate, cleanPrice)
            accd = bond._accruedInterest
            accd_days = bond._accruedDays

            testCases.print("%18s" % maturityDt, "%8.4f" % coupon,
                            "%10.4f" % cleanPrice, "%6.0f" % accd_days,
                            "%10.4f" % accd, "%8.4f" % ytm)

    ###########################################################################
    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = FinDayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settlementDate = FinDate(19, 4, 1994)
    issueDate = FinDate(15, 7, 1990)
    maturityDate = FinDate(15, 7, 1997)
    coupon = 0.085
    face = ONE_MILLION
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    bond = FinBond(issueDate, maturityDate,
                   coupon, freqType, accrualConvention, face)

    testCases.header("FIELD", "VALUE")
    fullPrice = bond.fullPriceFromYTM(settlementDate, y)
    testCases.print("Full Price = ", fullPrice)
    cleanPrice = bond.cleanPriceFromYTM(settlementDate, y)
    testCases.print("Clean Price = ", cleanPrice)
    accd = bond._accruedInterest
    testCases.print("Accrued = ", accd)
    ytm = bond.yieldToMaturity(settlementDate, cleanPrice)
    testCases.print("Yield to Maturity = ", ytm)

    bump = 1e-4
    priceBumpedUp = bond.fullPriceFromYTM(settlementDate, y + bump)
    testCases.print("Price Bumped Up:", priceBumpedUp)

    priceBumpedDn = bond.fullPriceFromYTM(settlementDate, y - bump)
    testCases.print("Price Bumped Dn:", priceBumpedDn)

    durationByBump = -(priceBumpedUp - fullPrice) / bump
    testCases.print("Duration by Bump = ", durationByBump)

    duration = bond.dollarDuration(settlementDate, y)
    testCases.print("Dollar Duration = ", duration)
    testCases.print("Duration Difference:", duration - durationByBump)

    modifiedDuration = bond.modifiedDuration(settlementDate, y)
    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate, y)
    testCases.print("Macauley Duration = ", macauleyDuration)

    conv = bond.convexityFromYTM(settlementDate, y)
    testCases.print("Convexity = ", conv)

    # ASSET SWAP SPREAD

    # When the libor curve is the flat bond curve then the ASW is zero by
    # definition
    flatCurve = FinDiscountCurveFlat(settlementDate,
                                     ytm,
                                     FinFrequencyTypes.SEMI_ANNUAL)

    testCases.header("FIELD", "VALUE")

    cleanPrice = bond.cleanPriceFromYTM(settlementDate, ytm)
    asw = bond.assetSwapSpread(settlementDate, cleanPrice, flatCurve)
    testCases.print("Discounted on Bond Curve ASW:", asw * 10000)

    # When the libor curve is the Libor curve then the ASW is positive
    liborCurve = buildIborCurve(settlementDate)
    asw = bond.assetSwapSpread(settlementDate, cleanPrice, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, cleanPrice, liborCurve)
    testCases.print("Discounted on LIBOR Curve ASW:", asw * 10000)
    testCases.print("Discounted on LIBOR Curve OAS:", oas * 10000)

    p = 90.0
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Deep discount bond at 90 ASW:", asw * 10000)
    testCases.print("Deep discount bond at 90 OAS:", oas * 10000)

    p = 100.0
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Par bond at 100 ASW:", asw * 10000)
    testCases.print("Par bond at 100 OAS:", oas * 10000)

    p = 120.0
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Above par bond at 120 ASW:", asw * 10000)
    testCases.print("Above par bond at 120 OAS:", oas * 10000)

##########################################################################
# https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
# Page 10 TREASURY NOTE SCREENSHOT
##########################################################################

    testCases.banner("BLOOMBERG US TREASURY EXAMPLE")
    settlementDate = FinDate(21, 7, 2017)
    issueDate = FinDate(15, 5, 2010)
    maturityDate = FinDate(15, 5, 2027)
    coupon = 0.02375
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    face = 100.0

    bond = FinBond(issueDate,
                   maturityDate,
                   coupon,
                   freqType,
                   accrualType,
                   face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 99.7808417

    yld = bond.currentYield(cleanPrice)
    testCases.print("Current Yield = ", yld)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.UK_DMO)
    testCases.print("UK DMO Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.US_STREET)
    testCases.print("US STREET Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity = ", ytm)

    fullPrice = bond.fullPriceFromYTM(settlementDate, ytm)
    testCases.print("Full Price = ", fullPrice)

    cleanPrice = bond.cleanPriceFromYTM(settlementDate, ytm)
    testCases.print("Clean Price = ", cleanPrice)

    accd = bond._accruedInterest
    testCases.print("Accrued = ", accd)

    accddays = bond._accruedDays
    testCases.print("Accrued Days = ", accddays)

    duration = bond.dollarDuration(settlementDate, ytm)
    testCases.print("Dollar Duration = ", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate, ytm)
    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate, ytm)
    testCases.print("Macauley Duration = ", macauleyDuration)

    conv = bond.convexityFromYTM(settlementDate, ytm)
    testCases.print("Convexity = ", conv)

##########################################################################
# Page 11 APPLE NOTE SCREENSHOT
##########################################################################

    testCases.banner("BLOOMBERG APPLE CORP BOND EXAMPLE")
    settlementDate = FinDate(21, 7, 2017)
    issueDate = FinDate(13, 5, 2012)
    maturityDate = FinDate(13, 5, 2022)
    coupon = 0.027
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.THIRTY_E_360_ISDA
    face = 100.0

    bond = FinBond(issueDate, maturityDate,
                   coupon, freqType, accrualType, face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 101.581564

    yld = bond.currentYield(cleanPrice)
    testCases.print("Current Yield", yld)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.UK_DMO)
    testCases.print("UK DMO Yield To Maturity", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.US_STREET)
    testCases.print("US STREET Yield To Maturity", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYTMCalcType.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity", ytm)

    fullPrice = bond.fullPriceFromYTM(settlementDate, ytm)
    testCases.print("Full Price", fullPrice)

    cleanPrice = bond.cleanPriceFromYTM(settlementDate, ytm)
    testCases.print("Clean Price", cleanPrice)

    accddays = bond._accruedDays
    testCases.print("Accrued Days", accddays)

    accd = bond._accruedInterest
    testCases.print("Accrued", accd)

    duration = bond.dollarDuration(settlementDate, ytm)
    testCases.print("Dollar Duration", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate, ytm)
    testCases.print("Modified Duration", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate, ytm)
    testCases.print("Macauley Duration", macauleyDuration)

    conv = bond.convexityFromYTM(settlementDate, ytm)
    testCases.print("Convexity", conv)

###############################################################################


def test_FinBondExDividend():

    issueDate = FinDate(7, 9, 2000)
    maturityDate = FinDate(7, 9, 2020)
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    face = 100.0
    exDivDays = 7
    testCases.header("LABEL", "VALUE")

    calendarType = FinCalendarTypes.UNITED_KINGDOM
    bond = FinBond(issueDate, maturityDate, coupon,
                   freqType, accrualType, face)
    settlementDate = FinDate(7, 9, 2003)
    accrued = bond.calcAccruedInterest(settlementDate, exDivDays, calendarType)
    testCases.print("SettlementDate:", settlementDate)
    testCases.print("Accrued:", accrued)

    ###########################################################################
    testCases.banner("=======================================================")
    testCases.header("SETTLEMENT", "ACCRUED")

    issueDate = FinDate(7, 9, 2000)
    maturityDate = FinDate(7, 9, 2020)
    coupon = 0.05
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    face = 100.0
    exDivDays = 7

    calendarType = FinCalendarTypes.UNITED_KINGDOM
    bond = FinBond(issueDate, maturityDate, coupon,
                   freqType, accrualType, face)

    settlementDate = FinDate(25, 8, 2010)

    for _ in range(0, 13):
        settlementDate = settlementDate.addDays(1)
        accrued = bond.calcAccruedInterest(
            settlementDate, exDivDays, calendarType)
        testCases.print(settlementDate, accrued)

###############################################################################


test_FinBond()
test_FinBondExDividend()
testCases.compareTestCases()
