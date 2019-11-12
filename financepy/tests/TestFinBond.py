# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:43:40 2019

@author: Dominic
"""

import datetime as dt

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve
from financepy.market.curves.FinFlatCurve import FinFlatCurve, FinCompoundingMethods
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond, FinYieldConventions
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

    liborCurve = FinLiborOneCurve("USD_LIBOR",
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


def test_FinBond():

    import pandas as pd
    bondDataFrame = pd.read_csv('./data/giltbondprices.txt', sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    settlement = FinDate(2012, 9, 19)

    for accrualType in FinDayCountTypes:

        testCases.header("MATURITY", "COUPON", "CLEAN_PRICE", "ACCD_DAYS",
                         "ACCRUED", "YTM")

        for index, bond in bondDataFrame.iterrows():

            dateString = bond['maturity']
            matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
            maturityDt = FinDate.fromDatetime(matDatetime)
            coupon = bond['coupon']/100.0
            cleanPrice = bond['mid']
            bond = FinBond(maturityDt, coupon, frequencyType, accrualType)

            ytm = bond.yieldToMaturity(settlement, cleanPrice)
            accd = bond._accrued
            accd_days = bond._accruedDays

            testCases.print("%18s" % maturityDt, "%8.4f" % coupon,
                            "%10.4f" % cleanPrice, "%6.0f" % accd_days,
                            "%10.4f" % accd, "%8.4f" % ytm)

    # EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = FinDayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settlementDate = FinDate(1994, 4, 19)
    maturityDate = FinDate(1997, 7, 15)
    coupon = 0.085
    face = 1.0
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    bond = FinBond(maturityDate, coupon, freqType, accrualConvention, face)

    testCases.header("FIELD", "VALUE")
    fullPrice = bond.fullPriceFromYield(settlementDate, y)
    testCases.print("Full Price = ", fullPrice)
    cleanPrice = bond.cleanPriceFromYield(settlementDate, y)
    testCases.print("Clean Price = ", cleanPrice)
    accd = bond._accrued
    testCases.print("Accrued = ", accd)
    ytm = bond.yieldToMaturity(settlementDate, cleanPrice)
    testCases.print("Yield to Maturity = ", ytm)

    bump = 1e-4
    priceBumpedUp = bond.fullPriceFromYield(settlementDate, y + bump)
    testCases.print("Price Bumped Up:", priceBumpedUp)

    priceBumpedDn = bond.fullPriceFromYield(settlementDate, y - bump)
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

    conv = bond.convexityFromYield(settlementDate, y)
    testCases.print("Convexity = ", conv)

    # ASSET SWAP SPREAD

    # When the libor curve is the flat bond curve then the ASW is zero by
    # definition
    flatCurve = FinFlatCurve(settlementDate,
                             ytm,
                             FinCompoundingMethods.SEMI_ANNUAL)

    testCases.header("FIELD", "VALUE")

    cleanPrice = bond.cleanPriceFromYield(settlementDate, ytm)
    asw = bond.assetSwapSpread(settlementDate, cleanPrice, flatCurve)
    testCases.print("Discounted on Bond Curve ASW:", asw * 10000)

    # When the libor curve is the Libor curve then the ASW is positive
    liborCurve = buildLiborCurve(settlementDate)
    asw = bond.assetSwapSpread(settlementDate, cleanPrice, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, cleanPrice, liborCurve)
    testCases.print("Discounted on LIBOR Curve ASW:", asw * 10000)
    testCases.print("Discounted on LIBOR Curve OAS:", oas * 10000)

    p = 0.90
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Deep discount bond at 90 ASW:", asw * 10000)
    testCases.print("Deep discount bond at 90 OAS:", oas * 10000)

    p = 1.00
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Par bond at 100 ASW:", asw * 10000)
    testCases.print("Par bond at 100 OAS:", oas * 10000)

    p = 1.20
    asw = bond.assetSwapSpread(settlementDate, p, liborCurve)
    oas = bond.optionAdjustedSpread(settlementDate, p, liborCurve)
    testCases.print("Above par bond at 120 ASW:", asw * 10000)
    testCases.print("Above par bond at 120 OAS:", oas * 10000)

##########################################################################
# https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
# Page 10 TREASURY NOTE SCREENSHOT
##########################################################################

    testCases.banner("BLOOMBERG US TREASURY EXAMPLE")
    settlementDate = FinDate(2017, 7, 21)
    maturityDate = FinDate(2027, 5, 15)
    coupon = 0.02375
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    face = 100.0

    bond = FinBond(maturityDate,
                   coupon,
                   freqType,
                   accrualType,
                   face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 99.780842

    yld = bond.currentYield(cleanPrice)
    testCases.print("Current Yield = ", yld)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.UK_DMO)
    testCases.print("UK DMO Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.US_STREET)
    testCases.print("US STREET Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity = ", ytm)

    fullPrice = bond.fullPriceFromYield(settlementDate, ytm)
    testCases.print("Full Price = ", fullPrice)

    cleanPrice = bond.cleanPriceFromYield(settlementDate, ytm)
    testCases.print("Clean Price = ", cleanPrice)

    accd = bond._accrued
    testCases.print("Accrued = ", accd)

    accddays = bond._accruedDays
    testCases.print("Accrued Days = ", accddays)

    duration = bond.dollarDuration(settlementDate, ytm)
    testCases.print("Dollar Duration = ", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate, ytm)
    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate, ytm)
    testCases.print("Macauley Duration = ", macauleyDuration)

    conv = bond.convexityFromYield(settlementDate, ytm)
    testCases.print("Convexity = ", conv)

##########################################################################
# Page 11 APPLE NOTE SCREENSHOT
##########################################################################

    testCases.banner("BLOOMBERG APPLE CORP BOND EXAMPLE")
    settlementDate = FinDate(2017, 7, 21)
    maturityDate = FinDate(2022, 5, 13)
    coupon = 0.027
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    face = 100.0

    bond = FinBond(maturityDate,
                   coupon,
                   freqType,
                   accrualType,
                   face)

    testCases.header("FIELD", "VALUE")
    cleanPrice = 101.581564

    yld = bond.currentYield(cleanPrice)
    testCases.print("Current Yield = ", yld)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.UK_DMO)
    testCases.print("UK DMO Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.US_STREET)
    testCases.print("US STREET Yield To Maturity = ", ytm)

    ytm = bond.yieldToMaturity(settlementDate, cleanPrice,
                               FinYieldConventions.US_TREASURY)
    testCases.print("US TREASURY Yield To Maturity = ", ytm)

    fullPrice = bond.fullPriceFromYield(settlementDate, ytm)
    testCases.print("Full Price = ", fullPrice)

    cleanPrice = bond.cleanPriceFromYield(settlementDate, ytm)
    testCases.print("Clean Price = ", cleanPrice)

    # I GET 69 DAYS BUT BBG GETS 68 - CANNOT EXPLAIN!!
    accddays = bond._accruedDays
    testCases.print("Accrued Days = ", accddays)

    accd = bond._accrued
    testCases.print("Accrued = ", accd)

    duration = bond.dollarDuration(settlementDate, ytm)
    testCases.print("Dollar Duration = ", duration)

    modifiedDuration = bond.modifiedDuration(settlementDate, ytm)
    testCases.print("Modified Duration = ", modifiedDuration)

    macauleyDuration = bond.macauleyDuration(settlementDate, ytm)
    testCases.print("Macauley Duration = ", macauleyDuration)

    conv = bond.convexityFromYield(settlementDate, ytm)
    testCases.print("Convexity = ", conv)

##########################################################################


test_FinBond()
testCases.compareTestCases()
