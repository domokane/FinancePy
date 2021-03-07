###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import time

import sys
sys.path.append("..")

from financepy.products.credit.FinCDSTranche import FinLossDistributionBuilder
from financepy.products.credit.FinCDSIndexPortfolio import FinCDSIndexPortfolio
from financepy.products.credit.FinCDSTranche import FinCDSTranche
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################

def buildIborCurve(tradeDate):

    valuationDate = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlementDate = valuationDate

    maturityDate = settlementDate.addMonths(12)
    swap1 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(24)
    swap2 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(36)
    swap3 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(48)
    swap4 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(60)
    swap5 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    liborCurve = FinIborSingleCurve(valuationDate, depos, fras, swaps)
    return liborCurve

##############################################################################

def loadHomogeneousCDSCurves(valuationDate,
                             liborCurve,
                             cdsSpread3Y,
                             cdsSpread5Y,
                             cdsSpread7Y,
                             cdsSpread10Y,
                             numCredits):

    maturity3Y = valuationDate.nextCDSDate(36)
    maturity5Y = valuationDate.nextCDSDate(60)
    maturity7Y = valuationDate.nextCDSDate(84)
    maturity10Y = valuationDate.nextCDSDate(120)

    recoveryRate = 0.40

    cds3Y = FinCDS(valuationDate, maturity3Y, cdsSpread3Y)
    cds5Y = FinCDS(valuationDate, maturity5Y, cdsSpread5Y)
    cds7Y = FinCDS(valuationDate, maturity7Y, cdsSpread7Y)
    cds10Y = FinCDS(valuationDate, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuerCurve = FinCDSCurve(valuationDate,
                              contracts,
                              liborCurve,
                              recoveryRate)

    issuerCurves = []
    for _ in range(0, numCredits):
        issuerCurves.append(issuerCurve)

    return issuerCurves

##########################################################################


def loadHeterogeneousSpreadCurves(valuationDate, liborCurve):

    maturity3Y = valuationDate.nextCDSDate(36)
    maturity5Y = valuationDate.nextCDSDate(60)
    maturity7Y = valuationDate.nextCDSDate(84)
    maturity10Y = valuationDate.nextCDSDate(120)
    path = os.path.join(os.path.dirname(__file__), './/data//CDX_NA_IG_S7_SPREADS.csv')
    f = open(path, 'r')
    data = f.readlines()
    f.close()
    issuerCurves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recoveryRate = float(splitRow[5])

        cds3Y = FinCDS(valuationDate, maturity3Y, spd3Y)
        cds5Y = FinCDS(valuationDate, maturity5Y, spd5Y)
        cds7Y = FinCDS(valuationDate, maturity7Y, spd7Y)
        cds10Y = FinCDS(valuationDate, maturity10Y, spd10Y)
        cdsContracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuerCurve = FinCDSCurve(valuationDate,
                                  cdsContracts,
                                  liborCurve,
                                  recoveryRate)

        issuerCurves.append(issuerCurve)

    return issuerCurves

##########################################################################


def test_FinCDSTranche():

    tradeDate = FinDate(1, 3, 2007)
    stepInDate = tradeDate.addDays(1)
    valuationDate = tradeDate.addDays(1)

    testCases.header("DATE")
    testCases.print(str((tradeDate)))
    testCases.print(str((stepInDate)))
    testCases.print(str((valuationDate)))

    liborCurve = buildIborCurve(tradeDate)

    trancheMaturity = FinDate(20, 12, 2011)
    tranche1 = FinCDSTranche(valuationDate, trancheMaturity, 0.00, 0.03)
    tranche2 = FinCDSTranche(valuationDate, trancheMaturity, 0.03, 0.06)
    tranche3 = FinCDSTranche(valuationDate, trancheMaturity, 0.06, 0.09)
    tranche4 = FinCDSTranche(valuationDate, trancheMaturity, 0.09, 0.12)
    tranche5 = FinCDSTranche(valuationDate, trancheMaturity, 0.12, 0.22)
    tranche6 = FinCDSTranche(valuationDate, trancheMaturity, 0.22, 0.60)
    tranche7 = FinCDSTranche(valuationDate, trancheMaturity, 0.00, 0.60)
    tranches = [
        tranche1,
        tranche2,
        tranche3,
        tranche4,
        tranche5,
        tranche6,
        tranche7]

    corr1 = 0.30
    corr2 = 0.35
    upfront = 0.0
    spd = 0.0

    cdsIndex = FinCDSIndexPortfolio()

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "====================== HOMOGENEOUS CURVE ==========================")
    testCases.banner(
        "===================================================================")
    numCredits = 125
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046

    issuerCurves = loadHomogeneousCDSCurves(valuationDate,
                                            liborCurve,
                                            spd3Y,
                                            spd5Y,
                                            spd7Y,
                                            spd10Y,
                                            numCredits)

    intrinsicSpd = cdsIndex.intrinsicSpread(valuationDate,
                                            stepInDate,
                                            trancheMaturity,
                                            issuerCurves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    testCases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    testCases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for numPoints in [40]:
                start = time.time()
                v = tranche.valueBC(
                    valuationDate,
                    issuerCurves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    numPoints,
                    method)
                end = time.time()
                period = (end - start)
                testCases.print(
                    method,
                    period,
                    numPoints,
                    tranche._k1,
                    tranche._k2,
                    v[3] * 10000)

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "=================== HETEROGENEOUS CURVES ==========================")
    testCases.banner(
        "===================================================================")

    issuerCurves = loadHeterogeneousSpreadCurves(valuationDate,
                                                 liborCurve)

    intrinsicSpd = cdsIndex.intrinsicSpread(valuationDate,
                                            stepInDate,
                                            trancheMaturity,
                                            issuerCurves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    testCases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    testCases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for numPoints in [40]:
                start = time.time()
                v = tranche.valueBC(
                    valuationDate,
                    issuerCurves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    numPoints,
                    method)
                end = time.time()
                period = (end - start)
                testCases.print(
                    method,
                    period,
                    numPoints,
                    tranche._k1,
                    tranche._k2,
                    v[3] * 10000)

    testCases.banner(
        "===================================================================")

##########################################################################


test_FinCDSTranche()
testCases.compareTestCases()
