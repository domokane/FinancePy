# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.credit.FinCDSIndexPortfolio import FinCDSIndexPortfolio
from financepy.products.credit.FinCDSBasket import FinCDSBasket
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve
from financepy.market.curves.FinCDSCurve import FinCDSCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinMath import corrMatrixGenerator
from financepy.finutils.FinDate import FinDate
import time
import numpy as np
from os.path import dirname, join
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


def buildLiborCurve(tradeDate):

    valuationDate = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlementDate = valuationDate

    maturityDate = settlementDate.addMonths(12)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(24)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(36)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(48)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(60)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    liborCurve = FinLiborOneCurve(
        "USD_LIBOR", settlementDate, depos, fras, swaps)

    return liborCurve

##########################################################################


def loadHomogeneousSpreadCurves(valuationDate,
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
    for iCredit in range(0, numCredits):
        issuerCurves.append(issuerCurve)

    return issuerCurves

##########################################################################


def loadHeterogeneousSpreadCurves(valuationDate, liborCurve):

    maturity3Y = valuationDate.nextCDSDate(36)
    maturity5Y = valuationDate.nextCDSDate(60)
    maturity7Y = valuationDate.nextCDSDate(84)
    maturity10Y = valuationDate.nextCDSDate(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, 'r')

    data = f.readlines()
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


def test_FinCDSBasket():

    tradeDate = FinDate(2007, 3, 1)
    stepInDate = tradeDate.addDays(1)
    valuationDate = tradeDate.addDays(1)

    liborCurve = buildLiborCurve(tradeDate)

    basketMaturity = FinDate(2011, 12, 20)

    cdsIndex = FinCDSIndexPortfolio()

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "====================== INHOMOGENEOUS CURVE ==========================")
    testCases.banner(
        "===================================================================")

    numCredits = 5
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046

    testCases.header("LABELS", "VALUE")

    if 1 == 0:
        issuerCurves = loadHomogeneousSpreadCurves(valuationDate,
                                                   liborCurve,
                                                   spd3Y,
                                                   spd5Y,
                                                   spd7Y,
                                                   spd10Y,
                                                   numCredits)
    else:
        issuerCurves = loadHeterogeneousSpreadCurves(valuationDate, liborCurve)
        issuerCurves = issuerCurves[0:numCredits]

    intrinsicSpd = cdsIndex.intrinsicSpread(valuationDate,
                                            stepInDate,
                                            basketMaturity,
                                            issuerCurves) * 10000.0

    testCases.print("INTRINSIC SPD BASKET MATURITY", intrinsicSpd)

    totalSpd = cdsIndex.totalSpread(valuationDate,
                                    stepInDate,
                                    basketMaturity,
                                    issuerCurves) * 10000.0

    testCases.print("SUMMED UP SPD BASKET MATURITY", totalSpd)

    minSpd = cdsIndex.minSpread(valuationDate,
                                stepInDate,
                                basketMaturity,
                                issuerCurves) * 10000.0

    testCases.print("MINIMUM SPD BASKET MATURITY", minSpd)

    maxSpd = cdsIndex.maxSpread(valuationDate,
                                stepInDate,
                                basketMaturity,
                                issuerCurves) * 10000.0

    testCases.print("MAXIMUM SPD BASKET MATURITY", maxSpd)

    seed = 1967
    basket = FinCDSBasket(valuationDate,
                          basketMaturity)

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "======================= GAUSSIAN COPULA ===========================")
    testCases.banner(
        "===================================================================")

    testCases.header("TIME", "Trials", "RHO", "NTD", "SPRD", "SPRD_HOMO")

    for ntd in range(1, numCredits + 1):
        for beta in [0.0, 0.5]:
            rho = beta * beta
            betaVector = np.ones(numCredits) * beta
            corrMatrix = corrMatrixGenerator(rho, numCredits)
            for numTrials in [1000]:  # [1000,5000,10000,20000,50000,100000]:
                start = time.time()

                v1 = basket.valueGaussian_MC(valuationDate,
                                             ntd,
                                             issuerCurves,
                                             corrMatrix,
                                             liborCurve,
                                             numTrials,
                                             seed)

                v2 = basket.value1FGaussian_Homo(valuationDate,
                                                 ntd,
                                                 issuerCurves,
                                                 betaVector,
                                                 liborCurve)

                end = time.time()
                period = (end - start)
                testCases.print(
                    period,
                    numTrials,
                    rho,
                    ntd,
                    v1[2] * 10000,
                    v2[3] * 10000)

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "==================== STUDENT'S-T CONVERGENCE ======================")
    testCases.banner(
        "===================================================================")

    testCases.header("TIME", "TRIALS", "RHO", "DOF", "NTD", "SPRD")

    for beta in [0.0, 0.5]:
        rho = beta ** 2
        corrMatrix = corrMatrixGenerator(rho, numCredits)
        for ntd in range(1, numCredits + 1):
            for doF in [3, 10]:
                start = time.time()

                v = basket.valueStudentT_MC(valuationDate,
                                            ntd,
                                            issuerCurves,
                                            corrMatrix,
                                            doF,
                                            liborCurve,
                                            numTrials,
                                            seed)

                end = time.time()
                period = (end - start)
                testCases.print(period, numTrials, rho, doF, ntd, v[2] * 10000)

            start = time.time()
            v = basket.valueGaussian_MC(
                valuationDate,
                ntd,
                issuerCurves,
                corrMatrix,
                liborCurve,
                numTrials,
                seed)
            end = time.time()
            period = (end - start)

            testCases.print(period, numTrials, rho, "GC", ntd, v[2] * 10000)

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "=================== STUDENT'S T WITH DOF = 5 ======================")
    testCases.banner(
        "===================================================================")
    doF = 5
    testCases.header("TIME", "NUMTRIALS", "RHO", "NTD", "SPD")
    for beta in [0.0, 0.5]:
        rho = beta ** 2
        corrMatrix = corrMatrixGenerator(rho, numCredits)
        for ntd in range(1, numCredits + 1):
            for numTrials in [1000]:
                start = time.time()

                v = basket.valueStudentT_MC(valuationDate,
                                            ntd,
                                            issuerCurves,
                                            corrMatrix,
                                            doF,
                                            liborCurve,
                                            numTrials,
                                            seed)
                end = time.time()
                period = (end - start)
                testCases.print(period, numTrials, rho, ntd, v[2] * 10000)

##########################################################################


test_FinCDSBasket()
testCases.compareTestCases()
