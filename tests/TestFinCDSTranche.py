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
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.Date import Date
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################

def buildIborCurve(tradeDate):

    valuation_date = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlement_date = valuation_date

    maturity_date = settlement_date.addMonths(12)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(24)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(36)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(48)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(60)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    libor_curve = FinIborSingleCurve(valuation_date, depos, fras, swaps)
    return libor_curve

##############################################################################

def loadHomogeneousCDSCurves(valuation_date,
                             libor_curve,
                             cdsSpread3Y,
                             cdsSpread5Y,
                             cdsSpread7Y,
                             cdsSpread10Y,
                             numCredits):

    maturity3Y = valuation_date.nextCDSDate(36)
    maturity5Y = valuation_date.nextCDSDate(60)
    maturity7Y = valuation_date.nextCDSDate(84)
    maturity10Y = valuation_date.nextCDSDate(120)

    recovery_rate = 0.40

    cds3Y = FinCDS(valuation_date, maturity3Y, cdsSpread3Y)
    cds5Y = FinCDS(valuation_date, maturity5Y, cdsSpread5Y)
    cds7Y = FinCDS(valuation_date, maturity7Y, cdsSpread7Y)
    cds10Y = FinCDS(valuation_date, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuer_curve = FinCDSCurve(valuation_date,
                              contracts,
                              libor_curve,
                              recovery_rate)

    issuer_curves = []
    for _ in range(0, numCredits):
        issuer_curves.append(issuer_curve)

    return issuer_curves

##########################################################################


def loadHeterogeneousSpreadCurves(valuation_date, libor_curve):

    maturity3Y = valuation_date.nextCDSDate(36)
    maturity5Y = valuation_date.nextCDSDate(60)
    maturity7Y = valuation_date.nextCDSDate(84)
    maturity10Y = valuation_date.nextCDSDate(120)
    path = os.path.join(os.path.dirname(__file__), './/data//CDX_NA_IG_S7_SPREADS.csv')
    f = open(path, 'r')
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = FinCDS(valuation_date, maturity3Y, spd3Y)
        cds5Y = FinCDS(valuation_date, maturity5Y, spd5Y)
        cds7Y = FinCDS(valuation_date, maturity7Y, spd7Y)
        cds10Y = FinCDS(valuation_date, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = FinCDSCurve(valuation_date,
                                  cds_contracts,
                                  libor_curve,
                                  recovery_rate)

        issuer_curves.append(issuer_curve)

    return issuer_curves

##########################################################################


def test_FinCDSTranche():

    tradeDate = Date(1, 3, 2007)
    step_in_date = tradeDate.addDays(1)
    valuation_date = tradeDate.addDays(1)

    testCases.header("DATE")
    testCases.print(str((tradeDate)))
    testCases.print(str((step_in_date)))
    testCases.print(str((valuation_date)))

    libor_curve = buildIborCurve(tradeDate)

    trancheMaturity = Date(20, 12, 2011)
    tranche1 = FinCDSTranche(valuation_date, trancheMaturity, 0.00, 0.03)
    tranche2 = FinCDSTranche(valuation_date, trancheMaturity, 0.03, 0.06)
    tranche3 = FinCDSTranche(valuation_date, trancheMaturity, 0.06, 0.09)
    tranche4 = FinCDSTranche(valuation_date, trancheMaturity, 0.09, 0.12)
    tranche5 = FinCDSTranche(valuation_date, trancheMaturity, 0.12, 0.22)
    tranche6 = FinCDSTranche(valuation_date, trancheMaturity, 0.22, 0.60)
    tranche7 = FinCDSTranche(valuation_date, trancheMaturity, 0.00, 0.60)
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

    issuer_curves = loadHomogeneousCDSCurves(valuation_date,
                                            libor_curve,
                                            spd3Y,
                                            spd5Y,
                                            spd7Y,
                                            spd10Y,
                                            numCredits)

    intrinsicSpd = cdsIndex.intrinsicSpread(valuation_date,
                                            step_in_date,
                                            trancheMaturity,
                                            issuer_curves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    testCases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    testCases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for num_points in [40]:
                start = time.time()
                v = tranche.valueBC(
                    valuation_date,
                    issuer_curves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    num_points,
                    method)
                end = time.time()
                period = (end - start)
                testCases.print(
                    method,
                    period,
                    num_points,
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

    issuer_curves = loadHeterogeneousSpreadCurves(valuation_date,
                                                 libor_curve)

    intrinsicSpd = cdsIndex.intrinsicSpread(valuation_date,
                                            step_in_date,
                                            trancheMaturity,
                                            issuer_curves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    testCases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    testCases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for num_points in [40]:
                start = time.time()
                v = tranche.valueBC(
                    valuation_date,
                    issuer_curves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    num_points,
                    method)
                end = time.time()
                period = (end - start)
                testCases.print(
                    method,
                    period,
                    num_points,
                    tranche._k1,
                    tranche._k2,
                    v[3] * 10000)

    testCases.banner(
        "===================================================================")

##########################################################################


test_FinCDSTranche()
testCases.compareTestCases()
