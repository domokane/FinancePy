###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from os.path import dirname, join

import sys
sys.path.append("..")

from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.products.credit.cds_basket import CDSBasket
from financepy.products.credit.cds import CDS
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.products.credit.cds_curve import CDSCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.math import corrMatrixGenerator
from financepy.utils.date import Date
from financepy.models.gbm_process_simulator import getPathsAssets
from financepy.utils.global_types import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


def buildIborCurve(tradeDate):

    valuation_date = tradeDate.addDays(1)
    dcType = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
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

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve

##########################################################################


def loadHomogeneousSpreadCurves(valuation_date,
                                libor_curve,
                                cdsSpread3Y,
                                cdsSpread5Y,
                                cdsSpread7Y,
                                cdsSpread10Y,
                                num_credits):

    maturity3Y = valuation_date.nextCDSDate(36)
    maturity5Y = valuation_date.nextCDSDate(60)
    maturity7Y = valuation_date.nextCDSDate(84)
    maturity10Y = valuation_date.nextCDSDate(120)

    recovery_rate = 0.40

    cds3Y = CDS(valuation_date, maturity3Y, cdsSpread3Y)
    cds5Y = CDS(valuation_date, maturity5Y, cdsSpread5Y)
    cds7Y = CDS(valuation_date, maturity7Y, cdsSpread7Y)
    cds10Y = CDS(valuation_date, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuer_curve = CDSCurve(valuation_date,
                            contracts,
                            libor_curve,
                            recovery_rate)

    issuer_curves = []
    for _ in range(0, num_credits):
        issuer_curves.append(issuer_curve)

    return issuer_curves

##########################################################################


def loadHeterogeneousSpreadCurves(valuation_date, libor_curve):

    maturity3Y = valuation_date.nextCDSDate(36)
    maturity5Y = valuation_date.nextCDSDate(60)
    maturity7Y = valuation_date.nextCDSDate(84)
    maturity10Y = valuation_date.nextCDSDate(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, 'r')

    data = f.readlines()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = CDS(valuation_date, maturity3Y, spd3Y)
        cds5Y = CDS(valuation_date, maturity5Y, spd5Y)
        cds7Y = CDS(valuation_date, maturity7Y, spd7Y)
        cds10Y = CDS(valuation_date, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(valuation_date,
                                cds_contracts,
                                libor_curve,
                                recovery_rate)

        issuer_curves.append(issuer_curve)

    return issuer_curves

##########################################################################


def test_FinCDSBasket():

    tradeDate = Date(1, 3, 2007)
    step_in_date = tradeDate.addDays(1)
    valuation_date = tradeDate.addDays(1)

    libor_curve = buildIborCurve(tradeDate)

    basketMaturity = Date(20, 12, 2011)

    cdsIndex = CDSIndexPortfolio()

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "====================== INHOMOGENEOUS CURVE ==========================")
    testCases.banner(
        "===================================================================")

    num_credits = 5
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046

    testCases.header("LABELS", "VALUE")

    if 1 == 0:
        issuer_curves = loadHomogeneousSpreadCurves(valuation_date,
                                                   libor_curve,
                                                   spd3Y,
                                                   spd5Y,
                                                   spd7Y,
                                                   spd10Y,
                                                   num_credits)
    else:
        issuer_curves = loadHeterogeneousSpreadCurves(valuation_date, libor_curve)
        issuer_curves = issuer_curves[0:num_credits]

    intrinsicSpd = cdsIndex.intrinsicSpread(valuation_date,
                                            step_in_date,
                                            basketMaturity,
                                            issuer_curves) * 10000.0

    testCases.print("INTRINSIC SPD BASKET MATURITY", intrinsicSpd)

    totalSpd = cdsIndex.totalSpread(valuation_date,
                                    step_in_date,
                                    basketMaturity,
                                    issuer_curves) * 10000.0

    testCases.print("SUMMED UP SPD BASKET MATURITY", totalSpd)

    minSpd = cdsIndex.minSpread(valuation_date,
                                step_in_date,
                                basketMaturity,
                                issuer_curves) * 10000.0

    testCases.print("MINIMUM SPD BASKET MATURITY", minSpd)

    maxSpd = cdsIndex.maxSpread(valuation_date,
                                step_in_date,
                                basketMaturity,
                                issuer_curves) * 10000.0

    testCases.print("MAXIMUM SPD BASKET MATURITY", maxSpd)

    seed = 1967
    basket = CDSBasket(valuation_date,
                       basketMaturity)

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "======================= GAUSSIAN COPULA ===========================")
    testCases.banner(
        "===================================================================")

    testCases.header("TIME", "Trials", "RHO", "NTD", "SPRD", "SPRD_HOMO")

    for ntd in range(1, num_credits + 1):
        for beta in [0.0, 0.5]:
            rho = beta * beta
            betaVector = np.ones(num_credits) * beta
            corrMatrix = corrMatrixGenerator(rho, num_credits)
            for num_trials in [1000]:  # [1000,5000,10000,20000,50000,100000]:
                start = time.time()

                v1 = basket.valueGaussian_MC(valuation_date,
                                             ntd,
                                             issuer_curves,
                                             corrMatrix,
                                             libor_curve,
                                             num_trials,
                                             seed)

                v2 = basket.value1FGaussian_Homo(valuation_date,
                                                 ntd,
                                                 issuer_curves,
                                                 betaVector,
                                                 libor_curve)

                end = time.time()
                period = (end - start)
                testCases.print(
                    period,
                    num_trials,
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
        corrMatrix = corrMatrixGenerator(rho, num_credits)
        for ntd in range(1, num_credits + 1):
            for doF in [3, 6]:
                start = time.time()

                v = basket.valueStudentT_MC(valuation_date,
                                            ntd,
                                            issuer_curves,
                                            corrMatrix,
                                            doF,
                                            libor_curve,
                                            num_trials,
                                            seed)

                end = time.time()
                period = (end - start)
                testCases.print(period, num_trials, rho, doF, ntd, v[2] * 10000)

            start = time.time()
            v = basket.valueGaussian_MC(
                valuation_date,
                ntd,
                issuer_curves,
                corrMatrix,
                libor_curve,
                num_trials,
                seed)
            end = time.time()
            period = (end - start)

            testCases.print(period, num_trials, rho, "GC", ntd, v[2] * 10000)

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
        corrMatrix = corrMatrixGenerator(rho, num_credits)
        for ntd in range(1, num_credits + 1):
            for num_trials in [1000]:
                start = time.time()

                v = basket.valueStudentT_MC(valuation_date,
                                            ntd,
                                            issuer_curves,
                                            corrMatrix,
                                            doF,
                                            libor_curve,
                                            num_trials,
                                            seed)
                end = time.time()
                period = (end - start)
                testCases.print(period, num_trials, rho, ntd, v[2] * 10000)

###############################################################################


def testFinGBMProcess():

    num_assets = 3
    num_paths = 5
    num_time_steps = 1
    t = 1.0
    mus = 0.03 * np.ones(num_assets)
    stock_prices = 100.0 * np.ones(num_assets)
    volatilities = 0.2 * np.ones(num_assets)
    rho = 0.8
    corrMatrix = corrMatrixGenerator(rho, num_assets)
    seed = 1912

    _ = getPathsAssets(num_assets, num_paths, num_time_steps, t,
                       mus, stock_prices, volatilities,
                       corrMatrix, seed)

###############################################################################


testFinGBMProcess()
test_FinCDSBasket()
testCases.compareTestCases()
