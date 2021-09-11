###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_tranche import CDSTranche
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.products.credit.cds_tranche import FinLossDistributionBuilder
import os
import time

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################

def build_Ibor_Curve(tradeDate):

    valuation_date = tradeDate.add_days(1)
    dcType = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    settlement_date = valuation_date

    maturity_date = settlement_date.add_months(12)
    swap1 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.add_months(24)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.add_months(36)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.add_months(48)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.add_months(60)
    swap5 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)
    return libor_curve

##############################################################################


def loadHomogeneousCDSCurves(valuation_date,
                             libor_curve,
                             cdsSpread3Y,
                             cdsSpread5Y,
                             cdsSpread7Y,
                             cdsSpread10Y,
                             num_credits):

    maturity3Y = valuation_date.next_cds_date(36)
    maturity5Y = valuation_date.next_cds_date(60)
    maturity7Y = valuation_date.next_cds_date(84)
    maturity10Y = valuation_date.next_cds_date(120)

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

    maturity3Y = valuation_date.next_cds_date(36)
    maturity5Y = valuation_date.next_cds_date(60)
    maturity7Y = valuation_date.next_cds_date(84)
    maturity10Y = valuation_date.next_cds_date(120)
    path = os.path.join(os.path.dirname(__file__),
                        './/data//CDX_NA_IG_S7_SPREADS.csv')
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


def test_FinCDSTranche():

    tradeDate = Date(1, 3, 2007)
    step_in_date = tradeDate.add_days(1)
    valuation_date = tradeDate.add_days(1)

    testCases.header("DATE")
    testCases.print(str((tradeDate)))
    testCases.print(str((step_in_date)))
    testCases.print(str((valuation_date)))

    libor_curve = build_Ibor_Curve(tradeDate)

    trancheMaturity = Date(20, 12, 2011)
    tranche1 = CDSTranche(valuation_date, trancheMaturity, 0.00, 0.03)
    tranche2 = CDSTranche(valuation_date, trancheMaturity, 0.03, 0.06)
    tranche3 = CDSTranche(valuation_date, trancheMaturity, 0.06, 0.09)
    tranche4 = CDSTranche(valuation_date, trancheMaturity, 0.09, 0.12)
    tranche5 = CDSTranche(valuation_date, trancheMaturity, 0.12, 0.22)
    tranche6 = CDSTranche(valuation_date, trancheMaturity, 0.22, 0.60)
    tranche7 = CDSTranche(valuation_date, trancheMaturity, 0.00, 0.60)
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

    cdsIndex = CDSIndexPortfolio()

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner(
        "====================== HOMOGENEOUS CURVE ==========================")
    testCases.banner(
        "===================================================================")
    num_credits = 125
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
                                             num_credits)

    intrinsicSpd = cdsIndex.intrinsic_spread(valuation_date,
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
                v = tranche.value_bc(
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

    intrinsicSpd = cdsIndex.intrinsic_spread(valuation_date,
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
                v = tranche.value_bc(
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
