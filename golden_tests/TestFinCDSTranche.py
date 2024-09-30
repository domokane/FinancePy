###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import time

import sys

sys.path.append("..")

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


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################


def build_Ibor_Curve(tradeDate):

    value_dt = tradeDate.add_days(1)
    dc_type = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    settle_dt = value_dt

    maturity_dt = settle_dt.add_months(12)
    swap1 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(24)
    swap2 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(36)
    swap3 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(48)
    swap4 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(60)
    swap5 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type
    )
    swaps.append(swap5)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)
    return libor_curve


##############################################################################


def loadHomogeneousCDSCurves(
    value_dt,
    libor_curve,
    cdsSpread3Y,
    cdsSpread5Y,
    cdsSpread7Y,
    cdsSpread10Y,
    num_credits,
):

    maturity3Y = value_dt.next_cds_date(36)
    maturity5Y = value_dt.next_cds_date(60)
    maturity7Y = value_dt.next_cds_date(84)
    maturity10Y = value_dt.next_cds_date(120)

    recovery_rate = 0.40

    cds3Y = CDS(value_dt, maturity3Y, cdsSpread3Y)
    cds5Y = CDS(value_dt, maturity5Y, cdsSpread5Y)
    cds7Y = CDS(value_dt, maturity7Y, cdsSpread7Y)
    cds10Y = CDS(value_dt, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuer_curve = CDSCurve(value_dt, contracts, libor_curve, recovery_rate)

    issuer_curves = []
    for _ in range(0, num_credits):
        issuer_curves.append(issuer_curve)

    return issuer_curves


##########################################################################


def loadHeterogeneousSpreadCurves(value_dt, libor_curve):

    maturity3Y = value_dt.next_cds_date(36)
    maturity5Y = value_dt.next_cds_date(60)
    maturity7Y = value_dt.next_cds_date(84)
    maturity10Y = value_dt.next_cds_date(120)
    path = os.path.join(
        os.path.dirname(__file__), ".//data//CDX_NA_IG_S7_SPREADS.csv"
    )
    f = open(path, "r")
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

        cds3Y = CDS(value_dt, maturity3Y, spd3Y)
        cds5Y = CDS(value_dt, maturity5Y, spd5Y)
        cds7Y = CDS(value_dt, maturity7Y, spd7Y)
        cds10Y = CDS(value_dt, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, recovery_rate
        )

        issuer_curves.append(issuer_curve)

    return issuer_curves


##########################################################################


def test_FinCDSTranche():

    tradeDate = Date(1, 3, 2007)
    step_in_dt = tradeDate.add_days(1)
    value_dt = tradeDate.add_days(1)

    test_cases.header("DATE")
    test_cases.print(str((tradeDate)))
    test_cases.print(str((step_in_dt)))
    test_cases.print(str((value_dt)))

    libor_curve = build_Ibor_Curve(tradeDate)

    trancheMaturity = Date(20, 12, 2011)
    tranche1 = CDSTranche(value_dt, trancheMaturity, 0.00, 0.03)
    tranche2 = CDSTranche(value_dt, trancheMaturity, 0.03, 0.06)
    tranche3 = CDSTranche(value_dt, trancheMaturity, 0.06, 0.09)
    tranche4 = CDSTranche(value_dt, trancheMaturity, 0.09, 0.12)
    tranche5 = CDSTranche(value_dt, trancheMaturity, 0.12, 0.22)
    tranche6 = CDSTranche(value_dt, trancheMaturity, 0.22, 0.60)
    tranche7 = CDSTranche(value_dt, trancheMaturity, 0.00, 0.60)
    tranches = [
        tranche1,
        tranche2,
        tranche3,
        tranche4,
        tranche5,
        tranche6,
        tranche7,
    ]

    corr1 = 0.30
    corr2 = 0.35
    upfront = 0.0
    spd = 0.0

    cdsIndex = CDSIndexPortfolio()

    ##########################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "====================== HOMOGENEOUS CURVE =========================="
    )
    test_cases.banner(
        "==================================================================="
    )
    num_credits = 125
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046

    issuer_curves = loadHomogeneousCDSCurves(
        value_dt, libor_curve, spd3Y, spd5Y, spd7Y, spd10Y, num_credits
    )

    intrinsicSpd = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, trancheMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.header("LABEL", "VALUE")
    test_cases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    test_cases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    test_cases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for num_points in [40]:
                start = time.time()
                v = tranche.value_bc(
                    value_dt,
                    issuer_curves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    num_points,
                    method,
                )
                end = time.time()
                period = end - start
                test_cases.print(
                    method,
                    period,
                    num_points,
                    tranche.k1,
                    tranche.k2,
                    v[3] * 10000,
                )

    ##########################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "=================== HETEROGENEOUS CURVES =========================="
    )
    test_cases.banner(
        "==================================================================="
    )

    issuer_curves = loadHeterogeneousSpreadCurves(value_dt, libor_curve)

    intrinsicSpd = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, trancheMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.header("LABEL", "VALUE")
    test_cases.print("INTRINSIC SPD TRANCHE MATURITY", intrinsicSpd)
    adjustedSpd = intrinsicSpd / 0.6
    test_cases.print("ADJUSTED  SPD TRANCHE MATURITY", adjustedSpd)

    test_cases.header("METHOD", "TIME", "NumPoints", "K1", "K2", "Sprd")

    for method in FinLossDistributionBuilder:
        for tranche in tranches:
            for num_points in [40]:
                start = time.time()
                v = tranche.value_bc(
                    value_dt,
                    issuer_curves,
                    upfront,
                    spd,
                    corr1,
                    corr2,
                    num_points,
                    method,
                )
                end = time.time()
                period = end - start
                test_cases.print(
                    method,
                    period,
                    num_points,
                    tranche.k1,
                    tranche.k2,
                    v[3] * 10000,
                )

    test_cases.banner(
        "==================================================================="
    )


##########################################################################


test_FinCDSTranche()
test_cases.compareTestCases()
