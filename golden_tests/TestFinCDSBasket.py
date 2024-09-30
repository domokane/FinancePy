###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from os.path import dirname, join

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
from financepy.products.credit.cds_basket import CDSBasket
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.utils.math import corr_matrix_generator


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
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


##########################################################################


def loadHomogeneousSpreadCurves(
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

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, "r")

    data = f.readlines()
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


def test_FinCDSBasket():

    tradeDate = Date(1, 3, 2007)
    step_in_dt = tradeDate.add_days(1)
    value_dt = tradeDate.add_days(1)

    libor_curve = build_Ibor_Curve(tradeDate)

    basketMaturity = Date(20, 12, 2011)

    cdsIndex = CDSIndexPortfolio()

    ##########################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "====================== INHOMOGENEOUS CURVE =========================="
    )
    test_cases.banner(
        "==================================================================="
    )

    num_credits = 5
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046

    test_cases.header("LABELS", "VALUE")

    if 1 == 1:
        issuer_curves = loadHomogeneousSpreadCurves(
            value_dt, libor_curve, spd3Y, spd5Y, spd7Y, spd10Y, num_credits
        )
    else:
        issuer_curves = loadHeterogeneousSpreadCurves(value_dt, libor_curve)
        issuer_curves = issuer_curves[0:num_credits]

    intrinsicSpd = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, basketMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.print("INTRINSIC SPD BASKET MATURITY", intrinsicSpd)

    totalSpd = (
        cdsIndex.total_spread(
            value_dt, step_in_dt, basketMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.print("SUMMED UP SPD BASKET MATURITY", totalSpd)

    minSpd = (
        cdsIndex.min_spread(
            value_dt, step_in_dt, basketMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.print("MINIMUM SPD BASKET MATURITY", minSpd)

    maxSpd = (
        cdsIndex.max_spread(
            value_dt, step_in_dt, basketMaturity, issuer_curves
        )
        * 10000.0
    )

    test_cases.print("MAXIMUM SPD BASKET MATURITY", maxSpd)

    seed = 1967
    basket = CDSBasket(value_dt, basketMaturity)

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "======================= GAUSSIAN COPULA ==========================="
    )
    test_cases.banner(
        "==================================================================="
    )

    test_cases.header("TIME", "Trials", "RHO", "NTD", "SPRD", "SPRD_HOMO")

    for ntd in range(1, num_credits + 1):
        for beta in [0.0, 0.5]:
            rho = beta * beta
            beta_vector = np.ones(num_credits) * beta
            corr_matrix = corr_matrix_generator(rho, num_credits)
            for num_trials in [1000]:  # [1000,5000,10000,20000,50000,100000]:
                start = time.time()

                v1 = basket.value_gaussian_mc(
                    value_dt,
                    ntd,
                    issuer_curves,
                    corr_matrix,
                    libor_curve,
                    num_trials,
                    seed,
                )

                v2 = basket.value_1f_gaussian_homo(
                    value_dt, ntd, issuer_curves, beta_vector, libor_curve
                )

                end = time.time()
                period = end - start
                test_cases.print(
                    period, num_trials, rho, ntd, v1[2] * 10000, v2[3] * 10000
                )

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "==================== STUDENT'S-T CONVERGENCE ======================"
    )
    test_cases.banner(
        "==================================================================="
    )

    test_cases.header("TIME", "TRIALS", "RHO", "DOF", "NTD", "SPRD")

    for beta in [0.0, 0.5]:
        rho = beta**2
        corr_matrix = corr_matrix_generator(rho, num_credits)
        for ntd in range(1, num_credits + 1):
            for doF in [3, 4]:
                start = time.time()

                v = basket.value_student_t_mc(
                    value_dt,
                    ntd,
                    issuer_curves,
                    corr_matrix,
                    doF,
                    libor_curve,
                    num_trials,
                    seed,
                )

                end = time.time()
                period = end - start
                test_cases.print(
                    period, num_trials, rho, doF, ntd, v[2] * 10000
                )

            start = time.time()
            v = basket.value_gaussian_mc(
                value_dt,
                ntd,
                issuer_curves,
                corr_matrix,
                libor_curve,
                num_trials,
                seed,
            )
            end = time.time()
            period = end - start

            test_cases.print(period, num_trials, rho, "GC", ntd, v[2] * 10000)

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner(
        "=================== STUDENT'S T WITH DOF = 5 ======================"
    )
    test_cases.banner(
        "==================================================================="
    )
    doF = 5
    test_cases.header("TIME", "NUMTRIALS", "RHO", "NTD", "SPD")
    for beta in [0.0, 0.5]:
        rho = beta**2
        corr_matrix = corr_matrix_generator(rho, num_credits)
        for ntd in range(1, num_credits + 1):
            for num_trials in [1000]:
                start = time.time()

                v = basket.value_student_t_mc(
                    value_dt,
                    ntd,
                    issuer_curves,
                    corr_matrix,
                    doF,
                    libor_curve,
                    num_trials,
                    seed,
                )
                end = time.time()
                period = end - start
                test_cases.print(period, num_trials, rho, ntd, v[2] * 10000)


###############################################################################

test_FinCDSBasket()
test_cases.compareTestCases()
