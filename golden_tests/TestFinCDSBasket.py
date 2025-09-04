# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import numpy as np

import add_fp_to_path

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


from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

# TO DO

########################################################################################


def build_ibor_curve(trade_dt):

    value_dt = trade_dt.add_days(1)
    dc_type = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    settle_dt = value_dt

    maturity_dt = settle_dt.add_months(12)
    swap1 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(24)
    swap2 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(36)
    swap3 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type)
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(48)
    swap4 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(60)
    swap5 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type)
    swaps.append(swap5)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


########################################################################################


def load_homogeneous_spread_curves(
    value_dt,
    libor_curve,
    cds_spd_3yr,
    cds_spd_5yr,
    cds_spd_7yr,
    cds_spd_10yr,
    num_credits,
):

    maturity_3yr = value_dt.next_cds_date(36)
    maturity_5yr = value_dt.next_cds_date(60)
    maturity_7yr = value_dt.next_cds_date(84)
    maturity_10yr = value_dt.next_cds_date(120)

    recovery_rate = 0.40

    cds_3yr = CDS(value_dt, maturity_3yr, cds_spd_3yr)
    cds_5yr = CDS(value_dt, maturity_5yr, cds_spd_5yr)
    cds_7yr = CDS(value_dt, maturity_7yr, cds_spd_7yr)
    cds_10yr = CDS(value_dt, maturity_10yr, cds_spd_10yr)

    contracts = [cds_3yr, cds_5yr, cds_7yr, cds_10yr]

    issuer_curve = CDSCurve(value_dt, contracts, libor_curve, recovery_rate)

    issuer_curves = []
    for _ in range(0, num_credits):
        issuer_curves.append(issuer_curve)

    return issuer_curves


########################################################################################


def load_hetero_spread_curves(value_dt, libor_curve):

    maturity_3yr = value_dt.next_cds_date(36)
    maturity_5yr = value_dt.next_cds_date(60)
    maturity_7yr = value_dt.next_cds_date(84)
    maturity_10yr = value_dt.next_cds_date(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, "r")

    data = f.readlines()
    issuer_curves = []

    for row in data[1:]:

        split_row = row.split(",")
        spd_3yr = float(split_row[1]) / 10000.0
        spd_5yr = float(split_row[2]) / 10000.0
        spd_7yr = float(split_row[3]) / 10000.0
        spd_10yr = float(split_row[4]) / 10000.0
        recovery_rate = float(split_row[5])

        cds_3yr = CDS(value_dt, maturity_3yr, spd_3yr)
        cds_5yr = CDS(value_dt, maturity_5yr, spd_5yr)
        cds_7yr = CDS(value_dt, maturity_7yr, spd_7yr)
        cds_10yr = CDS(value_dt, maturity_10yr, spd_10yr)
        cds_contracts = [cds_3yr, cds_5yr, cds_7yr, cds_10yr]

        issuer_curve = CDSCurve(value_dt, cds_contracts, libor_curve, recovery_rate)

        issuer_curves.append(issuer_curve)

    return issuer_curves


########################################################################################


def test_fin_cds_basket():

    trade_dt = Date(1, 3, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = trade_dt.add_days(1)

    libor_curve = build_ibor_curve(trade_dt)

    basket_maturity = Date(20, 12, 2011)

    cds_index = CDSIndexPortfolio()

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
    spd_3yr = 0.0012
    spd_5yr = 0.0025
    spd_7yr = 0.0034
    spd_10yr = 0.0046

    test_cases.header("LABELS", "VALUE")

    if 1 == 1:
        issuer_curves = load_homogeneous_spread_curves(
            value_dt, libor_curve, spd_3yr, spd_5yr, spd_7yr, spd_10yr, num_credits
        )
    else:
        issuer_curves = load_hetero_spread_curves(value_dt, libor_curve)
        issuer_curves = issuer_curves[0:num_credits]

    intrinsic_spd = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )

    test_cases.print("INTRINSIC SPD BASKET MATURITY", intrinsic_spd)

    total_spd = (
        cds_index.total_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )

    test_cases.print("SUMMED UP SPD BASKET MATURITY", total_spd)

    min_spd = (
        cds_index.min_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )

    test_cases.print("MINIMUM SPD BASKET MATURITY", min_spd)

    max_spd = (
        cds_index.max_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )

    test_cases.print("MAXIMUM SPD BASKET MATURITY", max_spd)

    seed = 1967
    basket = CDSBasket(value_dt, basket_maturity)

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
            for do_f in [3, 4]:
                start = time.time()

                v = basket.value_student_t_mc(
                    value_dt,
                    ntd,
                    issuer_curves,
                    corr_matrix,
                    do_f,
                    libor_curve,
                    num_trials,
                    seed,
                )

                end = time.time()
                period = end - start
                test_cases.print(period, num_trials, rho, do_f, ntd, v[2] * 10000)

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
    do_f = 5
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
                    do_f,
                    libor_curve,
                    num_trials,
                    seed,
                )
                end = time.time()
                period = end - start
                test_cases.print(period, num_trials, rho, ntd, v[2] * 10000)


########################################################################################

test_fin_cds_basket()
test_cases.compare_test_cases()
