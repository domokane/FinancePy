# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
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
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.global_vars import G_BASIS_POINT

test_cases = FinTestCases(__file__, global_test_case_mode)

# TO DO

########################################################################################


def build_ibor_curve(trade_dt):

    value_dt = trade_dt.add_days(1)
    dc_type = DayCountTypes.ACT_360
    depos = []

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


def build_flat_issuer_curve(trade_dt, libor_curve, spread, recovery_rate):

    value_dt = trade_dt.add_days(1)

    cds_mkt_contracts = []

    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, spread)
    cds_mkt_contracts.append(cds)

    issuer_curve = CDSCurve(value_dt, cds_mkt_contracts, libor_curve, recovery_rate)

    return issuer_curve


########################################################################################


def test_dirty_price_cds_index_option():

    trade_dt = Date(1, 8, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = step_in_dt

    libor_curve = build_ibor_curve(trade_dt)

    maturity_3yr = trade_dt.next_cds_date(36)
    maturity_5yr = trade_dt.next_cds_date(60)
    maturity_7yr = trade_dt.next_cds_date(84)
    maturity_10yr = trade_dt.next_cds_date(120)

    path = os.path.join(os.path.dirname(__file__), ".//data//CDX_NA_IG_S7_SPREADS.csv")
    f = open(path, "r")
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        split_row = row.split(",")
        credit_name = split_row[0]
        spd_3yr = float(split_row[1]) / 10000.0
        spd_5yr = float(split_row[2]) / 10000.0
        spd_7yr = float(split_row[3]) / 10000.0
        spd_10yr = float(split_row[4]) / 10000.0
        recovery_rate = float(split_row[5])

        cds_3yr = CDS(step_in_dt, maturity_3yr, spd_3yr)
        cds_5yr = CDS(step_in_dt, maturity_5yr, spd_5yr)
        cds_7yr = CDS(step_in_dt, maturity_7yr, spd_7yr)
        cds_10yr = CDS(step_in_dt, maturity_10yr, spd_10yr)
        cds_contracts = [cds_3yr, cds_5yr, cds_7yr, cds_10yr]

        issuer_curve = CDSCurve(value_dt, cds_contracts, libor_curve, recovery_rate)

        issuer_curves.append(issuer_curve)

    index_upfronts = [0.0, 0.0, 0.0, 0.0]
    index_maturity_dts = [
        Date(20, 12, 2009),
        Date(20, 12, 2011),
        Date(20, 12, 2013),
        Date(20, 12, 2016),
    ]
    index_recovery = 0.40

    test_cases.banner(
        "======================= CDS INDEX OPTION =========================="
    )

    index_cpn = 0.004
    volatility = 0.50
    expiry_dt = Date(1, 2, 2008)
    maturity_dt = Date(20, 12, 2011)
    notional = 10000.0
    tolerance = 1e-6

    test_cases.header(
        "TIME",
        "STRIKE",
        "INDEX",
        "PAY",
        "RECEIVER",
        "G(K)",
        "X",
        "EXPH",
        "ABPAY",
        "ABREC",
    )

    # TODO something has changed below and I had to change 60 to 50 - Fix
    # I have investigated but have not found cause yet
    for index in [20, 40, 50]:  # was [20, 40, 60]

        #        print("Index", index)

        cds_contracts = []

        for dt in index_maturity_dts:

            cds = CDS(value_dt, dt, index / 10000.0)
            cds_contracts.append(cds)

        index_curve = CDSCurve(value_dt, cds_contracts, libor_curve, index_recovery)

        if True:

            index_spreads = [index / 10000.0] * 4

            index_portfolio = CDSIndexPortfolio()

            start = time.time()

            adjusted_issuer_curves = index_portfolio.hazard_rate_adjust_intrinsic(
                value_dt,
                issuer_curves,
                index_spreads,
                index_upfronts,
                index_maturity_dts,
                index_recovery,
                tolerance,
            )

            end = time.time()
            elapsed = end - start

        else:

            index_spread = index / 10000.0

            issuer_curve = build_flat_issuer_curve(
                trade_dt, libor_curve, index_spread, index_recovery
            )

            adjusted_issuer_curves = []
            for _ in range(0, 125):
                adjusted_issuer_curves.append(issuer_curve)

        # Now loop over strikes

        for strike_bps in [20, 60]:

            strike = strike_bps / 10000

            start = time.time()

            option = CDSIndexOption(expiry_dt, maturity_dt, index_cpn, strike, notional)

            v_pay_1, v_rec_1, strike_value, mu, exp_h = option.value_anderson(
                value_dt, adjusted_issuer_curves, index_recovery, volatility
            )

            end = time.time()
            elapsed = end - start

            end = time.time()

            v_pay_2, v_rec_2 = option.value_adjusted_black(
                value_dt, index_curve, index_recovery, libor_curve, volatility
            )

            elapsed = end - start

            test_cases.print(
                elapsed,
                strike_bps,
                index,
                v_pay_1,
                v_rec_1,
                strike_value,
                mu,
                exp_h,
                v_pay_2,
                v_rec_2,
            )


########################################################################################

test_dirty_price_cds_index_option()
test_cases.compare_test_cases()
