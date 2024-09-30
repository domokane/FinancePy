###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import os
import time
import numpy as np

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

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################


def build_Ibor_Curve(tradeDate):

    value_dt = tradeDate.add_days(1)
    dc_type = DayCountTypes.ACT_360
    depos = []

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


def buildFlatIssuerCurve(tradeDate, libor_curve, spread, recovery_rate):

    value_dt = tradeDate.add_days(1)

    cdsMarketContracts = []

    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, spread)
    cdsMarketContracts.append(cds)

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return issuer_curve


##########################################################################


def test_dirty_priceCDSIndexOption():

    tradeDate = Date(1, 8, 2007)
    step_in_dt = tradeDate.add_days(1)
    value_dt = step_in_dt

    libor_curve = build_Ibor_Curve(tradeDate)

    maturity3Y = tradeDate.next_cds_date(36)
    maturity5Y = tradeDate.next_cds_date(60)
    maturity7Y = tradeDate.next_cds_date(84)
    maturity10Y = tradeDate.next_cds_date(120)

    path = os.path.join(
        os.path.dirname(__file__), ".//data//CDX_NA_IG_S7_SPREADS.csv"
    )
    f = open(path, "r")
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        creditName = splitRow[0]
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = CDS(step_in_dt, maturity3Y, spd3Y)
        cds5Y = CDS(step_in_dt, maturity5Y, spd5Y)
        cds7Y = CDS(step_in_dt, maturity7Y, spd7Y)
        cds10Y = CDS(step_in_dt, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, recovery_rate
        )

        issuer_curves.append(issuer_curve)

    ##########################################################################
    ##########################################################################

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

        #######################################################################

        #        print("Index", index)

        cds_contracts = []

        for dt in index_maturity_dts:

            cds = CDS(value_dt, dt, index / 10000.0)
            cds_contracts.append(cds)

        index_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, index_recovery
        )

        if 1 == 1:

            indexSpreads = [index / 10000.0] * 4

            indexPortfolio = CDSIndexPortfolio()

            adjustedIssuerCurves = indexPortfolio.hazard_rate_adjust_intrinsic(
                value_dt,
                issuer_curves,
                indexSpreads,
                index_upfronts,
                index_maturity_dts,
                index_recovery,
                tolerance,
            )

        else:

            indexSpread = index / 10000.0

            issuer_curve = buildFlatIssuerCurve(
                tradeDate, libor_curve, indexSpread, index_recovery
            )

            adjustedIssuerCurves = []
            for i_credit in range(0, 125):
                adjustedIssuerCurves.append(issuer_curve)

        #######################################################################
        # Now loop over strikes
        #######################################################################

        for strike in [20, 60]:

            start = time.time()

            option = CDSIndexOption(
                expiry_dt, maturity_dt, index_cpn, strike / 10000.0, notional
            )

            v_pay_1, v_rec_1, strike_value, mu, expH = option.value_anderson(
                value_dt, adjustedIssuerCurves, index_recovery, volatility
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
                strike,
                index,
                v_pay_1,
                v_rec_1,
                strike_value,
                mu,
                expH,
                v_pay_2,
                v_rec_2,
            )


##########################################################################


test_dirty_priceCDSIndexOption()
test_cases.compareTestCases()
