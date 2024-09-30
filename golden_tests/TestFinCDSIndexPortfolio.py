###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os

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
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

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


##############################################################################


def buildIssuerCurve(tradeDate, libor_curve):

    value_dt = tradeDate.add_days(1)

    cdsMarketContracts = []

    cds_cpn = 0.0048375
    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return issuer_curve


##########################################################################


def test_CDSIndexPortfolio():

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
    # Now determine the average spread of the index
    ##########################################################################

    cds_index = CDSIndexPortfolio()

    averageSpd3Y = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity3Y, issuer_curves
        )
        * 10000.0
    )

    averageSpd5Y = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity5Y, issuer_curves
        )
        * 10000.0
    )

    averageSpd7Y = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity7Y, issuer_curves
        )
        * 10000.0
    )

    averageSpd10Y = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity10Y, issuer_curves
        )
        * 10000.0
    )

    test_cases.header("LABEL", "VALUE")
    test_cases.print("AVERAGE SPD 3Y", averageSpd3Y)
    test_cases.print("AVERAGE SPD 5Y", averageSpd5Y)
    test_cases.print("AVERAGE SPD 7Y", averageSpd7Y)
    test_cases.print("AVERAGE SPD 10Y", averageSpd10Y)

    ##########################################################################
    # Now determine the intrinsic spread of the index to the same maturity
    # dates. As the single name CDS contracts
    ##########################################################################

    cds_index = CDSIndexPortfolio()

    intrinsicSpd3Y = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity3Y, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd5Y = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity5Y, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd7Y = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity7Y, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd10Y = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity10Y, issuer_curves
        )
        * 10000.0
    )

    ##########################################################################
    ##########################################################################

    test_cases.header("LABEL", "VALUE")
    test_cases.print("INTRINSIC SPD 3Y", intrinsicSpd3Y)
    test_cases.print("INTRINSIC SPD 5Y", intrinsicSpd5Y)
    test_cases.print("INTRINSIC SPD 7Y", intrinsicSpd7Y)
    test_cases.print("INTRINSIC SPD 10Y", intrinsicSpd10Y)


test_CDSIndexPortfolio()
test_cases.compareTestCases()
