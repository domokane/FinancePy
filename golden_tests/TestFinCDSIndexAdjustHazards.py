########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

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
from os.path import dirname, join

test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


def build_Ibor_Curve(trade_dt):

    value_dt = trade_dt.add_days(1)
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


def buildIssuerCurve(trade_dt, libor_curve):

    value_dt = trade_dt.add_days(1)

    cds_mkt_contracts = []

    cds_cpn = 0.0048375
    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cds_mkt_contracts, libor_curve, recovery_rate
    )

    return issuer_curve


##########################################################################


def test_performCDSIndexHazardRateAdjustment():

    trade_dt = Date(1, 8, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = step_in_dt

    libor_curve = build_Ibor_Curve(trade_dt)

    maturity_3yr = trade_dt.next_cds_date(36)
    maturity_5yr = trade_dt.next_cds_date(60)
    maturity_7yr = trade_dt.next_cds_date(84)
    maturity_10yr = trade_dt.next_cds_date(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, "r")

    data = f.readlines()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd_3yr = float(splitRow[1]) / 10000.0
        spd_5yr = float(splitRow[2]) / 10000.0
        spd_7yr = float(splitRow[3]) / 10000.0
        spd_10yr = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds_3yr = CDS(step_in_dt, maturity_3yr, spd_3yr)
        cds_5yr = CDS(step_in_dt, maturity_5yr, spd_5yr)
        cds_7yr = CDS(step_in_dt, maturity_7yr, spd_7yr)
        cds_10yr = CDS(step_in_dt, maturity_10yr, spd_10yr)
        cds_contracts = [cds_3yr, cds_5yr, cds_7yr, cds_10yr]

        issuer_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, recovery_rate
        )

        issuer_curves.append(issuer_curve)

    ##########################################################################
    # Now determine the average spread of the index
    ##########################################################################

    cdsIndex = CDSIndexPortfolio()

    averageSpd3Y = (
        cdsIndex.average_spread(
            value_dt, step_in_dt, maturity_3yr, issuer_curves
        )
        * 10000.0
    )

    averageSpd5Y = (
        cdsIndex.average_spread(
            value_dt, step_in_dt, maturity_5yr, issuer_curves
        )
        * 10000.0
    )

    averageSpd7Y = (
        cdsIndex.average_spread(
            value_dt, step_in_dt, maturity_7yr, issuer_curves
        )
        * 10000.0
    )

    averageSpd10Y = (
        cdsIndex.average_spread(
            value_dt, step_in_dt, maturity_10yr, issuer_curves
        )
        * 10000.0
    )

    test_cases.header("LABEL", "VALUE")
    test_cases.print("AVERAGE SPD 3Y", averageSpd3Y)
    test_cases.print("AVERAGE SPD 5Y", averageSpd5Y)
    test_cases.print("AVERAGE SPD 7Y", averageSpd7Y)
    test_cases.print("AVERAGE SPD 10Y", averageSpd10Y)
    test_cases.banner(
        "==================================================================="
    )

    ##########################################################################
    # Now determine the intrinsic spread of the index to same maturity dates
    # As the single name CDS contracts
    ##########################################################################

    cdsIndex = CDSIndexPortfolio()

    intrinsicSpd3Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, maturity_3yr, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd5Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, maturity_5yr, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd7Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, maturity_7yr, issuer_curves
        )
        * 10000.0
    )

    intrinsicSpd10Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, maturity_10yr, issuer_curves
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
    test_cases.banner(
        "==================================================================="
    )

    ##########################################################################
    ##########################################################################

    index_cpns = [0.002, 0.0037, 0.0050, 0.0063]
    index_upfronts = [0.0, 0.0, 0.0, 0.0]

    index_maturity_dts = [
        Date(20, 12, 2009),
        Date(20, 12, 2011),
        Date(20, 12, 2013),
        Date(20, 12, 2016),
    ]

    index_recoveryRate = 0.40

    tolerance = 1e-4  # should be smaller

    import time

    start = time.time()

    indexPortfolio = CDSIndexPortfolio()
    adjustedIssuerCurves = indexPortfolio.hazard_rate_adjust_intrinsic(
        value_dt,
        issuer_curves,
        index_cpns,
        index_upfronts,
        index_maturity_dts,
        index_recoveryRate,
        tolerance,
    )

    end = time.time()
    test_cases.header("TIME")
    test_cases.print(end - start)

    #    num_credits = len(issuer_curves)
    #    test_cases.print("#","MATURITY","CDS_UNADJ","CDS_ADJ")
    #    for m in range(0,num_credits):
    #        for cds in cds_contracts:
    #            unadjustedSpread = cds.par_spread(value_dt,issuer_curves[m])
    #            adjustedSpread = cds.par_spread(value_dt,adjustedIssuerCurves[m])
    #            test_cases.print(m,str(cds.maturity_dt),"%10.3f"%(unadjustedSpread*10000),"%10.3f" %(adjustedSpread*10000))

    cdsIndex = CDSIndexPortfolio()

    intrinsicSpd3Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, index_maturity_dts[0], adjustedIssuerCurves
        )
        * 10000.0
    )

    intrinsicSpd5Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, index_maturity_dts[1], adjustedIssuerCurves
        )
        * 10000.0
    )

    intrinsicSpd7Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, index_maturity_dts[2], adjustedIssuerCurves
        )
        * 10000.0
    )

    intrinsicSpd10Y = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, index_maturity_dts[3], adjustedIssuerCurves
        )
        * 10000.0
    )

    # If the adjustment works then this should equal the index spreads
    test_cases.header("LABEL", "VALUE")
    test_cases.print("ADJUSTED INTRINSIC SPD 3Y", intrinsicSpd3Y)
    test_cases.print("ADJUSTED INTRINSIC SPD 5Y", intrinsicSpd5Y)
    test_cases.print("ADJUSTED INTRINSIC SPD 7Y", intrinsicSpd7Y)
    test_cases.print("ADJUSTED INTRINSIC SPD 10Y", intrinsicSpd10Y)


########################################################################################


test_performCDSIndexHazardRateAdjustment()
test_cases.compareTestCases()
