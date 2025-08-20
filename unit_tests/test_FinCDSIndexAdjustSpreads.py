########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .helpers import build_Ibor_Curve
from financepy.utils.date import Date
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from os.path import dirname, join


def test_CDSIndexAdjustSpreads():

    trade_dt = Date(1, 8, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = trade_dt.add_days(1)

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

    # Now determine the average spread of the index
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

    assert round(averageSpd3Y, 4) == 19.8221
    assert round(averageSpd5Y, 4) == 36.0357
    assert round(averageSpd7Y, 4) == 50.1336
    assert round(averageSpd10Y, 4) == 63.6622

    # Now determine the intrinsic spread of the index to the same maturity dates
    # As the single name CDS contracts
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

    assert round(intrinsicSpd3Y, 4) == 19.6789
    assert round(intrinsicSpd5Y, 4) == 35.5394
    assert round(intrinsicSpd7Y, 4) == 49.0121
    assert round(intrinsicSpd10Y, 4) == 61.4140
