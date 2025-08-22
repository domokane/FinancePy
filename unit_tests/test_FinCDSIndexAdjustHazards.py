# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .helpers import build_ibor_curve
from financepy.utils.date import Date
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from os.path import dirname, join

########################################################################################


def test_perform_cds_index_hazard_rate_adjustment():

    trade_dt = Date(1, 8, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = step_in_dt

    libor_curve = build_ibor_curve(trade_dt)

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

        split_row = row.split(",")
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

    # Now determine the average spread of the index
    cds_index = CDSIndexPortfolio()

    avg_spd_3yr = (
        cds_index.average_spread(value_dt, step_in_dt, maturity_3yr, issuer_curves)
        * 10000.0
    )

    avg_spd_5yr = (
        cds_index.average_spread(value_dt, step_in_dt, maturity_5yr, issuer_curves)
        * 10000.0
    )

    avg_spd_7yr = (
        cds_index.average_spread(value_dt, step_in_dt, maturity_7yr, issuer_curves)
        * 10000.0
    )

    avg_spd_10yr = (
        cds_index.average_spread(value_dt, step_in_dt, maturity_10yr, issuer_curves)
        * 10000.0
    )

    assert round(avg_spd_3yr, 4) == 19.8221
    assert round(avg_spd_5yr, 4) == 36.0357
    assert round(avg_spd_7yr, 4) == 50.1336
    assert round(avg_spd_10yr, 4) == 63.6622

    # Now determine the intrinsic spread of the index to the same maturity dates
    # As the single name CDS contracts
    cds_index = CDSIndexPortfolio()

    intrinsic_spd_3yr = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_3yr, issuer_curves)
        * 10000.0
    )

    intrinsic_spd_5yr = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_5yr, issuer_curves)
        * 10000.0
    )

    intrinsic_spd_7yr = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_7yr, issuer_curves)
        * 10000.0
    )

    intrinsic_spd_10yr = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_10yr, issuer_curves)
        * 10000.0
    )

    assert round(intrinsic_spd_3yr, 4) == 19.6789
    assert round(intrinsic_spd_5yr, 4) == 35.5394
    assert round(intrinsic_spd_7yr, 4) == 49.0121
    assert round(intrinsic_spd_10yr, 4) == 61.4140
