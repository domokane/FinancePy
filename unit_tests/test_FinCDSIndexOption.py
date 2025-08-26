# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os

from financepy.utils.date import Date
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from .helpers import build_ibor_curve

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

    index_cpn = 0.004
    volatility = 0.50
    expiry_dt = Date(1, 2, 2008)
    maturity_dt = Date(20, 12, 2011)
    notional = 10000.0
    tolerance = 1e-6

    index_strike_results = [
        (20.0, 20.0, [16.0, 6.2, -70.8, 22.9, -60.8, 16.1, 6.1]),
        (25.0, 30.0, [11.8, 16.9, -35.3, 28.6, -40.5, 11.8, 16.8]),
        (50.0, 40.0, [63.3, 4.7, 0.0, 57.3, 60.2, 63.2, 4.7]),
    ]

    for index, strike, results in index_strike_results:

        cds_contracts = []
        for dt in index_maturity_dts:
            cds = CDS(value_dt, dt, index / 10000.0)
            cds_contracts.append(cds)

        index_curve = CDSCurve(value_dt, cds_contracts, libor_curve, index_recovery)

        index_spreads = [index / 10000.0] * 4

        index_portfolio = CDSIndexPortfolio()
        adjusted_issuer_curves = index_portfolio.hazard_rate_adjust_intrinsic(
            value_dt,
            issuer_curves,
            index_spreads,
            index_upfronts,
            index_maturity_dts,
            index_recovery,
            tolerance,
        )

        option = CDSIndexOption(
            expiry_dt, maturity_dt, index_cpn, strike / 10000.0, notional
        )

        v_pay_1, v_rec_1, strike_value, mu, exp_h = option.value_anderson(
            value_dt, adjusted_issuer_curves, index_recovery, volatility
        )

        v_pay_2, v_rec_2 = option.value_adjusted_black(
            value_dt, index_curve, index_recovery, libor_curve, volatility
        )

        assert round(v_pay_1, 1) == results[0]
        assert round(v_rec_1, 1) == results[1]
        assert round(strike_value, 1) == results[2]
        assert round(mu, 1) == results[3]
        assert round(exp_h, 1) == results[4]
        assert round(v_pay_2, 1) == results[5]
        assert round(v_rec_2, 1) == results[6]
