# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from os.path import dirname, join

from financepy.utils.date import Date
from financepy.market.curves.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from .helpers import build_ibor_curve

########################################################################################


def test_cds_index_adjust_spreads():

    trade_dt = Date(1, 8, 2007)
    step_in_dt = trade_dt.add_days(1)
    value_dt = trade_dt.add_days(1)

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

        issuer_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, recovery_rate
        )

        issuer_curves.append(issuer_curve)

    # Now determine the average spread of the index
    cds_index = CDSIndexPortfolio()

    avg_spd_3yr = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity_3yr, issuer_curves
        )
        * 10000.0
    )

    avg_spd_5yr = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity_5yr, issuer_curves
        )
        * 10000.0
    )

    avg_spd_7yr = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity_7yr, issuer_curves
        )
        * 10000.0
    )

    avg_spd_10yr = (
        cds_index.average_spread(
            value_dt, step_in_dt, maturity_10yr, issuer_curves
        )
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
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity_3yr, issuer_curves
        )
        * 10000.0
    )

    intrinsic_spd_5yr = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity_5yr, issuer_curves
        )
        * 10000.0
    )

    intrinsic_spd_7yr = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity_7yr, issuer_curves
        )
        * 10000.0
    )

    intrinsic_spd_10yr = (
        cds_index.intrinsic_spread(
            value_dt, step_in_dt, maturity_10yr, issuer_curves
        )
        * 10000.0
    )

    assert round(intrinsic_spd_3yr, 3) == 19.679
    assert round(intrinsic_spd_5yr, 3) == 35.538
    assert round(intrinsic_spd_7yr, 3) == 49.008
    assert round(intrinsic_spd_10yr, 3) == 61.407


def test_cds_index_spread_adjust_intrinsic_reprices_index():
    """The spread adjustment must move the intrinsic spread of the portfolio
    onto the index spread for each index maturity (it used to be a no-op)."""

    import numpy as np

    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.day_count import DayCountTypes
    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

    value_dt = Date(20, 12, 2027)
    step_in_dt = value_dt.add_days(1)
    libor_curve = FlatDiscountCurve(
        value_dt, 0.03, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    maturity_dts = [Date(20, 12, 2028), Date(20, 12, 2031)]

    rng = np.random.default_rng(1)
    issuer_curves = []
    for spd in rng.uniform(0.005, 0.03, 10):
        cds_contracts = [CDS(value_dt, mat_dt, spd) for mat_dt in maturity_dts]
        issuer_curves.append(CDSCurve(value_dt, cds_contracts, libor_curve, 0.4))

    index_cpns = [0.0100, 0.0150]
    index_upfronts = [0.0, 0.0]
    cds_index = CDSIndexPortfolio()

    adjusted_curves = cds_index.spread_adjust_intrinsic(
        value_dt, issuer_curves, index_cpns, index_upfronts, maturity_dts, 0.4
    )

    for mat_dt, cpn in zip(maturity_dts, index_cpns):
        before = cds_index.intrinsic_spread(value_dt, step_in_dt, mat_dt, issuer_curves)
        after = cds_index.intrinsic_spread(value_dt, step_in_dt, mat_dt, adjusted_curves)
        assert abs(before - cpn) > 0.0020
        assert abs(after - cpn) < 0.0001

    # With a non-zero upfront the adjusted curves must reprice the index value
    index_upfronts = [0.01, -0.005]
    adjusted_curves = cds_index.spread_adjust_intrinsic(
        value_dt, issuer_curves, index_cpns, index_upfronts, maturity_dts, 0.4
    )
    for mat_dt, cpn, upfront in zip(maturity_dts, index_cpns, index_upfronts):
        index_cds = CDS(value_dt, mat_dt, cpn)
        clean_values = [index_cds.value(value_dt, c, 0.4)[1] for c in adjusted_curves]
        assert abs(np.mean(clean_values) / index_cds.notional - upfront) < 1e-6
