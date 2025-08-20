########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from helpers import (
    build_Ibor_Curve,
    loadHeterogeneousSpreadCurves,
    loadHomogeneousCDSCurves,
)
from financepy.utils.date import Date
from financepy.products.credit.cds_tranche import CDSTranche
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.products.credit.cds_tranche import FinLossDistributionBuilder


trade_dt = Date(1, 3, 2007)
step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt.add_days(1)

libor_curve = build_Ibor_Curve(trade_dt)

trancheMaturity = Date(20, 12, 2011)
tranche1 = CDSTranche(value_dt, trancheMaturity, 0.00, 0.03)
tranche2 = CDSTranche(value_dt, trancheMaturity, 0.03, 0.06)
tranche3 = CDSTranche(value_dt, trancheMaturity, 0.06, 0.09)
tranche4 = CDSTranche(value_dt, trancheMaturity, 0.09, 0.12)
tranche5 = CDSTranche(value_dt, trancheMaturity, 0.12, 0.22)
tranche6 = CDSTranche(value_dt, trancheMaturity, 0.22, 0.60)
tranche7 = CDSTranche(value_dt, trancheMaturity, 0.00, 0.60)
tranches = [
    tranche1,
    tranche2,
    tranche3,
    tranche4,
    tranche5,
    tranche6,
    tranche7,
]

corr1 = 0.30
corr2 = 0.35
upfront = 0.0
spd = 0.0

cdsIndex = CDSIndexPortfolio()


def test_homogeneous():

    num_credits = 125
    spd3Y = 0.0012
    spd5Y = 0.0025
    spd7Y = 0.0034
    spd10Y = 0.0046
    num_points = 40

    issuer_curves = loadHomogeneousCDSCurves(
        value_dt, libor_curve, spd3Y, spd5Y, spd7Y, spd10Y, num_credits
    )

    intrinsicSpd = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, trancheMaturity, issuer_curves
        )
        * 10000.0
    )

    assert round(intrinsicSpd, 4) == 23.9767

    method = FinLossDistributionBuilder.RECURSION
    v = tranche1.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 582.3189

    method = FinLossDistributionBuilder.ADJUSTED_BINOMIAL
    v = tranche3.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 29.9804

    method = FinLossDistributionBuilder.GAUSSIAN
    v = tranche5.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 4.5835

    method = FinLossDistributionBuilder.LHP
    v = tranche7.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 39.9617


def test_heterogeneous():
    num_points = 40

    issuer_curves = loadHeterogeneousSpreadCurves(value_dt, libor_curve)

    intrinsicSpd = (
        cdsIndex.intrinsic_spread(
            value_dt, step_in_dt, trancheMaturity, issuer_curves
        )
        * 10000.0
    )

    assert round(intrinsicSpd, 4) == 34.3326

    method = FinLossDistributionBuilder.RECURSION
    v = tranche1.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 868.1327

    method = FinLossDistributionBuilder.ADJUSTED_BINOMIAL
    v = tranche2.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 173.4304

    method = FinLossDistributionBuilder.GAUSSIAN
    v = tranche4.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 16.1785

    method = FinLossDistributionBuilder.LHP
    v = tranche6.value_bc(
        value_dt, issuer_curves, upfront, spd, corr1, corr2, num_points, method
    )
    assert round(v[3] * 10000, 4) == 0.3386
