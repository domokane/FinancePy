# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.date import Date
from financepy.utils.math import corr_matrix_generator
from financepy.products.credit.cds_basket import CDSBasket
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from .helpers import build_ibor_curve
from .helpers import load_hetero_spread_curves

trade_dt = Date(1, 3, 2007)
step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt.add_days(1)

libor_curve = build_ibor_curve(trade_dt)

basket_maturity = Date(20, 12, 2011)

cds_index = CDSIndexPortfolio()

num_credits = 5
issuer_curves = load_hetero_spread_curves(value_dt, libor_curve)
issuer_curves = issuer_curves[0:num_credits]

seed = 1967
basket = CDSBasket(value_dt, basket_maturity)

########################################################################################


def test_inhomogeneous_curve():

    intrinsic_spd = (
        cds_index.intrinsic_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )
    assert round(intrinsic_spd, 4) == 32.0963

    total_spd = (
        cds_index.total_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )
    assert round(total_spd, 4) == 161.3163

    min_spd = (
        cds_index.min_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )
    assert round(min_spd, 4) == 10.6722

    max_spd = (
        cds_index.max_spread(value_dt, step_in_dt, basket_maturity, issuer_curves)
        * 10000.0
    )
    assert round(max_spd, 4) == 81.1462


########################################################################################


def test_gaussian_copula():

    num_trials = 1000

    ntd = 1
    beta = 0.0
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

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

    assert round(v1[2] * 10000, 3) == 152.398
    assert round(v2[3] * 10000, 4) == 160.0294

    ntd = 2
    beta = 0.5
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

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

    assert round(v1[2] * 10000, 4) == 15.6482
    assert round(v2[3] * 10000, 4) == 16.6387


########################################################################################


def test_student_t():

    num_trials = 1000
    do_f = 5

    ntd = 1
    beta = 0.0
    rho = beta * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

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

    assert round(v[2] * 10000, 4) == 134.1416

    ntd = 2
    beta = 0.5
    rho = beta * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

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

    assert round(v[2] * 10000, 4) == 26.6214
