###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from helpers import build_Ibor_Curve, loadHeterogeneousSpreadCurves
from financepy.utils.date import Date
from financepy.utils.math import corr_matrix_generator
from financepy.products.credit.cds_basket import CDSBasket
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
import numpy as np
from os.path import dirname, join


tradeDate = Date(1, 3, 2007)
step_in_date = tradeDate.add_days(1)
value_date = tradeDate.add_days(1)

libor_curve = build_Ibor_Curve(tradeDate)

basketMaturity = Date(20, 12, 2011)

cdsIndex = CDSIndexPortfolio()

num_credits = 5
issuer_curves = loadHeterogeneousSpreadCurves(
    value_date, libor_curve)
issuer_curves = issuer_curves[0:num_credits]

seed = 1967
basket = CDSBasket(value_date,
                   basketMaturity)


def test_inhomogeneous_curve():
    intrinsicSpd = cdsIndex.intrinsic_spread(value_date,
                                             step_in_date,
                                             basketMaturity,
                                             issuer_curves) * 10000.0
    assert round(intrinsicSpd, 4) == 32.0971

    totalSpd = cdsIndex.total_spread(value_date,
                                     step_in_date,
                                     basketMaturity,
                                     issuer_curves) * 10000.0
    assert round(totalSpd, 4) == 161.3169

    minSpd = cdsIndex.min_spread(value_date,
                                 step_in_date,
                                 basketMaturity,
                                 issuer_curves) * 10000.0
    assert round(minSpd, 4) == 10.6722

    maxSpd = cdsIndex.max_spread(value_date,
                                 step_in_date,
                                 basketMaturity,
                                 issuer_curves) * 10000.0
    assert round(maxSpd, 4) == 81.1466


def test_gaussian_copula():
    num_trials = 1000

    ntd = 1
    beta = 0.0
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v1 = basket.value_gaussian_mc(value_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  libor_curve,
                                  num_trials,
                                  seed)

    v2 = basket.value_1f_gaussian_homo(value_date,
                                       ntd,
                                       issuer_curves,
                                       beta_vector,
                                       libor_curve)

    assert round(v1[2] * 10000, 4) == 151.7163
    assert round(v2[3] * 10000, 4) == 160.0121

    ntd = 2
    beta = 0.5
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v1 = basket.value_gaussian_mc(value_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  libor_curve,
                                  num_trials,
                                  seed)

    v2 = basket.value_1f_gaussian_homo(value_date,
                                       ntd,
                                       issuer_curves,
                                       beta_vector,
                                       libor_curve)

    assert round(v1[2] * 10000, 4) == 15.6402
    assert round(v2[3] * 10000, 4) == 16.6395


def test_student_t():
    num_trials = 1000
    doF = 5

    ntd = 1
    beta = 0.0
    rho = beta * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v = basket.value_student_t_mc(value_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  doF,
                                  libor_curve,
                                  num_trials,
                                  seed)

    assert round(v[2] * 10000, 4) == 133.6284

    ntd = 2
    beta = 0.5
    rho = beta * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v = basket.value_student_t_mc(value_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  doF,
                                  libor_curve,
                                  num_trials,
                                  seed)

    assert round(v[2] * 10000, 4) == 26.5991
