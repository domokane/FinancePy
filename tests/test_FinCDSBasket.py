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
valuation_date = tradeDate.add_days(1)

libor_curve = build_Ibor_Curve(tradeDate)

basketMaturity = Date(20, 12, 2011)

cdsIndex = CDSIndexPortfolio()

num_credits = 5
issuer_curves = loadHeterogeneousSpreadCurves(
    valuation_date, libor_curve)
issuer_curves = issuer_curves[0:num_credits]

seed = 1967
basket = CDSBasket(valuation_date,
                   basketMaturity)


def test_inhomogeneous_curve():
    intrinsicSpd = cdsIndex.intrinsic_spread(valuation_date,
                                             step_in_date,
                                             basketMaturity,
                                             issuer_curves) * 10000.0
    assert round(intrinsicSpd, 4) == 32.0971

    totalSpd = cdsIndex.total_spread(valuation_date,
                                     step_in_date,
                                     basketMaturity,
                                     issuer_curves) * 10000.0
    assert round(totalSpd, 4) == 161.3219

    minSpd = cdsIndex.min_spread(valuation_date,
                                 step_in_date,
                                 basketMaturity,
                                 issuer_curves) * 10000.0
    assert round(minSpd, 4) == 10.6725

    maxSpd = cdsIndex.max_spread(valuation_date,
                                 step_in_date,
                                 basketMaturity,
                                 issuer_curves) * 10000.0
    assert round(maxSpd, 4) == 81.1492


def test_gaussian_copula():
    num_trials = 1000

    ntd = 1
    beta = 0.0
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v1 = basket.value_gaussian_mc(valuation_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  libor_curve,
                                  num_trials,
                                  seed)

    v2 = basket.value_1f_gaussian_homo(valuation_date,
                                       ntd,
                                       issuer_curves,
                                       beta_vector,
                                       libor_curve)

    assert round(v1[2] * 10000, 4) == 151.7163
    assert round(v2[3] * 10000, 4) == 159.021

    ntd = 2
    beta = 0.5
    rho = beta * beta
    beta_vector = np.ones(num_credits) * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v1 = basket.value_gaussian_mc(valuation_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  libor_curve,
                                  num_trials,
                                  seed)

    v2 = basket.value_1f_gaussian_homo(valuation_date,
                                       ntd,
                                       issuer_curves,
                                       beta_vector,
                                       libor_curve)

    assert round(v1[2] * 10000, 4) == 15.4566
    assert round(v2[3] * 10000, 4) == 16.4308


def test_student_t():
    num_trials = 1000
    doF = 5

    ntd = 1
    beta = 0.0
    rho = beta * beta
    corr_matrix = corr_matrix_generator(rho, num_credits)

    v = basket.value_student_t_mc(valuation_date,
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

    v = basket.value_student_t_mc(valuation_date,
                                  ntd,
                                  issuer_curves,
                                  corr_matrix,
                                  doF,
                                  libor_curve,
                                  num_trials,
                                  seed)

    assert round(v[2] * 10000, 4) == 26.2872
