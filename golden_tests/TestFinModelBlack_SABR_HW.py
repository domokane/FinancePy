# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import numpy as np
from financepy.models.sabr import SABR
from financepy.models.sabr_shifted import SABRShifted
from FinTestCases import FinTestCases, global_test_case_mode

import matplotlib.pyplot as plt

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

########################################################################################


def test_fin_sabr():

    strikes = np.linspace(0.01, 0.06, 10)

    alpha = 0.060277
    beta = 0.5
    rho = 0.2097
    nu = 0.75091
    model1 = SABR(alpha, beta, rho, nu)

    alpha = 0.058484
    beta = 0.5
    rho = 0.20568
    nu = 0.79647
    model2 = SABR(alpha, beta, rho, nu)

    f = 0.0350
    t = 1.0

    vols1 = model1.black_vol(f, strikes, t)
    vols2 = model2.black_vol(f, strikes, t)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols1)
        plt.plot(strikes, vols2)
        plt.title("SABR")


########################################################################################


def test_fin_shifted_sabr_simple():

    strikes = np.linspace(0.01, 0.06, 10)

    alpha = 0.060277
    beta = 0.5
    rho = 0.2097
    nu = 0.75091
    model1 = SABRShifted(alpha, beta, rho, nu, 0.0)

    alpha = 0.058484
    beta = 0.5
    rho = 0.20568
    nu = 0.79647
    model2 = SABRShifted(alpha, beta, rho, nu, 0.0)

    f = 0.0350
    t = 1.0

    vols1 = model1.black_vol(f, strikes, t)
    vols2 = model2.black_vol(f, strikes, t)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols1)
        plt.plot(strikes, vols2)
        plt.title("Shifted SIMPLE SABR")


########################################################################################


def test_fin_shifted_sabr():

    strikes = np.linspace(-0.006, 0.016, 10)

    alpha = 0.013345
    beta = 0.5
    rho = 0.46698
    nu = 0.49861
    shift = 0.008

    model = SABRShifted(alpha, beta, rho, nu, shift)

    f = 0.0006384
    t = 1.0

    vols = model.black_vol(f, strikes, t)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols)
        plt.title("SHIFTED SABR")


########################################################################################

test_fin_sabr()
test_fin_shifted_sabr_simple()
test_fin_shifted_sabr()

test_cases.compare_test_cases()
