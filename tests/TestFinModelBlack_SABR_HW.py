###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.models.FinModelSABR import FinModelSABR
from financepy.models.FinModelSABRShifted import FinModelSABRShifted

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinSABR():

    strikes = np.linspace(0.01, 0.06, 10)

    alpha = 0.060277
    beta = 0.5
    rho = 0.2097
    nu = 0.75091
    model1 = FinModelSABR(alpha, beta, rho, nu)

    alpha = 0.058484
    beta = 0.5
    rho = 0.20568
    nu = 0.79647
    model2 = FinModelSABR(alpha, beta, rho, nu)

    f = 0.0350
    T = 1.0

    vols1 = model1.blackVol(f, strikes, T)
    vols2 = model2.blackVol(f, strikes, T)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols1)
        plt.plot(strikes, vols2)
        plt.title("SABR")

###############################################################################


def test_FinShiftedSABRSimple():

    strikes = np.linspace(0.01, 0.06, 10)

    alpha = 0.060277
    beta = 0.5
    rho = 0.2097
    nu = 0.75091
    model1 = FinModelSABRShifted(alpha, beta, rho, nu, 0.0)

    alpha = 0.058484
    beta = 0.5
    rho = 0.20568
    nu = 0.79647
    model2 = FinModelSABRShifted(alpha, beta, rho, nu, 0.0)

    f = 0.0350
    T = 1.0

    vols1 = model1.blackVol(f, strikes, T)
    vols2 = model2.blackVol(f, strikes, T)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols1)
        plt.plot(strikes, vols2)
        plt.title("Shifted SIMPLE SABR")

###############################################################################


def test_FinShiftedSABR():

    strikes = np.linspace(-0.006, 0.016, 10)

    alpha = 0.013345
    beta = 0.5
    rho = 0.46698
    nu = 0.49861
    shift = 0.008

    model = FinModelSABRShifted(alpha, beta, rho, nu, shift)

    f = 0.0006384
    T = 1.0

    vols = model.blackVol(f, strikes, T)

    if PLOT_GRAPHS:
        plt.figure()
        plt.plot(strikes, vols)
        plt.title("SHIFTED SABR")

###############################################################################


test_FinSABR()
test_FinShiftedSABRSimple()
test_FinShiftedSABR()

testCases.compareTestCases()
