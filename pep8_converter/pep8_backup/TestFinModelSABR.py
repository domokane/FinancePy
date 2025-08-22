########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
from financepy.models.sabr import vol_function_sabr
from financepy.models.sabr import SABR
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, global_test_case_mode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################

alpha = 0.28
beta = 1.0
rho = -0.09
nu = 0.21

f = 0.043
k = 0.050
t = 2.0

########################################################################################


def test_SABR():

    test_cases.header("ALPHA", "BETA", "RHO", "VOL")

    for alpha in [0.1, 0.2, 0.3]:
        for beta in [0.5, 1.0, 2.0]:
            for rho in [-0.8, 0.0, 0.8]:
                params = np.array([alpha, beta, rho, nu])
                vol = vol_function_sabr(params, f, k, t)
                test_cases.print(alpha, beta, rho, vol)


########################################################################################


def test_SABR_Calibration():

    beta = 0.5
    rho = -0.09
    nu = 0.1

    strikeVol = 0.1

    f = 0.043
    k = 0.050
    r = 0.03
    t_exp = 2.0

    call_optionType = OptionTypes.EUROPEAN_CALL
    put_optionType = OptionTypes.EUROPEAN_PUT

    df = np.exp(-r * t_exp)

    test_cases.header("TEST", "CALIBRATION ERROR")

    # Make SABR equivalent to lognormal (Black) model
    # (i.e. alpha = 0, beta = 1, rho = 0, nu = 0, shift = 0)
    modelSABR_01 = SABR(0.0, 1.0, 0.0, 0.0)
    modelSABR_01.set_alpha_from_black_vol(strikeVol, f, k, t_exp)

    impliedLognormalVol = modelSABR_01.black_vol(f, k, t_exp)
    impliedATMLognormalVol = modelSABR_01.black_vol(k, k, t_exp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol

    assert impliedLognormalSmile == 0.0, "In lognormal model, smile should be flat"
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    test_cases.print("LOGNORMAL CASE", calibrationError)

    # Volatility: pure SABR dynamics
    modelSABR_02 = SABR(alpha, beta, rho, nu)
    modelSABR_02.set_alpha_from_black_vol(strikeVol, f, k, t_exp)

    impliedLognormalVol = modelSABR_02.black_vol(f, k, t_exp)
    impliedATMLognormalVol = modelSABR_02.black_vol(k, k, t_exp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    test_cases.print("SABR CASE", calibrationError)

    # Valuation: pure SABR dynamics
    valueCall = modelSABR_02.value(f, k, t_exp, df, call_optionType)
    valuePut = modelSABR_02.value(f, k, t_exp, df, put_optionType)
    assert round(valueCall - valuePut, 12) == round(
        df * (f - k), 12
    ), "The method called 'value()' doesn't comply with Call-Put parity"


########################################################################################


test_SABR()
test_SABR_Calibration()

test_cases.compare_test_cases()
