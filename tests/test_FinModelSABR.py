###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import OptionTypes
from financepy.models.sabr import SABR
from financepy.models.sabr import vol_function_sabr
import numpy as np


def test_SABR():
    nu = 0.21
    f = 0.043
    k = 0.050
    t = 2.0

    alpha = 0.2
    beta = 0.5
    rho = -0.8
    params = np.array([alpha, beta, rho, nu])
    vol = vol_function_sabr(params, f, k, t)
    assert round(vol, 4) == 0.8971

    alpha = 0.3
    beta = 1.0
    rho = 0.0
    params = np.array([alpha, beta, rho, nu])
    vol = vol_function_sabr(params, f, k, t)
    assert round(vol, 4) == 0.3028

    alpha = 0.1
    beta = 2.0
    rho = 0.8
    params = np.array([alpha, beta, rho, nu])
    vol = vol_function_sabr(params, f, k, t)
    assert round(vol, 4) == 0.0148


def test_SABR_Calibration():
    alpha = 0.28
    beta = 0.5
    rho = -0.09
    nu = 0.1

    strikeVol = 0.1

    f = 0.043
    k = 0.050
    r = 0.03
    texp = 2.0

    call_optionType = OptionTypes.EUROPEAN_CALL
    put_optionType = OptionTypes.EUROPEAN_PUT

    df = np.exp(-r * texp)

    # Make SABR equivalent to lognormal (Black) model
    # (i.e. alpha = 0, beta = 1, rho = 0, nu = 0, shift = 0)
    modelSABR_01 = SABR(0.0, 1.0, 0.0, 0.0)
    modelSABR_01.set_alpha_from_black_vol(strikeVol, f, k, texp)

    impliedLognormalVol = modelSABR_01.black_vol(f, k, texp)
    impliedATMLognormalVol = modelSABR_01.black_vol(k, k, texp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol

    assert impliedLognormalSmile == 0.0, "In lognormal model, smile should be flat"
    calibrationError = round(strikeVol - impliedLognormalVol, 6)
    assert calibrationError == 0.0

    # Volatility: pure SABR dynamics
    modelSABR_02 = SABR(alpha, beta, rho, nu)
    modelSABR_02.set_alpha_from_black_vol(strikeVol, f, k, texp)

    impliedLognormalVol = modelSABR_02.black_vol(f, k, texp)
    impliedATMLognormalVol = modelSABR_02.black_vol(k, k, texp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    calibrationError = round(strikeVol - impliedLognormalVol, 6)
    assert calibrationError == 0.0

    # Valuation: pure SABR dynamics
    valueCall = modelSABR_02.value(f, k, texp, df, call_optionType)
    valuePut = modelSABR_02.value(f, k, texp, df, put_optionType)
    assert round(valueCall - valuePut, 12) == round(df*(f - k), 12), \
        "The method called 'value()' doesn't comply with Call-Put parity"
