########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
from financepy.models.sabr_shifted import SABRShifted
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

########################################################################################

########################################################################################


def test_ShiftedSABR():

    test_cases.header("TEST", "CALIBRATION ERROR")

    alpha = 0.0
    beta = 0.5
    rho = -0.09
    nu = 0.1
    shift = 0.02

    strikeVol = 0.1

    f = 0.043
    k = 0.050
    r = 0.03
    t_exp = 2.0

    call_optionType = OptionTypes.EUROPEAN_CALL
    put_optionType = OptionTypes.EUROPEAN_PUT

    df = np.exp(-r * t_exp)

    # SABR equivalent to lognormal (Black) model (i.e. beta = 1, rho = 0, nu = 0, shift = 0)
    modelSABR_01 = SABRShifted(0.0, 1.0, 0.0, 0.0, 0.0)
    modelSABR_01.set_alpha_from_black_vol(strikeVol, f, k, t_exp)

    impliedLognormalVol = modelSABR_01.black_vol(f, k, t_exp)
    impliedATMLognormalVol = modelSABR_01.black_vol(k, k, t_exp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol

    assert (
        impliedLognormalSmile == 0.0
    ), "In lognormal model, smile should be flat"
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    test_cases.print("LOGNORMAL CASE", calibrationError)

    # Volatility: pure SABR dynamics
    modelSABR_02 = SABRShifted(alpha, beta, rho, nu, shift)
    modelSABR_02.set_alpha_from_black_vol(strikeVol, f, k, t_exp)

    impliedLognormalVol = modelSABR_02.black_vol(f, k, t_exp)
    impliedATMLognormalVol = modelSABR_02.black_vol(k, k, t_exp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    test_cases.print("SABR CASE", calibrationError)

    # Valuation: pure SABR dynamics
    valueCall = modelSABR_02.value(f, k, t_exp, df, call_optionType)
    valuePut = modelSABR_02.value(f, k, t_exp, df, put_optionType)
    assert round((valueCall - valuePut), 12) == round(
        df * (f - k), 12
    ), "The method called 'value()' doesn't comply with Call-Put parity"

    # TODO: adding Call-Put parity test for all sensitivities


########################################################################################


test_ShiftedSABR()
test_cases.compareTestCases()
