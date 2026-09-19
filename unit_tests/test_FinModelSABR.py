# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.models.sabr import SABR
from financepy.models.sabr import vol_function_sabr

########################################################################################


def test_sabr():

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


########################################################################################


def test_sabr__calibration():

    alpha = 0.28
    beta = 0.5
    rho = -0.09
    nu = 0.1

    strike_vol = 0.1

    f = 0.043
    k = 0.050
    r = 0.03
    t_exp = 2.0

    call_option_type = OptionTypes.EUROPEAN_CALL
    put_option_type = OptionTypes.EUROPEAN_PUT

    df = np.exp(-r * t_exp)

    # Make SABR equivalent to lognormal (Black) model
    # (i.e. alpha = 0, beta = 1, rho = 0, nu = 0, shift = 0)
    model_sabr_01 = SABR(0.0, 1.0, 0.0, 0.0)
    model_sabr_01.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

    implied_lognormal_vol = model_sabr_01.black_vol(f, k, t_exp)
    implied_atm_lognormal_vol = model_sabr_01.black_vol(k, k, t_exp)
    implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol

    assert implied_lognormal_smile == 0.0, "In lognormal model, smile should be flat"
    calibration_error = round(strike_vol - implied_lognormal_vol, 6)
    assert calibration_error == 0.0

    # Volatility: pure SABR dynamics
    model_sabr_02 = SABR(alpha, beta, rho, nu)
    model_sabr_02.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

    implied_lognormal_vol = model_sabr_02.black_vol(f, k, t_exp)
    implied_atm_lognormal_vol = model_sabr_02.black_vol(k, k, t_exp)
    implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol
    calibration_error = round(strike_vol - implied_lognormal_vol, 6)
    assert calibration_error == 0.0

    # Valuation: pure SABR dynamics
    value_call = model_sabr_02.value(f, k, t_exp, df, call_option_type)
    value_put = model_sabr_02.value(f, k, t_exp, df, put_option_type)
    assert round(value_call - value_put, 12) == round(
        df * (f - k), 12
    ), "The method called 'value()' doesn't comply with Call-Put parity"
