# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

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
    assert round(vol, 4) == 0.8969

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


########################################################################################


def _hagan_black_vol(alpha, beta, rho, nu, f, k, t):
    """Hagan et al. (2002) equation (2.17a), written out independently."""
    one_minus_beta = 1.0 - beta
    fk = f * k
    log_fk = np.log(f / k)
    denom = fk ** (one_minus_beta / 2.0) * (
        1.0 + one_minus_beta**2 / 24.0 * log_fk**2 + one_minus_beta**4 / 1920.0 * log_fk**4
    )
    z = nu / alpha * fk ** (one_minus_beta / 2.0) * log_fk
    x = np.log((np.sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho))
    z_over_x = 1.0 if abs(z) < 1e-12 else z / x
    correction = 1.0 + (
        one_minus_beta**2 / 24.0 * alpha**2 / fk**one_minus_beta
        + 0.25 * rho * beta * nu * alpha / fk ** (one_minus_beta / 2.0)
        + (2.0 - 3.0 * rho**2) / 24.0 * nu**2
    ) * t
    return alpha / denom * z_over_x * correction


def test_sabr_matches_hagan_expansion_away_from_the_money():
    f = 0.03
    t = 5.0
    rho = -0.3
    nu = 0.4

    for beta in [0.0, 0.3, 0.5, 0.7, 1.0]:
        alpha = 0.25 * f ** (1.0 - beta)
        for k in [0.01, 0.02, 0.03, 0.05, 0.08]:
            params = np.array([alpha, beta, rho, nu])
            vol = vol_function_sabr(params, f, k, t)
            expected = _hagan_black_vol(alpha, beta, rho, nu, f, k, t)
            assert abs(vol - expected) < 1e-12


def test_sabr_normal_limit_matches_bachelier():
    # With beta = 0 and vanishing vol of vol, SABR is the Bachelier model
    # with normal volatility alpha. The Black vol that reprices the
    # Bachelier call must agree with the expansion to within its error.
    f = 0.03
    t = 1.0
    alpha = 0.0075
    model = SABR(alpha, 0.0, 0.0, 1e-8)

    for k in [0.01, 0.02, 0.045, 0.08]:
        d = (f - k) / (alpha * np.sqrt(t))
        bachelier = (f - k) * norm.cdf(d) + alpha * np.sqrt(t) * norm.pdf(d)
        black_vol = brentq(lambda v: _black_call(f, k, t, v) - bachelier, 1e-4, 5.0)
        assert abs(model.black_vol(f, k, t) - black_vol) < 5e-4


def _black_call(f, k, t, v):
    d1 = (np.log(f / k) + 0.5 * v * v * t) / (v * np.sqrt(t))
    return f * norm.cdf(d1) - k * norm.cdf(d1 - v * np.sqrt(t))
