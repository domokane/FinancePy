##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..utils.error import FinError

from ..utils.global_types import OptionTypes
from ..utils.math import normcdf

from .equity_asian_option_mc import equity_asian_value_mc_fast_cv_numba
from .equity_asian_option_mc import equity_asian_value_mc_fast_numba
from .equity_asian_option_mc import equity_asian_value_mc_numba
from .equity_asian_option_mc import error_str


def value_geometric(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price,
    r,
    q,
    model,
    accrued_average,
):
    """This option valuation is based on paper by Kemna and Vorst 1990. It
    calculates the Geometric Asian option price which is a lower bound on
    the Arithmetic option price. This should not be used as a valuation
    model for the Arithmetic Average option but can be used as a control
    variate for other approaches."""

    # the years to the start of the averaging period
    tau = t_exp - t_avg

    volatility = model.volatility

    n = num_obs
    s0 = stock_price

    multiplier = 1.0

    if t_avg < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t_avg) / t_exp
        # the number of options is rescaled also
        multiplier = t_exp / tau
        # there is no pre-averaging time
        t_avg = 0.0
        # the number of observations is scaled
        n = n * t_exp / tau

    sig_sq = volatility**2
    mean_geo = (r - q - sig_sq / 2.0) * (t_avg + (t_exp - t_avg) / 2.0)
    var_geo = sig_sq * (t_avg + (t_exp - t_avg) * (2 * n - 1) / (6 * n))
    eg = s0 * np.exp(mean_geo + var_geo / 2.0)

    if np.abs(var_geo) < 1e-10:
        raise FinError("Asian option geometric variance is zero.")

    d1 = (mean_geo + np.log(s0 / k) + var_geo) / np.sqrt(var_geo)
    d2 = d1 - np.sqrt(var_geo)

    # the Geometric price is the lower bound
    call_g = np.exp(-r * t_exp) * (eg * normcdf(d1) - k * normcdf(d2))

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        v = call_g
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        put_g = call_g - (eg - k) * np.exp(-r * t_exp)
        v = put_g
    else:
        raise FinError("Unknown OPTION_TYPE " + str(opt_type_value))

    v = v * multiplier
    return v


####################################################################################


def value_curran(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price,
    r,
    q,
    model,
    accrued_average,
):
    """Valuation of an Asian option using the result by Vorst."""

    tau = t_exp - t_avg

    multiplier = 1.0

    volatility = model.volatility

    s0 = stock_price
    b = r - q
    sigma2 = volatility**2
    n = num_obs

    if t_avg < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t_avg) / t_exp
        # the number of options is rescaled also
        multiplier = t_exp / tau
        # there is no pre-averaging time
        t_avg = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t_exp / tau + 0.5) + 1

    h = (t_exp - t_avg) / (n - 1)
    u = (1.0 - np.exp(b * h * n)) / (1.0 - np.exp(b * h))
    w = (1.0 - np.exp((2 * b + sigma2) * h * n)) / (1.0 - np.exp((2 * b + sigma2) * h))

    fa = (s0 / n) * np.exp(b * t_avg) * u
    ea2 = (s0 * s0 / n / n) * np.exp((2.0 * b + sigma2) * t_avg)
    ea2 = ea2 * (w + 2.0 / (1.0 - np.exp((b + sigma2) * h)) * (u - w))
    sigma_aa = np.sqrt((np.log(ea2) - 2.0 * np.log(fa)) / t_exp)

    d1 = (np.log(fa / k) + sigma_aa * sigma_aa * t_exp / 2.0) / (sigma_aa * np.sqrt(t_exp))
    d2 = d1 - sigma_aa * np.sqrt(t_exp)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        v = np.exp(-r * t_exp) * (fa * normcdf(d1) - k * normcdf(d2))
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        v = np.exp(-r * t_exp) * (k * normcdf(-d2) - fa * normcdf(-d1))
    else:
        return None

    v = v * multiplier
    return v


####################################################################################


def value_turnbull_wakeman(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price,
    r,
    q,
    model,
    accrued_average,
):
    """Asian option valuation based on paper by Turnbull and Wakeman 1991
    which uses the edgeworth expansion to find the first two moments of the
    arithmetic average."""

    tau = t_exp - t_avg

    multiplier = 1.0

    volatility = model.volatility

    if t_avg < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t_avg) / t_exp
        # the number of options is rescaled also
        multiplier = t_exp / tau
        # there is no pre-averaging time
        t_avg = 0.0
        # the number of observations is scaled and floored at 1
        # n = int(num_obs * t_exp / tau + 0.5) + 1

    # need to handle this
    b = r - q
    sigma2 = volatility**2
    a1 = b + sigma2
    a2 = 2 * b + sigma2
    s0 = stock_price

    dt = t_exp - t_avg

    if b == 0:
        m1 = 1.0
        m2 = 2.0 * np.exp(sigma2 * t_exp) - 2.0 * np.exp(sigma2 * t_avg) * (1.0 + sigma2 * dt)
        m2 = m2 / sigma2 / sigma2 / dt / dt
    else:
        m1 = s0 * (np.exp(b * t_exp) - np.exp(b * t_avg)) / (b * dt)
        m2 = np.exp(a2 * t_exp) / a1 / a2 / dt / dt + (np.exp(a2 * t_avg) / b / dt / dt) * (
            1.0 / a2 - np.exp(b * dt) / a1
        )
        m2 = 2.0 * m2 * s0 * s0

    f0 = m1
    sigma2 = (1.0 / t_exp) * np.log(m2 / m1 / m1)
    sigma = np.sqrt(sigma2)

    d1 = (np.log(f0 / k) + sigma2 * t_exp / 2) / sigma / np.sqrt(t_exp)
    d2 = d1 - sigma * np.sqrt(t_exp)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        call = np.exp(-r * t_exp) * (f0 * normcdf(d1) - k * normcdf(d2))
        v = call
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        put = np.exp(-r * t_exp) * (k * normcdf(-d2) - f0 * normcdf(-d1))
        v = put
    else:
        return None

    v = v * multiplier
    return v


####################################################################################


def value_mc(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price: float,
    r: float,
    q: float,
    model,
    num_paths: int,
    seed: int,
    accrued_average: float,
):
    """Monte Carlo valuation of the Asian Average option using standard
    Monte Carlo code enhanced by Numba. I have discontinued the use of this
    as it is both slow and has limited variance reduction."""

    volatility = model.volatility

    v = equity_asian_value_mc_numba(
        t_avg,
        t_exp,
        k,
        num_obs,
        opt_type_value,
        stock_price,
        r,
        q,
        volatility,
        num_paths,
        seed,
        accrued_average,
    )

    return v


####################################################################################


def value_mc_fast(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price,
    r: float,
    q: float,
    model,  # Model
    num_paths,  # Numpaths integer
    seed,
    accrued_average,
):
    """Monte Carlo valuation of the Asian Average option. This method uses
    a lot of Numpy vectorisation. It is also helped by Numba."""

    tau = t_exp - t_avg

    n = num_obs

    volatility = model.volatility

    v = equity_asian_value_mc_fast_numba(
        t_avg,
        t_exp,
        tau,
        k,
        n,
        opt_type_value,
        stock_price,
        r,
        q,
        volatility,
        num_paths,
        seed,
        accrued_average,
    )

    return v


####################################################################################


def value_mc_fast_vc_numba(
    t_avg,
    t_exp,
    k,
    num_obs,
    opt_type_value,
    stock_price: float,
    r: float,
    q: float,
    model,
    num_paths: int,
    seed: int,
    accrued_average: float,
):
    """Monte Carlo valuation of the Asian Average option using a control
    variate method that improves accuracy and reduces the variance of the
    price. This uses Numpy and Numba. This is the standard MC pricer."""

    tau = t_exp - t_avg
    n = num_obs

    volatility = model.volatility

    # For control variate we price a Geometric average option exactly
    v_g_exact = value_geometric(
        t_avg,
        t_exp,
        stock_price,
        r,
        q,
        model,
        accrued_average,
    )

    v = equity_asian_value_mc_fast_cv_numba(
        t_avg,
        t_exp,
        tau,
        k,
        n,
        opt_type_value,
        stock_price,
        r,
        q,
        volatility,
        num_paths,
        seed,
        accrued_average,
        v_g_exact,
    )

    return v
