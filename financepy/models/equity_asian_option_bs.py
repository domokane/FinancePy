##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..utils.error import FinError

from ..utils.global_types import OptionTypes
from ..utils.math import normcdf

from .equity_asian_option_mc import error_str


def value_geometric(
    t_avg,
    t_exp,
    k,
    num_obs_per_year,
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
    vol2 = volatility**2
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

    averaging_time = t_exp - t_avg
    n = max(1, int(averaging_time * num_obs_per_year + 0.5))
    dt = averaging_time / n

    mean_time = t_avg + dt * (n + 1) / 2.0

    variance_time = (
        t_avg
        + dt * (n + 1) * (2 * n + 1) / (6.0 * n)
    )

    mean_geo = (
        r - q - vol2 / 2.0
    ) * mean_time

    var_geo = vol2 * variance_time
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
    num_obs_per_year,
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

    if t_avg < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t_avg) / t_exp
        # the number of options is rescaled also
        multiplier = t_exp / tau
        # there is no pre-averaging time
        t_avg = 0.0

    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    # Observation times are:
    #
    #     t_avg + h, ..., t_avg + n*h = t_exp
    #
    # so the first observation occurs one time step after t_avg.

    h = averaging_time / n
    t0 = t_avg + h

    bh = b * h
    bsh = (b + sigma2) * h
    b2sh = (2.0 * b + sigma2) * h

    exp_bh = np.exp(bh)
    exp_bhn = np.exp(bh * n)

    exp_bsh = np.exp(bsh)

    exp_b2sh = np.exp(b2sh)
    exp_b2shn = np.exp(b2sh * n)

    u = (1.0 - exp_bhn) / (1.0 - exp_bh)
    w = (1.0 - exp_b2shn) / (1.0 - exp_b2sh)

    fa = (s0 / n) * np.exp(b * t0) * u

    ea2 = (
        (s0 * s0 / (n * n))
        * np.exp((2.0 * b + sigma2) * t0)
        * (
            w
            + 2.0 / (1.0 - exp_bsh) * (u - w)
        )
    )

    if ea2 < fa * fa:
        raise FinError(
            "Curran second moment is less than squared first moment."
        )

    var_aa = np.log(ea2 / (fa * fa))
    sqrt_var_aa = np.sqrt(var_aa)

    d1 = (np.log(fa / k) + 0.5 * var_aa) / sqrt_var_aa
    d2 = d1 - sqrt_var_aa

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
    num_obs_per_year,
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
        m1 = s0
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
    var_a = np.log(m2 / (m1 * m1))
    sqrt_var_a = np.sqrt(var_a)

    d1 = (
        np.log(f0 / k)
        + 0.5 * var_a
    ) / sqrt_var_a

    d2 = d1 - sqrt_var_a

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


# def value_mc(
#     t_avg,
#     t_exp,
#     k,
#     num_obs_per_year,
#     opt_type_value,
#     stock_price: float,
#     r: float,
#     q: float,
#     model,
#     num_paths: int,
#     seed: int,
#     accrued_average: float,
# ):
#     """Monte Carlo valuation of the Asian Average option using standard
#     Monte Carlo code enhanced by Numba. I have discontinued the use of this
#     as it is both slow and has limited variance reduction."""

#     volatility = model.volatility

#     v = equity_asian_value_mc_numba(
#         t_avg,
#         t_exp,
#         k,
#         num_obs_per_year,
#         opt_type_value,
#         stock_price,
#         r,
#         q,
#         volatility,
#         num_paths,
#         seed,
#         accrued_average,
#     )

#     return v


# ####################################################################################


# def value_mc_fast(
#     t_avg,
#     t_exp,
#     k,
#     num_obs_per_year,
#     opt_type_value,
#     stock_price,
#     r: float,
#     q: float,
#     volatility,  # Model
#     num_paths,  # Numpaths integer
#     seed,
#     accrued_average,
# ):
#     """Monte Carlo valuation of the Asian Average option. This method uses
#     a lot of Numpy vectorisation. It is also helped by Numba."""

#     v = equity_asian_value_mc_fast_numba(
#         t_avg,
#         t_exp,
#         k,
#         num_obs_per_year,
#         opt_type_value,
#         stock_price,
#         r,
#         q,
#         volatility,
#         num_paths,
#         seed,
#         accrued_average,
#     )

#     return v


# ####################################################################################


# def value_mc_fast_vc_numba(
#     t_avg,
#     t_exp,
#     k,
#     num_obs_per_year,
#     opt_type_value,
#     stock_price: float,
#     r: float,
#     q: float,
#     model,
#     num_paths: int,
#     seed: int,
#     accrued_average: float,
# ):
#     """Monte Carlo valuation of the Asian Average option using a control
#     variate method that improves accuracy and reduces the variance of the
#     price. This uses Numpy and Numba. This is the standard MC pricer."""

#     volatility = model.volatility

#     # For control variate we price a Geometric average option exactly
#     v_g_exact = value_geometric(
#         t_avg,
#         t_exp,
#         k,
#         num_obs_per_year,
#         opt_type_value,
#         stock_price,
#         r,
#         q,
#         model,
#         accrued_average,
#     )

#     v = equity_asian_value_mc_fast_cv_numba(
#         t_avg,
#         t_exp,
#         k,
#         num_obs_per_year,
#         opt_type_value,
#         stock_price,
#         r,
#         q,
#         volatility,
#         num_paths,
#         seed,
#         accrued_average,
#         v_g_exact,
#     )

#     return v
