##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from .asian_option_mc import error_str
from ..utils.math import normcdf
from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from numba import njit
import numpy as np

print("DO NOT USE THIS MODULE - USE ASIAN_OPTION_BS.py !!!")


def _geom_sum(x, n):
    """Return sum(exp(j*x), j=0,...,n-1) stably."""
    if abs(x) < 1e-10:
        return float(n)
    return np.expm1(n * x) / np.expm1(x)


def value_asian_kemna_vorst_geometric(
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
    """Price a discretely sampled geometric-average Asian option
    under Black-Scholes/GBM using the Kemna-Vorst result.

    For otherwise identical fixed-strike contracts, the geometric-average
    call is a lower bound on the arithmetic-average call, while the
    geometric-average put is an upper bound on the arithmetic-average put.
    The geometric price may also be used as a control variate for Monte
    Carlo pricing of arithmetic-average Asian options.
    """

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
    variance_time = (t_avg + dt * (n + 1) * (2 * n + 1) / (6.0 * n))

    mean_geo = (r - q - vol2 / 2.0) * mean_time
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


def value_asian_curran(
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
    """Price a discretely sampled arithmetic-average Asian option
    using Curran's conditioning approximation.

    Curran's approximation conditions the arithmetic average on the
    geometric average.

    Observation times are

        t_avg + h, ..., t_avg + n*h = t_exp.

    Reference:
        Curran (1994), "Valuing Asian and Portfolio Options by
        Conditioning on the Geometric Mean Price".
    """

    tau = t_exp - t_avg
    multiplier = 1.0

    sigma = model.volatility
    sigma2 = sigma * sigma

    s0 = stock_price
    b = r - q

    # ------------------------------------------------------------
    # Already inside the averaging period.
    # ------------------------------------------------------------

    if t_avg < 0.0:

        if accrued_average is None:
            raise FinError(error_str)

        k = (k * tau + accrued_average * t_avg) / t_exp
        multiplier = t_exp / tau
        t_avg = 0.0

    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    h = averaging_time / n

    # Observation times:
    #
    # t_i = t_avg + (i + 1) * h
    #
    # for i = 0, ..., n-1.

    times = np.empty(n)

    for i in range(n):
        times[i] = t_avg + (i + 1) * h

    # ------------------------------------------------------------
    # Distribution of Y = log(G)
    #
    # G = geometric average
    #
    # Y = (1/n) sum log(S_i)
    # ------------------------------------------------------------

    mean_t = 0.0

    for i in range(n):
        mean_t += times[i]

    mean_t /= n

    mu_y = (
        np.log(s0)
        + (b - 0.5 * sigma2) * mean_t
    )

    # ------------------------------------------------------------
    # Calculate
    #
    # Cov(log(S_i), Y)
    #
    # where
    #
    # Cov(log(S_i), log(S_j))
    #     = sigma^2 min(t_i, t_j).
    #
    # Since the observation times are ordered, this can be evaluated
    # in O(n) without constructing an n x n covariance matrix.
    # ------------------------------------------------------------

    cov_i_y = np.empty(n)

    for i in range(n):

        ti = times[i]

        # sum(times[0:i+1])
        #
        # times[j] = t_avg + (j+1) h
        #
        # so
        #
        # sum = (i+1)t_avg
        #       + h (i+1)(i+2)/2

        sum_to_i = (
            (i + 1) * t_avg
            + h * (i + 1) * (i + 2) / 2.0
        )

        # For j > i:
        #
        # min(t_i, t_j) = t_i.

        sum_min = (
            sum_to_i
            + (n - i - 1) * ti
        )

        cov_i_y[i] = sigma2 * sum_min / n

    # Var(Y) = average_i Cov(log(S_i), Y)

    var_y = 0.0

    for i in range(n):
        var_y += cov_i_y[i]

    var_y /= n

    # ------------------------------------------------------------
    # Deterministic limit.
    # ------------------------------------------------------------

    if var_y < 1.0e-14:

        expected_a = 0.0

        for i in range(n):
            expected_a += s0 * np.exp(b * times[i])

        expected_a /= n

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:

            payoff = max(expected_a - k, 0.0)

        elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:

            payoff = max(k - expected_a, 0.0)

        else:

            raise FinError(
                "Unknown OPTION_TYPE " + str(opt_type_value)
            )

        return (
            multiplier
            * np.exp(-r * t_exp)
            * payoff
        )

    sqrt_var_y = np.sqrt(var_y)

    # ------------------------------------------------------------
    # Conditional expectation:
    #
    # E[S_i | Y=y]
    #
    # = alpha_i exp(beta_i y)
    #
    # where
    #
    # beta_i = Cov(log(S_i), Y) / Var(Y).
    # ------------------------------------------------------------

    beta = np.empty(n)
    alpha = np.empty(n)

    for i in range(n):

        beta[i] = cov_i_y[i] / var_y

        mu_i = (
            np.log(s0)
            + (b - 0.5 * sigma2) * times[i]
        )

        cond_var_i = (
            sigma2 * times[i]
            - cov_i_y[i] * cov_i_y[i] / var_y
        )

        # Protect against tiny negative values caused by
        # floating-point rounding.

        if cond_var_i < 0.0:

            if cond_var_i > -1.0e-14:
                cond_var_i = 0.0
            else:
                raise FinError(
                    "Negative conditional variance in Curran model."
                )

        alpha[i] = np.exp(
            mu_i
            - beta[i] * mu_y
            + 0.5 * cond_var_i
        )

    # ------------------------------------------------------------
    # Find the critical geometric-average level y_star satisfying
    #
    # E[A | Y=y_star] = K.
    #
    # The conditional arithmetic average is monotonic in y.
    # ------------------------------------------------------------

    lo = mu_y - 10.0 * sqrt_var_y
    hi = mu_y + 10.0 * sqrt_var_y

    # Expand lower bound if necessary.

    for _ in range(100):

        conditional_a = 0.0

        for i in range(n):
            conditional_a += alpha[i] * np.exp(beta[i] * lo)

        conditional_a /= n

        if conditional_a <= k:
            break

        lo -= 5.0 * sqrt_var_y

    # Expand upper bound if necessary.

    for _ in range(100):

        conditional_a = 0.0

        for i in range(n):
            conditional_a += alpha[i] * np.exp(beta[i] * hi)

        conditional_a /= n

        if conditional_a >= k:
            break

        hi += 5.0 * sqrt_var_y

    # ------------------------------------------------------------
    # Bisection.
    # ------------------------------------------------------------

    for _ in range(100):

        mid = 0.5 * (lo + hi)

        conditional_a = 0.0

        for i in range(n):
            conditional_a += alpha[i] * np.exp(beta[i] * mid)

        conditional_a /= n

        if conditional_a > k:
            hi = mid
        else:
            lo = mid

    y_star = 0.5 * (lo + hi)

    # ------------------------------------------------------------
    # Curran call approximation.
    #
    # For each observation:
    #
    # E[S_i 1(Y > y_star)]
    #
    # = E[S_i] *
    #   N((mu_Y + Cov(log(S_i),Y) - y_star) / sigma_Y)
    # ------------------------------------------------------------

    weighted_sum = 0.0
    expected_a = 0.0

    for i in range(n):

        expected_si = s0 * np.exp(b * times[i])

        d_i = (
            mu_y
            + cov_i_y[i]
            - y_star
        ) / sqrt_var_y

        # FinancePy normcdf is scalar.
        weighted_sum += expected_si * normcdf(d_i)

        expected_a += expected_si

    weighted_sum /= n
    expected_a /= n

    d_k = (
        mu_y
        - y_star
    ) / sqrt_var_y

    df = np.exp(-r * t_exp)

    call = df * (
        weighted_sum
        - k * normcdf(d_k)
    )

    # ------------------------------------------------------------
    # Put-call parity for an arithmetic-average Asian:
    #
    # C - P = DF * (E[A] - K)
    # ------------------------------------------------------------

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:

        v = call

    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:

        v = call - df * (expected_a - k)

    else:

        raise FinError(
            "Unknown OPTION_TYPE " + str(opt_type_value)
        )

    return multiplier * v

####################################################################################


def value_asian_turnbull_wakeman_discrete(
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
    """Discrete-observation version of the Turnbull-Wakeman
    two-moment lognormal approximation.

    Exact first and second moments of the discretely sampled
    arithmetic average are matched to a lognormal distribution.

    Observation times are

        t_avg + h, ..., t_avg + n*h = t_exp.
    """

    tau = t_exp - t_avg
    multiplier = 1.0

    sigma = model.volatility
    sigma2 = sigma * sigma

    s0 = stock_price
    b = r - q

    if t_avg < 0.0:

        if accrued_average is None:
            raise FinError(error_str)

        k = (k * tau + accrued_average * t_avg) / t_exp
        multiplier = t_exp / tau
        t_avg = 0.0

    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    h = averaging_time / n

    # ---------------------------------------------------------
    # Exact first and second moments of the discrete
    # arithmetic average.
    # ---------------------------------------------------------

    m1 = 0.0
    m2 = 0.0

    for i in range(n):

        ti = t_avg + (i + 1) * h

        m1 += np.exp(b * ti)

        for j in range(n):

            tj = t_avg + (j + 1) * h

            m2 += np.exp(
                b * (ti + tj)
                + sigma2 * min(ti, tj)
            )

    m1 = s0 * m1 / n
    m2 = s0 * s0 * m2 / (n * n)

    # ---------------------------------------------------------
    # Match a lognormal distribution to m1 and m2.
    # ---------------------------------------------------------

    ratio = m2 / (m1 * m1)

    if ratio < 1.0:

        if ratio > 1.0 - 1.0e-12:
            ratio = 1.0
        else:
            raise FinError(
                "Asian second moment is less than "
                "squared first moment."
            )

    var_a = np.log(ratio)

    df = np.exp(-r * t_exp)

    # Deterministic limit.
    if var_a < 1.0e-14:

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            v = df * max(m1 - k, 0.0)

        elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
            v = df * max(k - m1, 0.0)

        else:
            raise FinError(
                "Unknown OPTION_TYPE " + str(opt_type_value)
            )

        return multiplier * v

    sqrt_var_a = np.sqrt(var_a)

    d1 = (
        np.log(m1 / k)
        + 0.5 * var_a
    ) / sqrt_var_a

    d2 = d1 - sqrt_var_a

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:

        v = df * (
            m1 * normcdf(d1)
            - k * normcdf(d2)
        )

    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:

        v = df * (
            k * normcdf(-d2)
            - m1 * normcdf(-d1)
        )

    else:

        raise FinError(
            "Unknown OPTION_TYPE " + str(opt_type_value)
        )

    return multiplier * v

####################################################################################


def value_asian_turnbull_wakeman(
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
    """Approximate a continuously averaged arithmetic Asian option using
    the Turnbull-Wakeman first-two-moment lognormal approximation.
    """

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

    if abs(b) < 1.0e-12:

        m1 = s0

        x = sigma2 * dt

        # Stable evaluation of:
        #
        # phi(x) = (exp(x) - 1 - x) / x^2
        #
        # phi(0) = 1/2

        if abs(x) < 1.0e-8:
            phi = (
                0.5
                + x / 6.0
                + x * x / 24.0
                + x * x * x / 120.0
            )
        else:
            phi = (np.expm1(x) - x) / (x * x)

        m2 = (
            2.0
            * s0
            * s0
            * np.exp(sigma2 * t_avg)
            * phi
        )

    else:

        m1 = (
            s0
            * (np.exp(b * t_exp) - np.exp(b * t_avg))
            / (b * dt)
        )

        m2 = (
            np.exp(a2 * t_exp)
            / (a1 * a2 * dt * dt)
            + np.exp(a2 * t_avg)
            / (b * dt * dt)
            * (
                1.0 / a2
                - np.exp(b * dt) / a1
            )
        )

        m2 = 2.0 * m2 * s0 * s0

    f0 = m1

    ratio = m2 / (m1 * m1)

    # Protect against small floating-point errors.
    if ratio < 1.0:
        if ratio > 1.0 - 1.0e-12:
            ratio = 1.0
        else:
            raise FinError(
                "Turnbull-Wakeman second moment is less than "
                "squared first moment."
            )

    var_a = np.log(ratio)

    # Deterministic / zero-variance limit.
    if var_a < 1.0e-14:

        df = np.exp(-r * t_exp)

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            v = df * max(f0 - k, 0.0)

        elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
            v = df * max(k - f0, 0.0)

        else:
            raise FinError(
                "Unknown OPTION_TYPE " + str(opt_type_value)
            )

        return multiplier * v

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
