from enum import Enum
from numba import njit
import numpy as np

from ..utils.error import FinError
from ..utils.math import phi2
from ..utils.global_vars import G_SMALL

from ..utils.global_types import OptionTypes
from ..utils.math import normcdf, normcdf_vect

from .finite_difference import black_scholes_fd
from .bs_european_analytic import european_value

from .equity_lsmc import equity_lsmc, BoundaryFitTypes
from .equity_crr_tree import crr_tree_val_avg
from .equity_crr_tree import crr_richardson
from .equity_crr_tree import leisen_reimer_tree


class BlackScholesAlgorithmTypes(Enum):

    DEFAULT = 0
    ANALYTICAL = 1
    CRR_TREE = 2
    BARONE_ADESI_WHALEY = 3
    LSMC = 4
    BJERKSUND_STENSLAND_1993 = 5
    BJERKSUND_STENSLAND_2002 = 6
    FINITE_DIFFERENCE = 7
    CRR_RICHARDSON = 8
    LEISEN_REIMER = 9
#    PSOR = 9 # VERY SLOW


########################################################################################


def american_value(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type: OptionTypes,
    algo_type: BlackScholesAlgorithmTypes = BlackScholesAlgorithmTypes.DEFAULT,
    **params,
) -> float:
    """Compute an option value using the selected Black-Scholes algorithm."""

    if not isinstance(algo_type, BlackScholesAlgorithmTypes):
        raise FinError(
            f"Invalid Black-Scholes algorithm type: {algo_type}"
        )

    if opt_type in (
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ):

        if algo_type is BlackScholesAlgorithmTypes.DEFAULT:
            algo_type = BlackScholesAlgorithmTypes.ANALYTICAL

        return _value_european(
            algo_type,
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type.value,
            **params,
        )

    if opt_type in (
        OptionTypes.AMERICAN_CALL,
        OptionTypes.AMERICAN_PUT,
    ):

        if algo_type is BlackScholesAlgorithmTypes.DEFAULT:
            algo_type = BlackScholesAlgorithmTypes.CRR_TREE

        return _value_american(
            algo_type,
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type.value,
            **params,
        )

    raise FinError(f"Unsupported option type: {opt_type}")


########################################################################################


def _value_european(
    algo_type: BlackScholesAlgorithmTypes,
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    **params,
) -> float:

    if algo_type is BlackScholesAlgorithmTypes.ANALYTICAL:
        return european_value(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
        )

    if algo_type is BlackScholesAlgorithmTypes.CRR_TREE:
        return _value_crr_tree(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.CRR_RICHARDSON:
        return _value_crr_tree_richardson(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.LEISEN_REIMER:
        return _value_leisen_reimer(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.LSMC:
        return _value_lsmc(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.FINITE_DIFFERENCE:
        return _value_fd(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    # if algo_type is BlackScholesAlgorithmTypes.PSOR:
    #     return _value_psor(
    #         spot_price,
    #         time_to_expiry,
    #         strike_price,
    #         risk_free_rate,
    #         dividend_rate,
    #         volatility,
    #         opt_type_value,
    #         **params,
    #     )

    raise FinError(
        f"Unsupported European Black-Scholes algorithm: {algo_type}"
    )


########################################################################################


def _value_american(
    algo_type: BlackScholesAlgorithmTypes,
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    **params,
) -> float:

    if algo_type is BlackScholesAlgorithmTypes.BARONE_ADESI_WHALEY:
        return baw_value(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
        )

    if algo_type is BlackScholesAlgorithmTypes.BJERKSUND_STENSLAND_1993:
        return bjerksund_stensland_1993_value(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
        )

    if algo_type is BlackScholesAlgorithmTypes.BJERKSUND_STENSLAND_2002:
        return bjerksund_stensland_2002_value(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
        )

    if algo_type is BlackScholesAlgorithmTypes.CRR_TREE:
        return _value_crr_tree(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.CRR_RICHARDSON:
        return _value_crr_tree_richardson(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.LEISEN_REIMER:
        return _value_leisen_reimer_tree(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.LSMC:
        return _value_lsmc(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    if algo_type is BlackScholesAlgorithmTypes.FINITE_DIFFERENCE:
        return _value_fd(
            spot_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_rate,
            volatility,
            opt_type_value,
            **params,
        )

    # if algo_type is BlackScholesAlgorithmTypes.PSOR:
    #     return _value_psor(
    #         spot_price,
    #         time_to_expiry,
    #         strike_price,
    #         risk_free_rate,
    #         dividend_rate,
    #         volatility,
    #         opt_type_value,
    #         **params,
    #     )

    raise FinError(
        f"Unsupported American Black-Scholes algorithm: {algo_type}"
    )


########################################################################################


def _value_crr_tree(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    num_steps_per_year: int = 200,
    do_averaging: bool = True,
    ** _,
) -> float:

    return crr_tree_val_avg(
        spot_price,
        risk_free_rate,
        dividend_rate,
        volatility,
        num_steps_per_year,
        time_to_expiry,
        opt_type_value,
        strike_price,
        do_averaging
    )["value"]


########################################################################################

def _value_crr_tree_richardson(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    num_steps_per_year: int = 200,
    ** _,
) -> float:

    return crr_richardson(
        spot_price,
        risk_free_rate,
        dividend_rate,
        volatility,
        num_steps_per_year,
        time_to_expiry,
        opt_type_value,
        strike_price,
    )["value"]

########################################################################################


def _value_leisen_reimer(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    num_steps_per_year: int = 200,
    ** _,
) -> float:

    return leisen_reimer_tree(
        spot_price,
        risk_free_rate,
        dividend_rate,
        volatility,
        num_steps_per_year,
        time_to_expiry,
        opt_type_value,
        strike_price,
    )["value"]


########################################################################################


def _value_lsmc(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    num_steps_per_year: int = 100,
    num_paths: int = 10_000,
    poly_degree: int = 3,
    fit_type: BoundaryFitTypes = BoundaryFitTypes.HERMITE_E,
    use_sobol: bool = False,
    seed: int = 42,
    **_,
) -> float:

    return equity_lsmc(
        spot_price=spot_price,
        risk_free_rate=risk_free_rate,
        dividend_yield=dividend_rate,
        sigma=volatility,
        num_steps_per_year=num_steps_per_year,
        num_paths=num_paths,
        time_to_expiry=time_to_expiry,
        opt_type_value=opt_type_value,
        strike_price=strike_price,
        poly_degree=poly_degree,
        fit_type_value=fit_type.value,
        use_sobol=use_sobol,
        seed=seed,
    )


########################################################################################


def _value_fd(
    spot_price: float,
    time_to_expiry: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_rate: float,
    volatility: float,
    opt_type_value: int,
    **params,
) -> float:

    return black_scholes_fd(
        spot_price=spot_price,
        time_to_expiry=time_to_expiry,
        strike_price=strike_price,
        risk_free_rate=risk_free_rate,
        dividend_yield=dividend_rate,
        volatility=volatility,
        opt_type_value=opt_type_value,
        **params,
    )


########################################################################################

# def _value_psor(
#     spot_price: float,
#     time_to_expiry: float,
#     strike_price: float,
#     risk_free_rate: float,
#     dividend_rate: float,
#     volatility: float,
#     opt_type_value: int,
#     **params,
# ) -> float:

#     return black_scholes_fd_psor(
#         spot_price=spot_price,
#         time_to_expiry=time_to_expiry,
#         strike_price=strike_price,
#         risk_free_rate=risk_free_rate,
#         dividend_yield=dividend_rate,
#         volatility=volatility,
#         opt_type_value=opt_type_value,
#         **params,
#     )

@njit(cache=True, fastmath=True)
def _baw_objective(s, t, k, r, q, v, phi):

    if phi == 1:
        return _fcall(s, t, k, r, q, v)

    return _fput(s, t, k, r, q, v)

########################################################################################


@njit(cache=True, fastmath=True)
def _baw_secant(x0, t, k, r, q, v, phi, tol=1.0e-7, maxiter=50):

    eps = 1.0e-4

    p0 = x0
    p1 = x0 * (1.0 + eps)

    if p1 > 0.0:
        p1 += eps
    else:
        p1 -= eps

    q0 = _baw_objective(p0, t, k, r, q, v, phi)
    q1 = _baw_objective(p1, t, k, r, q, v, phi)

    if np.abs(q1) < np.abs(q0):
        p0, p1, q0, q1 = p1, p0, q1, q0

    for _ in range(maxiter):

        if q1 == q0:
            p = (p1 + p0) / 2.0

            if np.abs(q1) <= tol or np.abs(p1 - p0) <= tol:
                return p

            raise FinError("Secant slope is zero before convergence")

        if np.abs(q1) > np.abs(q0):
            p = (-q0 / q1 * p1 + p0) / (1.0 - q0 / q1)
        else:
            p = (-q1 / q0 * p0 + p1) / (1.0 - q1 / q0)

        if np.abs(p - p1) < tol:
            return p

        p0, q0 = p1, q1
        p1 = p
        q1 = _baw_objective(p1, t, k, r, q, v, phi)

    raise FinError("Failed to converge")

########################################################################################


@njit(fastmath=True, cache=True)
def _baw_mm_over_kk(r, t, v2):
    """Evaluate (2r/v^2) / (1-exp(-rt)), including its r=0 limit."""

    if r == 0.0:
        return 2.0 / (v2 * t)

    return 2.0 * r / (v2 * -np.expm1(-r * t))


########################################################################################


@njit(fastmath=True, cache=True)
def _fcall(si, *args):
    """Function to determine ststar for pricing American call options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v * v

    ww = 2.0 * b / v2
    mm_over_kk = _baw_mm_over_kk(r, t, v2)

    q2 = (1.0 - ww + np.sqrt((ww - 1.0) ** 2 + 4.0 * mm_over_kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    obj_fn = si - k
    obj_fn = obj_fn - european_value(si, t, k, r, q, v, OptionTypes.EUROPEAN_CALL.value)
    obj_fn = obj_fn - (1.0 - np.exp(-q * t) * normcdf_vect(d1)) * si / q2
    return obj_fn


########################################################################################


@njit(fastmath=True, cache=True)
def _fput(si, *args):
    """Function to determine sstar for pricing American put options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v * v

    ww = 2.0 * b / v2
    mm_over_kk = _baw_mm_over_kk(r, t, v2)

    q1 = (1.0 - ww - np.sqrt((ww - 1.0) ** 2 + 4.0 * mm_over_kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    obj_fn = si - k
    obj_fn = obj_fn + european_value(
        si,
        t,
        k,
        r,
        q,
        v,
        OptionTypes.EUROPEAN_PUT.value,
    )
    obj_fn = obj_fn - (1.0 - np.exp(-q * t) * normcdf_vect(-d1)) * si / q1
    return obj_fn


########################################################################################

@njit(cache=True, fastmath=True)
def baw_value(s, t, k, r, q, v, opt_type_value):
    """American Option Pricing Approximation using the Barone-Adesi-wwhaley
    approximation for the Black-Scholes mmodel"""

    b = r - q
    G_SMALL = 1e-10

    if opt_type_value == OptionTypes.AMERICAN_CALL.value:

        if t <= G_SMALL:
            return max(s - k, 0.0)

        euro_type = OptionTypes.EUROPEAN_CALL.value

        if b >= r:
            return european_value(s, t, k, r, q, v, euro_type)

        sstar = _baw_secant(s, t, k, r, q, v, 1)

        v2 = v * v
        ww = 2.0 * b / v2
        mm_over_kk = _baw_mm_over_kk(r, t, v2)
        d1 = (np.log(sstar / k) + (b + v * v / 2.0) * t) / (v * np.sqrt(t))
        q2 = (
            -1.0 * (ww - 1.0)
            + np.sqrt((ww - 1.0) ** 2 + 4.0 * mm_over_kk)
        ) / 2.0
        a2 = (sstar / q2) * (1.0 - np.exp(-q * t) * normcdf_vect(d1))

        if s < sstar:
            return european_value(s, t, k, r, q, v, euro_type) + a2 * ((s / sstar) ** q2)

        return s - k

    elif opt_type_value == OptionTypes.AMERICAN_PUT.value:

        if t <= G_SMALL:
            return max(k - s, 0.0)

        euro_type = OptionTypes.EUROPEAN_PUT.value

        # With zero carry on the strike and a non-negative dividend yield,
        # the European put already dominates immediate exercise.
        if r == 0.0 and q >= 0.0:
            return european_value(s, t, k, r, q, v, euro_type)

        sstar = _baw_secant(s, t, k, r, q, v, -1)

        v2 = v * v

        ww = 2.0 * b / v2
        mm_over_kk = _baw_mm_over_kk(r, t, v2)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (
            -1.0 * (ww - 1.0)
            - np.sqrt((ww - 1.0) ** 2 + 4.0 * mm_over_kk)
        ) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * normcdf_vect(-d1))

        if s > sstar:
            bsv = european_value(s, t, k, r, q, v, euro_type)
            return bsv + a1 * ((s / sstar) ** q1)
        else:
            return k - s

    return np.nan


########################################################################################

@njit(cache=True, fastmath=True)
def _bs93_phi(s, t, gamma, h, i, r, b, v):
    """Equation (13) in the Bjerksund–Stensland 1993 approximation."""

    v2 = v * v
    sqrt_t = np.sqrt(t)

    lambda_ = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * v2) * t

    d = -(np.log(s / h) + (b + (gamma - 0.5) * v2) * t) / (v * sqrt_t)

    kappa = 2.0 * gamma - 1.0 + 2.0 * b / v2

    reflected_d = d - 2.0 * np.log(i / s) / (v * sqrt_t)

    return np.exp(lambda_) * s**gamma * (normcdf(d) - (i / s) ** kappa * normcdf(reflected_d))


########################################################################################


@njit(cache=True, fastmath=True)
def bjerksund_stensland_1993_value(
    s,
    t,
    k,
    r,
    q,
    v,
    opt_type_value,
):
    """Bjerksund–Stensland 1993 American-option approximation.

    This implementation assumes:
        s > 0
        k > 0
        t >= 0
        v > 0
        r >= 0
        q >= 0
    """

    is_call = opt_type_value == OptionTypes.AMERICAN_CALL.value
    is_put = opt_type_value == OptionTypes.AMERICAN_PUT.value

    if not is_call and not is_put:
        return np.nan

    # Handle expiry before performing the put-call transformation.
    if t <= G_SMALL:
        if is_call:
            return max(s - k, 0.0)

        return max(k - s, 0.0)

    if s <= 0.0 or k <= 0.0 or v <= 0.0:
        return np.nan

    # This simple implementation does not handle the negative-rate
    # double-boundary regimes.
    if r < 0.0 or q < 0.0:
        return np.nan

    if is_put:
        # American put-call transformation:
        #
        # P(S, K, r, q) = C(K, S, q, r)
        transformed_s = k
        transformed_k = s
        transformed_r = q
        transformed_q = r

        s = transformed_s
        k = transformed_k
        r = transformed_r
        q = transformed_q

    b = r - q
    v2 = v * v
    sqrt_t = np.sqrt(t)

    euro_call_type = OptionTypes.EUROPEAN_CALL.value

    # After the put-call transformation we are always pricing a call.
    #
    # For q <= 0, early exercise of the call is not optimal in the
    # standard non-negative-rate setting.
    if q <= 0.0:
        return european_value(
            s,
            t,
            k,
            r,
            q,
            v,
            euro_call_type,
        )

    beta = 0.5 - b / v2 + np.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)

    if not np.isfinite(beta) or beta <= 1.0:
        return np.nan

    b_infinity = k * beta / (beta - 1.0)

    # Correct formula:
    #
    # B0 = max(K, rK / (r-b))
    #
    # and r-b = q.
    b_0 = max(k, k * r / q)

    denominator = b_infinity - b_0

    if denominator <= 0.0 or not np.isfinite(denominator):
        return np.nan

    h_t = -(b * t + 2.0 * v * sqrt_t) * b_0 / denominator

    x_t = b_0 + denominator * (-np.expm1(h_t))

    if not np.isfinite(x_t) or x_t <= k:
        return np.nan

    # The continuation formula is only valid below the trigger.
    if s >= x_t:
        return s - k

    alpha = (x_t - k) * x_t ** (-beta)

    value = (
        alpha * s**beta
        - alpha * _bs93_phi(s, t, beta, x_t, x_t, r, b, v)
        + _bs93_phi(s, t, 1.0, x_t, x_t, r, b, v)
        - _bs93_phi(s, t, 1.0, k, x_t, r, b, v)
        - k * _bs93_phi(s, t, 0.0, x_t, x_t, r, b, v)
        + k * _bs93_phi(s, t, 0.0, k, x_t, r, b, v)
    )

    european_value_transformed = european_value(
        s,
        t,
        k,
        r,
        q,
        v,
        euro_call_type,
    )

    intrinsic_value = max(s - k, 0.0)

    # Numerical safeguards. Under the transformation, these are also
    # the correct European and intrinsic lower bounds for the put.
    return max(
        value,
        european_value_transformed,
        intrinsic_value,
    )


###############################################################################


# @vectorize(
#     fastmath=True,
#     cache=True,
# )
# def value(
#     s: float,
#     t: float,
#     k: float,
#     r: float,
#     q: float,
#     v: float,
#     opt_type_value: int,
# ) -> float:
#     """Price European or American vanilla options."""

#     BAW = False  # True
#     #    BAW = True

#     if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
#         return european_value(s, t, k, r, q, v, opt_type_value)

#     if opt_type_value == OptionTypes.EUROPEAN_PUT.value:
#         return european_value(s, t, k, r, q, v, opt_type_value)

#     if opt_type_value == OptionTypes.AMERICAN_CALL.value:
#         if BAW:
#             return baw_value(s, t, k, r, q, v, opt_type_value)
#         else:
#             return bjerksund_stensland_1993_value(s, t, k, r, q, v, opt_type_value)

#     if opt_type_value == OptionTypes.AMERICAN_PUT.value:
#         if BAW:
#             return baw_value(s, t, k, r, q, v, opt_type_value)
#         else:
#             return bjerksund_stensland_1993_value(s, t, k, r, q, v, opt_type_value)

#     return np.nan


########################################################################################
########################################################################################

@njit(cache=True, fastmath=True)
def _bs02_phi(s, t, gamma, h, i, r, b, v):
    """Phi function used by Bjerksund-Stensland (2002)."""

    v2 = v * v
    sqrt_t = np.sqrt(t)

    lam = (
        -r
        + gamma * b
        + 0.5 * gamma * (gamma - 1.0) * v2
    )

    kappa = 2.0 * b / v2 + 2.0 * gamma - 1.0

    d = -(
        np.log(s / h)
        + (b + (gamma - 0.5) * v2) * t
    ) / (v * sqrt_t)

    reflected_d = (
        d
        - 2.0 * np.log(i / s) / (v * sqrt_t)
    )

    return (
        np.exp(lam * t)
        * s**gamma
        * (
            normcdf(d)
            - (i / s)**kappa * normcdf(reflected_d)
        )
    )

########################################################################################


@njit(cache=True, fastmath=True)
def _bs02_psi(
    s,
    t,
    gamma,
    h,
    i2,
    i1,
    t1,
    r,
    b,
    v,
):
    """Psi function used by Bjerksund-Stensland (2002)."""

    v2 = v * v

    sqrt_t = np.sqrt(t)
    sqrt_t1 = np.sqrt(t1)

    mu = b + (gamma - 0.5) * v2

    lam = (
        -r
        + gamma * b
        + 0.5 * gamma * (gamma - 1.0) * v2
    )

    kappa = (
        2.0 * b / v2
        + 2.0 * gamma
        - 1.0
    )

    rho = np.sqrt(t1 / t)

    # ------------------------------------------------------------
    # t1 arguments
    # ------------------------------------------------------------

    e1 = (
        np.log(s / i1)
        + mu * t1
    ) / (v * sqrt_t1)

    e2 = (
        np.log(i2 * i2 / (s * i1))
        + mu * t1
    ) / (v * sqrt_t1)

    e3 = (
        np.log(s / i1)
        - mu * t1
    ) / (v * sqrt_t1)

    e4 = (
        np.log(i2 * i2 / (s * i1))
        - mu * t1
    ) / (v * sqrt_t1)

    # ------------------------------------------------------------
    # T arguments
    # ------------------------------------------------------------

    f1 = (
        np.log(s / h)
        + mu * t
    ) / (v * sqrt_t)

    f2 = (
        np.log(i2 * i2 / (s * h))
        + mu * t
    ) / (v * sqrt_t)

    f3 = (
        np.log(i1 * i1 / (s * h))
        + mu * t
    ) / (v * sqrt_t)

    f4 = (
        np.log(
            s * i1 * i1
            / (h * i2 * i2)
        )
        + mu * t
    ) / (v * sqrt_t)

    return (
        np.exp(lam * t)
        * s**gamma
        * (
            phi2(-e1, -f1, rho)

            - (i2 / s)**kappa
            * phi2(-e2, -f2, rho)

            - (i1 / s)**kappa
            * phi2(-e3, -f3, -rho)

            + (i1 / i2)**kappa
            * phi2(-e4, -f4, -rho)
        )
    )

########################################################################################


@njit(cache=True, fastmath=True)
def bjerksund_stensland_2002_value(
    s,
    t,
    k,
    r,
    q,
    v,
    opt_type_value,
):
    """Bjerksund-Stensland (2002) American option approximation.

    Assumes:
        s > 0
        k > 0
        t >= 0
        v > 0
        r >= 0
        q >= 0

    Puts are handled using the American put-call transformation.
    """

    is_call = opt_type_value == OptionTypes.AMERICAN_CALL.value
    is_put = opt_type_value == OptionTypes.AMERICAN_PUT.value

    if not is_call and not is_put:
        return np.nan

    # Expiry
    if t <= G_SMALL:
        if is_call:
            return max(s - k, 0.0)
        return max(k - s, 0.0)

    if s <= 0.0 or k <= 0.0 or v <= 0.0:
        return np.nan

    # Keep the same restricted domain as the BS93 implementation.
    if r < 0.0 or q < 0.0:
        return np.nan

    # ------------------------------------------------------------
    # American put-call transformation
    #
    # P(S,K,r,q) = C(K,S,q,r)
    # ------------------------------------------------------------

    if is_put:
        old_s = s
        old_k = k
        old_r = r
        old_q = q

        s = old_k
        k = old_s
        r = old_q
        q = old_r

    b = r - q

    # No early exercise of the transformed call when q <= 0
    # under the restricted non-negative-rate setup.
    if q <= 0.0:
        return european_value(
            s,
            t,
            k,
            r,
            q,
            v,
            OptionTypes.EUROPEAN_CALL.value,
        )

    v2 = v * v

    # ------------------------------------------------------------
    # beta
    # ------------------------------------------------------------

    beta = (
        0.5
        - b / v2
        + np.sqrt(
            (b / v2 - 0.5) ** 2
            + 2.0 * r / v2
        )
    )

    if not np.isfinite(beta) or beta <= 1.0:
        return np.nan

    # ------------------------------------------------------------
    # Asymptotic boundaries
    # ------------------------------------------------------------

    b_infinity = beta * k / (beta - 1.0)

    # r - b == q
    b0 = max(k, r * k / q)

    denominator = b_infinity - b0

    if denominator <= 0.0 or not np.isfinite(denominator):
        return np.nan

    # ------------------------------------------------------------
    # BS2002 time split
    # ------------------------------------------------------------

    t1 = 0.5 * (np.sqrt(5.0) - 1.0) * t

    sqrt_t1 = np.sqrt(t1)
    sqrt_t = np.sqrt(t)

    # ------------------------------------------------------------
    # Two exercise boundaries
    # ------------------------------------------------------------

    h1 = -(
        b * t1 + 2.0 * v * sqrt_t1
    ) * (k * k) / (
        denominator * b0
    )

    h2 = -(
        b * t + 2.0 * v * sqrt_t
    ) * (k * k) / (
        denominator * b0
    )

    i1 = b0 + denominator * (-np.expm1(h1))
    i2 = b0 + denominator * (-np.expm1(h2))

    if (
        not np.isfinite(i1)
        or not np.isfinite(i2)
        or i1 <= k
        or i2 <= k
    ):
        return np.nan

    # Immediate exercise region
    if s >= i2:
        return s - k

    # ------------------------------------------------------------
    # Alpha coefficients
    # ------------------------------------------------------------

    alpha1 = (i1 - k) * i1**(-beta)
    alpha2 = (i2 - k) * i2**(-beta)

    # ------------------------------------------------------------
    # Bjerksund-Stensland 2002 formula
    # ------------------------------------------------------------

    value = (
        alpha2 * s**beta

        - alpha2
        * _bs02_phi(
            s, t1, beta, i2, i2, r, b, v
        )

        + _bs02_phi(
            s, t1, 1.0, i2, i2, r, b, v
        )

        - _bs02_phi(
            s, t1, 1.0, i1, i2, r, b, v
        )

        - k
        * _bs02_phi(
            s, t1, 0.0, i2, i2, r, b, v
        )

        + k
        * _bs02_phi(
            s, t1, 0.0, i1, i2, r, b, v
        )

        + alpha1
        * _bs02_phi(
            s, t1, beta, i1, i2, r, b, v
        )

        - alpha1
        * _bs02_psi(
            s, t, beta,
            i1, i2, i1, t1,
            r, b, v,
        )

        + _bs02_psi(
            s, t, 1.0,
            i1, i2, i1, t1,
            r, b, v,
        )

        - _bs02_psi(
            s, t, 1.0,
            k, i2, i1, t1,
            r, b, v,
        )

        - k
        * _bs02_psi(
            s, t, 0.0,
            i1, i2, i1, t1,
            r, b, v,
        )

        + k
        * _bs02_psi(
            s, t, 0.0,
            k, i2, i1, t1,
            r, b, v,
        )
    )

    european = european_value(
        s,
        t,
        k,
        r,
        q,
        v,
        OptionTypes.EUROPEAN_CALL.value,
    )

    intrinsic = max(s - k, 0.0)

    return max(
        value,
        european,
        intrinsic,
    )
