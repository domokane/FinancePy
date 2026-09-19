# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Greeks added thanks to Guillaume Lefieux

import numpy as np
from numba import njit, float64, types

from ..utils.math import normcdf_vect, normcdf_prime_vect
from ..utils.global_vars import G_SMALL
from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes
from ..utils.global_types import BlackTypes
from ..utils.error import FinError
from ..utils.solver_1d import bisection, newton
from ..models.equity_crr_tree import crr_tree_val_avg



class Black:
    """Black model for European and American options on forwards/futures."""

    def __init__(
        self,
        volatility: float,
        implementation_type: BlackTypes = BlackTypes.ANALYTICAL,
        num_steps: int = 0,
    ) -> None:

        if volatility < 0.0:
            raise FinError("Volatility must be non-negative.")

        if num_steps < 0:
            raise FinError("num_steps must be non-negative.")

        if not isinstance(implementation_type, BlackTypes):
            raise FinError("Invalid Black implementation type.")

        self.volatility = volatility
        self.implementation_type = implementation_type
        self.num_steps = num_steps
        self.seed = 0
        self.param1 = 0
        self.param2 = 0

    def _validate_inputs(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
    ) -> None:

        if forward_rate < 0.0:
            raise FinError("forward_rate must be non-negative.")

        if strike_rate < 0.0:
            raise FinError("strike_rate must be non-negative.")

        if time_to_expiry < 0.0:
            raise FinError("time_to_expiry must be non-negative.")

        if df <= 0.0 or df > 1.0:
            raise FinError("df must be in (0, 1].")

    @staticmethod
    def _rate_from_df(df: float, time_to_expiry: float) -> float:
        if time_to_expiry <= 0.0:
            return 0.0
        return -np.log(df) / time_to_expiry

    def value(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
        opt_type,
    ) -> float:

        self._validate_inputs(forward_rate, strike_rate, time_to_expiry, df)

        f = forward_rate
        k = strike_rate
        t = time_to_expiry
        v = self.volatility
        r = self._rate_from_df(df, t)

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self.implementation_type is not BlackTypes.ANALYTICAL:
                raise FinError("Implementation not available for this product")

            return black_value(f, t, k, r, v, opt_type)

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            if self.implementation_type is not BlackTypes.CRR_TREE:
                raise FinError("Implementation not available for this product")

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )
            return results["value"]

        raise FinError("Option type must be a European/American Call or Put")

    def delta(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
        opt_type,
    ) -> float:

        self._validate_inputs(forward_rate, strike_rate, time_to_expiry, df)

        f = forward_rate
        k = strike_rate
        t = time_to_expiry
        v = self.volatility
        r = self._rate_from_df(df, t)

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self.implementation_type is not BlackTypes.ANALYTICAL:
                raise FinError("Implementation not available for this product")

            return black_delta(f, t, k, r, v, opt_type)

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            if self.implementation_type is not BlackTypes.CRR_TREE:
                raise FinError("Implementation not available for this product")

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )
            return results["delta"]

        raise FinError("Option type must be a European/American Call or Put")

    def gamma(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
        opt_type,
    ) -> float:

        self._validate_inputs(forward_rate, strike_rate, time_to_expiry, df)

        f = forward_rate
        k = strike_rate
        t = time_to_expiry
        v = self.volatility
        r = self._rate_from_df(df, t)

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self.implementation_type is not BlackTypes.ANALYTICAL:
                raise FinError("Implementation not available for this product")

            return black_gamma(f, t, k, r, v, opt_type)

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            if self.implementation_type is not BlackTypes.CRR_TREE:
                raise FinError("Implementation not available for this product")

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )
            return results["gamma"]

        raise FinError("Option type must be a European/American Call or Put")

    def theta(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
        opt_type,
    ) -> float:

        self._validate_inputs(forward_rate, strike_rate, time_to_expiry, df)

        f = forward_rate
        k = strike_rate
        t = time_to_expiry
        v = self.volatility
        r = self._rate_from_df(df, t)

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self.implementation_type is not BlackTypes.ANALYTICAL:
                raise FinError("Implementation not available for this product")

            return black_theta(f, t, k, r, v, opt_type)

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            if self.implementation_type is not BlackTypes.CRR_TREE:
                raise FinError("Implementation not available for this product")

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )
            return results["theta"]

        raise FinError("Option type must be a European/American Call or Put")

    def vega(
        self,
        forward_rate: float,
        strike_rate: float,
        time_to_expiry: float,
        df: float,
        opt_type,
    ) -> float:

        self._validate_inputs(forward_rate, strike_rate, time_to_expiry, df)

        f = forward_rate
        k = strike_rate
        t = time_to_expiry
        v = self.volatility
        r = self._rate_from_df(df, t)

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self.implementation_type is not BlackTypes.ANALYTICAL:
                raise FinError("Implementation not available for this product")

            return black_vega(f, t, k, r, v, opt_type)

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            if self.implementation_type is not BlackTypes.CRR_TREE:
                raise FinError("Implementation not available for this product")

            bump_size = 0.01

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )

            results_volshift = crr_tree_val_avg(
                f, 0.0, 0.0, v + bump_size, self.num_steps, t, opt_type.value, k
            )

            return (results_volshift["value"] - results["value"]) / bump_size

        raise FinError("Option type must be a European/American Call or Put")

    def __repr__(self):

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        s += label_to_string("IMPLEMENTATION", self.implementation_type)
        s += label_to_string("NUMSTEPS", self.num_steps)
        return s


def black_value(
    fwd: float,
    t: float,
    k: float,
    r: float,
    v: float,
    opt_type,
) -> float:

    if t == 0.0 or v == 0.0:
        df = np.exp(-r * t)

        if opt_type is OptionTypes.EUROPEAN_CALL:
            return df * max(fwd - k, 0.0)

        if opt_type is OptionTypes.EUROPEAN_PUT:
            return df * max(k - fwd, 0.0)

        raise FinError("Option type must be a European Call or Put")

    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    df = np.exp(-r * t)

    if opt_type is OptionTypes.EUROPEAN_CALL:
        return df * (fwd * normcdf_vect(d1) - k * normcdf_vect(d2))

    if opt_type is OptionTypes.EUROPEAN_PUT:
        return df * (k * normcdf_vect(-d2) - fwd * normcdf_vect(-d1))

    raise FinError("Option type must be a European Call or Put")


def black_delta(
    fwd: float,
    t: float,
    k: float,
    r: float,
    v: float,
    opt_type,
) -> float:

    df = np.exp(-r * t)

    if t == 0.0 or v == 0.0:
        if opt_type is OptionTypes.EUROPEAN_CALL:
            if fwd > k:
                return df
            if fwd < k:
                return 0.0
            return 0.5 * df

        if opt_type is OptionTypes.EUROPEAN_PUT:
            if fwd < k:
                return -df
            if fwd > k:
                return 0.0
            return -0.5 * df

        raise FinError("Option type must be a European Call or Put")

    d1, _ = calculate_d1_d2(fwd, t, k, v)

    if opt_type is OptionTypes.EUROPEAN_CALL:
        return df * normcdf_vect(d1)

    if opt_type is OptionTypes.EUROPEAN_PUT:
        return -df * normcdf_vect(-d1)

    raise FinError("Option type must be a European Call or Put")


def black_gamma(
    fwd: float,
    t: float,
    k: float,
    r: float,
    v: float,
    opt_type,
) -> float:

    if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        raise FinError("Option type must be a European Call or Put")

    if t == 0.0 or v == 0.0:
        return 0.0

    d1, _ = calculate_d1_d2(fwd, t, k, v)
    return np.exp(-r * t) * normcdf_prime_vect(d1) / (fwd * v * np.sqrt(t))


def black_vega(
    fwd: float,
    t: float,
    k: float,
    r: float,
    v: float,
    opt_type,
) -> float:

    if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        raise FinError("Option type must be a European Call or Put")

    if t == 0.0 or v == 0.0:
        return 0.0

    d1, _ = calculate_d1_d2(fwd, t, k, v)
    return np.exp(-r * t) * fwd * np.sqrt(t) * normcdf_prime_vect(d1)


def black_theta(
    fwd: float,
    t: float,
    k: float,
    r: float,
    v: float,
    opt_type,
) -> float:

    if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        raise FinError("Option type must be a European Call or Put")

    if t == 0.0 or v == 0.0:
        return 0.0

    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    df = np.exp(-r * t)

    if opt_type is OptionTypes.EUROPEAN_CALL:
        return df * (
            -(fwd * v * normcdf_prime_vect(d1)) / (2.0 * np.sqrt(t))
            + r * fwd * normcdf_vect(d1)
            - r * k * normcdf_vect(d2)
        )

    return df * (
        -(fwd * v * normcdf_prime_vect(d1)) / (2.0 * np.sqrt(t))
        - r * fwd * normcdf_vect(-d1)
        + r * k * normcdf_vect(-d2)
    )


@njit(
    types.UniTuple(float64, 2)(float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def calculate_d1_d2(
    f: float,
    t: float,
    k: float,
    v: float,
):

    t = max(t, G_SMALL)
    vol = max(v, G_SMALL)
    strike = max(k, G_SMALL)
    sqrt_t = np.sqrt(t)

    if f <= 0.0:
        raise ValueError("Forward must be positive.")

    if strike <= 0.0:
        raise ValueError("Strike must be positive.")

    d1 = (np.log(f / strike) + 0.5 * vol * vol * t) / (vol * sqrt_t)
    d2 = d1 - vol * sqrt_t

    return d1, d2


def implied_volatility(
    fwd: float,
    t: float,
    r: float,
    k: float,
    price: float,
    opt_type,
    debug_print: bool = False,
) -> float:
    """Calculate Black implied volatility using Newton with bisection fallback."""

    if fwd <= 0.0:
        raise FinError("fwd must be positive.")

    if k <= 0.0:
        raise FinError("strike must be positive.")

    if t <= 0.0:
        raise FinError("time_to_expiry must be positive.")

    if price < 0.0:
        raise FinError("price must be non-negative.")

    def _f_european(sigma, args):
        fwd, t, k, r, opt_type, price = args
        return black_value(fwd, t, k, r, sigma, opt_type) - price

    def _f_european_vega(sigma, args):
        fwd, t, k, r, opt_type, _ = args
        return black_vega(fwd, t, k, r, sigma, opt_type)

    def _f_american(sigma, args):
        fwd, t, k, _, opt_type, price = args
        num_steps = 200
        results = crr_tree_val_avg(
            fwd, 0.0, 0.0, sigma, num_steps, t, opt_type.value, k
        )
        return results["value"] - price

    def _f_american_vega(sigma, args):
        fwd, t, k, _, opt_type, _ = args
        bump_size = 0.01
        num_steps = 200

        results = crr_tree_val_avg(
            fwd, 0.0, 0.0, sigma, num_steps, t, opt_type.value, k
        )

        results_volshift = crr_tree_val_avg(
            fwd, 0.0, 0.0, sigma + bump_size, num_steps, t, opt_type.value, k
        )

        return (results_volshift["value"] - results["value"]) / bump_size

    def _estimate_vol_from_price(fwd, t, k, european_opt_type, european_price):

        if european_opt_type is OptionTypes.EUROPEAN_CALL:
            atm_price = european_price

        elif european_opt_type is OptionTypes.EUROPEAN_PUT:
            atm_price = european_price + np.exp(-r * t) * (fwd - k)

        else:
            raise FinError("Option type must be a European Call or Put")

        denom = 0.398 * np.sqrt(t) * fwd * np.exp(-r * t)

        if denom <= G_SMALL:
            return 0.20

        sigma_guess = atm_price / denom
        return float(np.clip(sigma_guess, 1.0e-4, 5.0))

    if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        _f = _f_european
        _f_vega = _f_european_vega
        sigma0 = _estimate_vol_from_price(fwd, t, k, opt_type, price)

    elif opt_type is OptionTypes.AMERICAN_CALL:
        _f = _f_american
        _f_vega = _f_american_vega
        sigma0 = _estimate_vol_from_price(fwd, t, k, OptionTypes.EUROPEAN_CALL, price)

    elif opt_type is OptionTypes.AMERICAN_PUT:
        _f = _f_american
        _f_vega = _f_american_vega
        sigma0 = _estimate_vol_from_price(fwd, t, k, OptionTypes.EUROPEAN_PUT, price)

    else:
        raise FinError("Option type must be a European/American Call or Put")

    args = fwd, t, k, r, opt_type, price
    tol = 1.0e-6

    sigma = newton(_f, sigma0, _f_vega, args, tol=tol)

    if sigma is None:
        sigma = bisection(_f, 1.0e-4, 10.0, args, xtol=tol)
        method = "Failed" if sigma is None else "Bisection"
    else:
        method = "Newton"

    if debug_print:
        print(
            "S: %7.2f K: %7.3f T:%5.3f V:%10.7f Sig0: %7.5f NW: %7.5f %10s"
            % (fwd, k, t, price, sigma0 * 100.0, sigma * 100.0, method)
        )

    return sigma
