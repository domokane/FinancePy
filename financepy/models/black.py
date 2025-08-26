# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Greeks added thanks to Guillaume Lefieux

# TODO Fix this
from enum import Enum

import numpy as np
from numba import njit, float64

from ..utils.math import normcdf_vect, normcdf_prime_vect
from ..utils.global_vars import G_SMALL
from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from ..utils.solver_1d import bisection, newton
from ..models.equity_crr_tree import crr_tree_val_avg

# TODO: Use Numba ?

########################################################################################


class BlackTypes(Enum):

    ANALYTICAL = 1
    CRR_TREE = 2


########################################################################################


class Black:
    """Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation."""

    ####################################################################################

    def __init__(
        self,
        volatility,
        implementation_type=BlackTypes.ANALYTICAL,
        num_steps=0,
    ):
        """Create FinModel black using parameters."""
        self.volatility = volatility
        self.implementation_type = implementation_type
        self.num_steps = num_steps
        self.seed = 0
        self.param1 = 0
        self.param2 = 0

    ####################################################################################

    def value(
        self,
        forward_rate,  # Forward rate F
        strike_rate,  # Strike Rate K
        time_to_expiry,  # Time to Expiry (years)
        df,  # df RFR to expiry date
        opt_type,
    ):  # Call or put
        """Price a derivative using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self.volatility
        r = -np.log(df) / t
        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            if self.implementation_type == BlackTypes.ANALYTICAL:
                value = black_value(f, t, k, r, v, opt_type)
            else:
                raise FinError("Implementation not available for this product")

        elif opt_type in (
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):

            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
            )
            value = results["value"]
        else:
            raise FinError("Option type must be a European/American Call or Put")
        return value

    ####################################################################################

    def delta(
        self,
        forward_rate,  # Forward rate
        strike_rate,  # Strike Rate
        time_to_expiry,  # Time to Expiry (years)
        df,  # Discount factor to expiry date
        opt_type,
    ):  # Call or put
        """Calculate delta using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self.volatility
        r = -np.log(df) / t

        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):

            if self.implementation_type == BlackTypes.ANALYTICAL:
                return black_delta(f, t, k, r, v, opt_type)

            raise FinError("Implementation not available for this product")

        elif opt_type in (
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):
            if self.implementation_type == BlackTypes.CRR_TREE:
                results = crr_tree_val_avg(
                    f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
                )
                return results["delta"]

            raise FinError("Implementation not available for this product")

        else:
            raise FinError("Option type must be a European/American Call or Put")

    ####################################################################################

    def gamma(
        self,
        forward_rate,  # Forward rate F
        strike_rate,  # Strike Rate K
        time_to_expiry,  # Time to Expiry (years)
        df,  # RFR to expiry date
        opt_type,
    ):  # Call or put
        """Calculate gamma using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self.volatility
        r = -np.log(df) / t

        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            if self.implementation_type == BlackTypes.ANALYTICAL:
                return black_gamma(f, t, k, r, v, opt_type)
            else:
                raise FinError("Implementation not available for this product")
        elif opt_type in (
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):
            if self.implementation_type == BlackTypes.CRR_TREE:
                results = crr_tree_val_avg(
                    f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
                )
                return results["gamma"]
            else:
                raise FinError("Implementation not available for this product")
        else:
            raise FinError("Option type must be a European/American Call or Put")

    ####################################################################################

    def theta(
        self,
        forward_rate,  # Forward rate F
        strike_rate,  # Strike Rate K
        time_to_expiry,  # Time to Expiry (years)
        df,  # Discount Factor to expiry date
        opt_type,
    ):  # Call or put
        """Calculate theta using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self.volatility
        r = -np.log(df) / t

        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            if self.implementation_type == BlackTypes.ANALYTICAL:
                theta = black_theta(f, t, k, r, v, opt_type)
            else:
                raise FinError("Implementation not available for this product")
        elif opt_type in (
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):
            if self.implementation_type == BlackTypes.CRR_TREE:
                results = crr_tree_val_avg(
                    f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
                )
                return results["theta"]
            else:
                raise FinError("Implementation not available for this product")
        else:
            raise FinError("Option type must be a European/American Call or Put")

        return theta

    ####################################################################################

    def vega(
        self,
        forward_rate,  # Forward rate F
        strike_rate,  # Strike Rate K
        time_to_expiry,  # Time to Expiry (years)
        df,  # df RFR to expiry date
        opt_type,
    ):  # Call or put
        """Calculate vega using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self.volatility
        r = -np.log(df) / t

        if opt_type in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            if self.implementation_type == BlackTypes.ANALYTICAL:
                vega = black_vega(f, t, k, r, v, opt_type)
            else:
                raise FinError("Implementation not available for this product")
        elif opt_type in (
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT,
        ):
            if self.implementation_type == BlackTypes.CRR_TREE:
                bump_size = 0.01
                results = crr_tree_val_avg(
                    f, 0.0, 0.0, v, self.num_steps, t, opt_type.value, k
                )
                results_volshift = crr_tree_val_avg(
                    f,
                    0.0,
                    0.0,
                    v + bump_size,
                    self.num_steps,
                    t,
                    opt_type.value,
                    k,
                )
                vega = (results_volshift["value"] - results["value"]) / bump_size
                return vega
            else:
                raise FinError("Implementation not available for this product")
        else:
            raise FinError("Option type must be a European/American Call or Put")
        return vega

    ####################################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        s += label_to_string("IMPLEMENTATION", self.implementation_type)
        s += label_to_string("NUMSTEPS", self.num_steps)
        return s


########################################################################################


def black_value(fwd, t, k, r, v, opt_type):
    """Price a derivative using Black model."""
    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    if opt_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r * t) * (fwd * normcdf_vect(d1) - k * normcdf_vect(d2))
    elif opt_type == OptionTypes.EUROPEAN_PUT:
        return np.exp(-r * t) * (k * normcdf_vect(-d2) - fwd * normcdf_vect(-d1))
    else:
        raise FinError("Option type must be a European Call or Put")


########################################################################################


def black_delta(fwd, t, k, r, v, opt_type):
    """Return delta of a derivative using Black model."""
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if opt_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r * t) * normcdf_vect(d1)
    elif opt_type == OptionTypes.EUROPEAN_PUT:
        return -np.exp(-r * t) * normcdf_vect(-d1)
    else:
        raise FinError("Option type must be a European Call or Put")


########################################################################################


def black_gamma(fwd, t, k, r, v, opt_type):
    """Return gamma of a derivative using Black model."""
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        return np.exp(-r * t) * normcdf_prime_vect(d1) / (fwd * v * np.sqrt(t))
    else:
        raise FinError("Option type must be a European Call or Put")


########################################################################################


def black_vega(fwd, t, k, r, v, opt_type):
    """Return vega of a derivative using Black model."""
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        return np.exp(-r * t) * fwd * np.sqrt(t) * normcdf_prime_vect(d1)
    else:
        raise FinError("Option type must be a European Call or Put")


########################################################################################


def black_theta(fwd, t, k, r, v, opt_type):
    """Return theta of a derivative using Black model."""
    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    if opt_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r * t) * (
            -(fwd * v * normcdf_prime_vect(d1)) / (2 * np.sqrt(t))
            + r * fwd * normcdf_vect(d1)
            - r * k * normcdf_vect(d2)
        )
    elif opt_type == OptionTypes.EUROPEAN_PUT:
        return np.exp(-r * t) * (
            -(fwd * v * normcdf_prime_vect(d1)) / (2 * np.sqrt(t))
            - r * fwd * normcdf_vect(-d1)
            + r * k * normcdf_vect(-d2)
        )
    else:
        raise FinError("Option type must be a European Call or Put")


########################################################################################


@njit(float64[:](float64, float64, float64, float64), fastmath=True, cache=True)
def calculate_d1_d2(f, t, k, v):
    """Calculate d1 and d2 for Black-Scholes model."""

    t = np.maximum(t, G_SMALL)
    vol = np.maximum(v, G_SMALL)
    k = np.maximum(k, G_SMALL)
    sqrt_t = np.sqrt(t)

    if f <= 0.0:
        raise FinError("Forward is zero.")

    if k <= 0.0:
        raise FinError("Strike is zero.")

    d1 = (np.log(f / k) + vol * vol * t / 2.0) / (vol * sqrt_t)
    d2 = d1 - vol * sqrt_t

    return np.array([d1, d2])


########################################################################################


def implied_volatility(fwd, t, r, k, price, opt_type, debug_print=True):
    """Calculate the Black implied volatility of a European/American
    options on futures contracts using Newton with
    a fallback to bisection."""

    ####################################################################################

    def _f_european(sigma, args):
        """Function to determine ststar for pricing
        European options on future contracts."""
        fwd, t, k, r, opt_type, price = args
        value = black_value(fwd, t, k, r, sigma, opt_type)
        obj = value - price
        return obj

    ####################################################################################

    def _f_european_vega(sigma, args):
        """Function to calculate the Vega of European
        options on future contracts."""
        fwd, t, k, r, opt_type, _ = args
        vega = black_vega(fwd, t, k, r, sigma, opt_type)
        return vega

    ####################################################################################

    def _f_american(sigma, args):
        """Function to determine ststar for pricing
        American options on future contracts."""
        fwd, t, k, _, opt_type, price = args
        num_steps = 200
        results = crr_tree_val_avg(
            fwd, 0.0, 0.0, sigma, num_steps, t, opt_type.value, k
        )
        obj = results["value"] - price
        return obj

    ####################################################################################

    def _f_american_vega(sigma, args):
        """Function to calculate the Vega of American
        options on future contracts."""
        fwd, t, k, _, opt_type, _ = args
        bump_size = 0.01
        num_steps = 200
        results = crr_tree_val_avg(
            fwd, 0.0, 0.0, sigma, num_steps, t, opt_type.value, k
        )
        results_volshift = crr_tree_val_avg(
            fwd,
            0.0,
            0.0,
            sigma + bump_size,
            num_steps,
            t,
            opt_type.value,
            k,
        )
        vega = (results_volshift["value"] - results["value"]) / bump_size
        return vega

    ####################################################################################

    def _estimate_vol_from_price(fwd, t, k, european_opt_type, european_price):

        # Brenner and Subrahmanyan (1988) and Feinstein
        # (1988) approximation for at-the-money forward call option. See Eq.(3) in
        # https://www.tandfonline.com/doi/abs/10.2469/faj.v44.n5.80?journalCode=ufaj20
        if european_opt_type == OptionTypes.EUROPEAN_CALL:
            price = european_price
        elif european_opt_type == OptionTypes.EUROPEAN_PUT:
            price = european_price + np.exp(-r * t) * (fwd - k)
        else:
            raise FinError("Option type must be a European Call or Put")
        sigma_guess = price / (0.398 * np.sqrt(t) * fwd * np.exp(-r * t))
        if t < G_SMALL:  # avoid zero division
            sigma_guess = 0.0
        return sigma_guess

    # Set objective function, its first derivative, and initial point of
    # implied volatility
    # A simple approximation is used to estimate implied volatility
    # ,which is used as the input of calibration.
    if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        _f = _f_european
        _f_vega = _f_european_vega
        sigma0 = _estimate_vol_from_price(fwd, t, k, opt_type, price)
    elif opt_type == OptionTypes.AMERICAN_CALL:
        _f = _f_american
        _f_vega = _f_american_vega
        # NOTE:
        # Instead of european price, american price is
        # passed to an approximation formula for european option.
        # But it's just an initial point, so has no affect on the calibration.
        sigma0 = _estimate_vol_from_price(fwd, t, k, OptionTypes.EUROPEAN_CALL, price)
    elif opt_type == OptionTypes.AMERICAN_PUT:
        _f = _f_american
        _f_vega = _f_american_vega
        # NOTE:
        # Same argument as american call's case
        sigma0 = _estimate_vol_from_price(fwd, t, k, OptionTypes.EUROPEAN_PUT, price)
    else:
        raise FinError("Option type must be a European Call or Put")

    # search implied volatility
    args = fwd, t, k, r, opt_type, price
    tol = 1e-6
    sigma = newton(_f, sigma0, _f_vega, args, tol=tol)
    if sigma is None:
        sigma = bisection(_f, 1e-4, 10.0, args, xtol=tol)
        if sigma is None:
            method = "Failed"
        else:
            method = "Bisection"
    else:
        method = "Newton"
    if debug_print:
        print(
            "S: %7.2f K: %7.3f T:%5.3f V:%10.7f Sig0: %7.5f NW: %7.5f %10s"
            % (fwd, k, t, price, sigma0 * 100.0, sigma * 100.0, method)
        )
    return sigma
