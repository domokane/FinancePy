##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Greeks added thanks to Guillaume Lefieux
##############################################################################

# TODO Fix this
from enum import Enum

import numpy as np
from numba import njit, float64

from ..utils.math import n_vect, n_prime_vect
from ..utils.global_vars import gSmall
from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from ..utils.solver_1d import bisection, newton
from ..models.equity_crr_tree import crr_tree_val_avg
###############################################################################
# TODO: Use Numba ?
###############################################################################


class BlackTypes(Enum):
    ANALYTICAL = 1
    CRR_TREE = 2


class Black():
    """ Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. """

    def __init__(self, volatility, implementation_type=BlackTypes.ANALYTICAL, num_steps=0):
        """ Create FinModel black using parameters. """
        self._volatility = volatility
        self._implementation_type = implementation_type
        self._num_steps = num_steps
        self._seed = 0
        self._param1 = 0
        self._param2 = 0

###############################################################################

    def value(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # df RFR to expiry date
              option_type):    # Call or put
        """ Price a derivative using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self._implementation_type == BlackTypes.ANALYTICAL:
                value = black_value(f, t, k, r, v, option_type)
            else:
                raise FinError("Implementation not available for this product")
        elif option_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self._num_steps, t, option_type.value, k)
            value = results['value']
        else:
            raise FinError(
                "Option type must be a European/American Call or Put")
        return value

###############################################################################

    def delta(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # RFR to expiry date
              option_type):    # Call or put
        """ Calculate delta using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self._implementation_type == BlackTypes.ANALYTICAL:
                delta = black_delta(f, t, k, r, v, option_type)
            else:
                raise FinError("Implementation not available for this product")
        elif option_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self._num_steps, t, option_type.value, k)
            delta = results['delta']
        else:
            raise FinError("Option type must be a European Call or Put")
        return delta

###############################################################################

    def gamma(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # RFR to expiry date
              option_type):    # Call or put
        """ Calculate gamma using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self._implementation_type == BlackTypes.ANALYTICAL:
                gamma = black_gamma(f, t, k, r, v, option_type)
            else:
                raise FinError("Implementation not available for this product")
        elif option_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self._num_steps, t, option_type.value, k)
            gamma = results['gamma']
        else:
            raise FinError("Option type must be a European Call or Put")
        return gamma

###############################################################################

    def theta(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # Discount Factor to expiry date
              option_type):    # Call or put
        """ Calculate theta using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self._implementation_type == BlackTypes.ANALYTICAL:
                theta = black_theta(f, t, k, r, v, option_type)
            else:
                raise FinError("Implementation not available for this product")
        elif option_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self._num_steps, t, option_type.value, k)
            theta = results['theta']
        else:
            raise FinError("Option type must be a European Call or Put")
        return theta

###############################################################################

    def vega(self,
             forward_rate,   # Forward rate F
             strike_rate,    # Strike Rate K
             time_to_expiry,  # Time to Expiry (years)
             df,  # df RFR to expiry date
             option_type):    # Call or put
        """ Calculate vega using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            if self._implementation_type == BlackTypes.ANALYTICAL:
                vega = black_vega(f, t, k, r, v, option_type)
            else:
                raise FinError("Implementation not available for this product")
        elif option_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            bump_size = 0.01
            results = crr_tree_val_avg(
                f, 0.0, 0.0, v, self._num_steps, t, option_type.value, k)
            results_volshfit = crr_tree_val_avg(
                f, 0.0, 0.0, v+bump_size, self._num_steps, t, option_type.value, k)
            vega = (results_volshfit['value'] - results['value']) / bump_size
        else:
            raise FinError("Option type must be a European Call or Put")
        return vega

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self._volatility)
        s += label_to_string("IMPLEMENTATION", self._implementation)
        s += label_to_string("NUMSTEPS", self._num_steps)
        return s

###############################################################################


def black_value(fwd, t, k, r, v, option_type):
    """ Price a derivative using Black model. """
    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    if option_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r*t) * (fwd * n_vect(d1) - k * n_vect(d2))
    elif option_type == OptionTypes.EUROPEAN_PUT:
        return np.exp(-r*t) * (k * n_vect(-d2) - fwd * n_vect(-d1))
    else:
        raise FinError("Option type must be a European Call or Put")


def black_delta(fwd, t, k, r, v, option_type):
    """Return delta of a derivative using Black model. """
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if option_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r*t) * n_vect(d1)
    elif option_type == OptionTypes.EUROPEAN_PUT:
        return - np.exp(-r*t) * n_vect(-d1)
    else:
        raise FinError("Option type must be a European Call or Put")


def black_gamma(fwd, t, k, r, v, option_type):
    """Return gamma of a derivative using Black model. """
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        return np.exp(-r*t) * n_prime_vect(d1) / (fwd * v * np.sqrt(t))
    else:
        raise FinError("Option type must be a European Call or Put")


def black_vega(fwd, t, k, r, v, option_type):
    """Return vega of a derivative using Black model. """
    d1, _ = calculate_d1_d2(fwd, t, k, v)
    if option_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        return np.exp(-r*t) * fwd * np.sqrt(t) * n_prime_vect(d1)
    else:
        raise FinError("Option type must be a European Call or Put")


def black_theta(fwd, t, k, r, v, option_type):
    """Return theta of a derivative using Black model. """
    d1, d2 = calculate_d1_d2(fwd, t, k, v)
    if option_type == OptionTypes.EUROPEAN_CALL:
        return np.exp(-r*t) * (-(fwd * v * n_prime_vect(d1)) / (2 * np.sqrt(t))
                               + r * fwd * n_vect(d1) - r * k * n_vect(d2))
    elif option_type == OptionTypes.EUROPEAN_PUT:
        return np.exp(-r*t) * (-(fwd * v * n_prime_vect(d1)) / (2 * np.sqrt(t))
                               - r * fwd * n_vect(-d1) + r * k * n_vect(-d2))
    else:
        raise FinError("Option type must be a European Call or Put")


@njit(float64[:](float64, float64, float64, float64), fastmath=True, cache=True)
def calculate_d1_d2(f, t, k, v):

    t = np.maximum(t, gSmall)
    vol = np.maximum(v, gSmall)
    k = np.maximum(k, gSmall)
    sqrtT = np.sqrt(t)

    if f <= 0.0:
        raise FinError("Forward is zero.")

    if k <= 0.0:
        raise FinError("Strike is zero.")

    d1 = (np.log(f/k) + vol * vol * t / 2.0) / (vol * sqrtT)
    d2 = d1 - vol * sqrtT

    return np.array([d1, d2])

###############################################################################


def implied_volatility(fwd, t, r, k, price, option_type, debug_print=True):
    """ Calculate the Black implied volatility of a European 
    options on futures contracts using Newton with 
    a fallback to bisection. """

    def _fcall(sigma, args):
        """Function to determine ststar for pricing 
        European call options on future contracts. """
        fwd, t, k, r, price = args
        d1, d2 = calculate_d1_d2(fwd, t, k, sigma)
        value = np.exp(-r*t) * (fwd * n_vect(d1) - k * n_vect(d2))
        obj = value - price
        return obj

    def _fcall_vega(sigma, args):
        """ Function to calculate the Vega of European 
        options on future contracts. """
        fwd, t, k, r, _ = args
        sqrtT = np.sqrt(t)
        d1, _ = calculate_d1_d2(fwd, t, k, sigma)
        vega = np.exp(-r*t) * fwd * sqrtT * n_prime_vect(d1)
        return vega
    # convert to call value
    if option_type == OptionTypes.EUROPEAN_CALL:
        pass
    elif option_type == OptionTypes.EUROPEAN_PUT:
        # put-call parity
        price += np.exp(-r*t) * (fwd - k)
    else:
        raise FinError("Option type must be a European Call or Put")

    # For faster convergence, initial point of inflexion is
    # estimated by Brenner and Subrahmanyan (1988) and Feinstein
    # (1988) approximation for at-the-money forward option. See Eq.(3) in
    # https://www.tandfonline.com/doi/abs/10.2469/faj.v44.n5.80?journalCode=ufaj20
    sigma_guess = price / (0.398 * np.sqrt(t) * fwd * np.exp(-r*t))
    if t < gSmall:  # avoid zero division
        sigma_guess = 0.0
    sigma0 = sigma_guess

    # search implied volatility
    args = fwd, t, k, r, price
    tol = 1e-6
    sigma = newton(_fcall, sigma0, _fcall_vega, args, tol=tol)
    if sigma is None:
        sigma = bisection(_f, 1e-4, 10.0, args, xtol=tol)
        if sigma is None:
            method = "Failed"
        else:
            method = "Bisection"
    else:
        method = "Newton"
    if debug_print:
        print("S: %7.2f K: %7.3f T:%5.3f V:%10.7f Sig0: %7.5f NW: %7.5f %10s" % (
            fwd, k, t, price, sigma0*100.0, sigma*100.0, method))
    return sigma
