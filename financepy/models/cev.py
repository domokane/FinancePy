# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 21:15:09 2026

@author: Dominic
"""

##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp

import numpy as np
from scipy.stats import ncx2


from ..utils.error import FinError
from ..utils.global_types import OptionTypes
from ..models.black_scholes_analytic import (
    european_value,
    implied_volatility,
)

########################################################################################


class CEV:

    def __init__(self, sigma: float, beta: float):

        if sigma <= 0.0:
            raise FinError("Sigma must be positive.")

        if beta <= 0.0:
            raise FinError("Beta must be positive.")

        if beta > 1.0:
            raise FinError("Analytical CEV pricing currently supports beta <= 1 only.")

        self._sigma = sigma
        self._beta = beta

    ####################################################################################

    def value(
        self,
        stock_price,
        t_exp,
        strike,
        option_type,
        interest_rate,
        dividend_yield,
    ):
        """
        Value a European call or put option under the CEV model.
        """

        self._validate_inputs(
            stock_price,
            t_exp,
            strike,
        )

        call_value = self.call_value(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
        )

        if option_type == OptionTypes.EUROPEAN_CALL.value:
            return call_value

        if option_type == OptionTypes.EUROPEAN_PUT.value:
            return call_value - stock_price * exp(-dividend_yield * t_exp) + strike * exp(-interest_rate * t_exp)

        raise FinError("Unsupported option type.")

    ####################################################################################

    def call_value(
        self,
        stock_price,
        t_exp,
        strike,
        interest_rate,
        dividend_yield,
    ):
        """
        Value a European call option under the CEV model.

        For beta < 1 the exact European CEV formula is expressed in
        terms of non-central chi-square distribution functions.

        For beta = 1 the model reduces to Black-Scholes.
        """

        self._validate_inputs(
            stock_price,
            t_exp,
            strike,
        )

        s = stock_price
        t = t_exp
        k = strike
        r = interest_rate
        q = dividend_yield

        sigma = self._sigma
        beta = self._beta

        # ------------------------------------------------------------------
        # Black-Scholes limiting case
        # ------------------------------------------------------------------

        if abs(beta - 1.0) < 1.0e-4:

            return european_value(
                s,
                t,
                k,
                r,
                q,
                sigma,
                OptionTypes.EUROPEAN_CALL.value,
            )

        # ------------------------------------------------------------------
        # CEV closed-form solution
        #
        # dS = (r-q) S dt + sigma S^beta dW
        #
        # Define
        #
        #     p = 2(1-beta)
        #
        # and map the transition distribution to a non-central chi-square.
        # ------------------------------------------------------------------

        p = 2.0 * (1.0 - beta)

        mu = r - q

        z = mu * p * t

        # The expression for kappa contains
        #
        #     mu / (exp(mu * p * t) - 1)
        #
        # which has a removable singularity at mu = 0.
        #
        # expm1 would also be suitable, but the explicit limiting
        # expression makes the mathematics clear.
        if abs(z) < 1.0e-10:

            kappa = 2.0 / (sigma * sigma * p * p * t)

        else:

            kappa = 2.0 * mu / (sigma * sigma * p * np.expm1(z))

        x = kappa * s**p * exp(z)

        y = kappa * k**p

        df_1 = 2.0 + 2.0 / p
        df_2 = 2.0 / p

        # Survival probability is preferable to 1 - CDF in the
        # far tail because it generally has better numerical accuracy.
        p_1 = ncx2.sf(
            2.0 * y,
            df_1,
            2.0 * x,
        )

        p_2 = ncx2.cdf(
            2.0 * x,
            df_2,
            2.0 * y,
        )

        call_value = s * exp(-q * t) * p_1 - k * exp(-r * t) * p_2

        return call_value

    ####################################################################################

    def put_value(
        self,
        stock_price,
        t_exp,
        strike,
        interest_rate,
        dividend_yield,
    ):
        """
        Value a European put option using put-call parity.
        """

        call_value = self.call_value(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
        )

        put_value = call_value - stock_price * exp(-dividend_yield * t_exp) + strike * exp(-interest_rate * t_exp)

        return put_value

    ####################################################################################

    def local_volatility(
        self,
        stock_price,
    ):
        """
        Return the instantaneous percentage local volatility

            sigma_loc(S) = sigma S^(beta - 1).
        """

        stock_price = np.asarray(
            stock_price,
            dtype=float,
        )

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be positive.")

        local_vol = self._sigma * stock_price ** (self._beta - 1.0)

        if local_vol.ndim == 0:
            return float(local_vol)

        return local_vol

    ####################################################################################

    def diffusion_coefficient(
        self,
        stock_price,
    ):
        """
        Return the absolute diffusion coefficient

            sigma S^beta

        appearing in

            dS = (r-q)S dt + sigma S^beta dW.
        """

        stock_price = np.asarray(
            stock_price,
            dtype=float,
        )

        if np.any(stock_price <= 0.0):
            raise FinError("Stock price must be positive.")

        diffusion = self._sigma * stock_price**self._beta

        if diffusion.ndim == 0:
            return float(diffusion)

        return diffusion

    ####################################################################################

    def implied_volatility(
        self,
        stock_price,
        t_exp,
        strike,
        interest_rate,
        dividend_yield,
    ):
        """
        Return the Black-Scholes implied volatility corresponding
        to the CEV European call price.
        """

        price = self.call_value(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
        )

        return implied_volatility(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
            price,
            OptionTypes.EUROPEAN_CALL.value,
        )

    ####################################################################################

    def implied_volatility_curve(
        self,
        stock_price,
        t_exp,
        strikes,
        interest_rate,
        dividend_yield,
    ):
        """
        Return the Black-Scholes implied-volatility curve across strikes.
        """

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive.")

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        implied_vols = np.empty(
            len(strikes),
            dtype=float,
        )

        for i, strike in enumerate(strikes):

            implied_vols[i] = self.implied_volatility(
                stock_price,
                t_exp,
                strike,
                interest_rate,
                dividend_yield,
            )

        return implied_vols

    ####################################################################################

    def implied_volatility_skew(
        self,
        stock_price,
        t_exp,
        strikes,
        interest_rate,
        dividend_yield,
    ):
        """
        Alias for implied_volatility_curve().

        For beta < 1 the basic CEV model generally produces a
        negative implied-volatility skew.
        """

        return self.implied_volatility_curve(
            stock_price,
            t_exp,
            strikes,
            interest_rate,
            dividend_yield,
        )

    ####################################################################################

    def value_mc(
        self,
        stock_price,
        t_exp,
        strike_price,
        option_type,
        interest_rate,
        dividend_yield,
        num_paths=10000,
        num_steps_per_year=252,
        seed=4242,
    ):
        """
        Monte Carlo valuation using a full-truncation Euler scheme.

        This is primarily intended as a numerical cross-check of the
        analytic European pricing formula.
        """

        self._validate_inputs(
            stock_price,
            t_exp,
            strike_price,
        )

        if num_paths <= 0:
            raise FinError("Number of paths must be positive.")

        if num_steps_per_year <= 0:
            raise FinError("Number of steps per year must be positive.")

        if option_type not in (
            OptionTypes.EUROPEAN_CALL.value,
            OptionTypes.EUROPEAN_PUT.value,
        ):
            raise FinError("Unsupported option type.")

        rng = np.random.default_rng(seed)

        num_steps = max(
            1,
            int(round(t_exp * num_steps_per_year)),
        )

        dt = t_exp / num_steps
        sqrt_dt = np.sqrt(dt)

        s = np.full(
            num_paths,
            stock_price,
            dtype=float,
        )

        mu = interest_rate - dividend_yield

        for _ in range(num_steps):

            z = rng.standard_normal(num_paths)

            s_plus = np.maximum(
                s,
                0.0,
            )

            diffusion = self._sigma * s_plus**self._beta

            s = s + mu * s_plus * dt + diffusion * sqrt_dt * z

            # For beta < 1 zero is an attainable boundary for some
            # parameter choices.  Prevent the Euler discretisation
            # from creating negative stock prices.
            s = np.maximum(
                s,
                0.0,
            )

        if option_type == OptionTypes.EUROPEAN_CALL.value:

            payoff = np.maximum(
                s - strike_price,
                0.0,
            )

        else:

            payoff = np.maximum(
                strike_price - s,
                0.0,
            )

        return exp(-interest_rate * t_exp) * np.mean(payoff)

    #####################################################################################

    def implied_volatility_surface(
        self,
        stock_price,
        expiries,
        strikes,
        interest_rate,
        dividend_yield,
    ):
        """
        Return the Black-Scholes implied-volatility surface generated
        by the CEV model.

        Parameters
        ----------
        stock_price : float
            Current stock price.

        expiries : array_like
            Times to expiry in years.

        strikes : array_like
            Strike grid.

        interest_rate : float
            Continuously compounded risk-free rate.

        dividend_yield : float
            Continuously compounded dividend yield.

        Returns
        -------
        vols : ndarray
            Matrix of implied volatilities with shape

                (len(expiries), len(strikes))

            so that

                vols[i, j]

            is the implied volatility for expiries[i] and strikes[j].
        """

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive.")

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if np.any(expiries <= 0.0):
            raise FinError("Expiries must be positive.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        num_expiries = len(expiries)
        num_strikes = len(strikes)

        vols = np.empty(
            (
                num_expiries,
                num_strikes,
            ),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            vols[i, :] = self.implied_volatility_curve(
                stock_price,
                t_exp,
                strikes,
                interest_rate,
                dividend_yield,
            )

        return vols

    ####################################################################################

    def _validate_inputs(
        self,
        stock_price,
        t_exp,
        strike,
    ):

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive.")

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

    ####################################################################################

    def __repr__(self):

        return "CEV(" f"sigma={self._sigma}, " f"beta={self._beta}" ")"


########################################################################################
