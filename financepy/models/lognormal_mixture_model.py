# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 17:59:06 2026

@author: Dominic
"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import least_squares, brentq


class LognormalMixtureModel:
    """
    Two-component lognormal mixture model for a single maturity.

    Parameters
    ----------
    F : float
        Market forward price.
    T : float
        Time to maturity.
    r : float
        Continuously compounded risk-free rate.

    Model parameters
    ----------------
    p : float
        Weight of first component.
    displacement : float
        Determines F1 through F1 = F * exp(displacement).
    sigma1 : float
        Volatility of first component.
    sigma2 : float
        Volatility of second component.

    The second component forward F2 is determined by

        p F1 + (1-p) F2 = F

    so that the mixture satisfies the risk-neutral forward condition.
    """

    def __init__(self, F, T, r=0.0):
        self.F = float(F)
        self.T = float(T)
        self.r = float(r)

        self.p = None
        self.displacement = None
        self.sigma1 = None
        self.sigma2 = None

        self.F1 = None
        self.F2 = None

        self.calibration_result = None

    # ============================================================
    # Black / Black-Scholes utilities
    # ============================================================

    @staticmethod
    def black_call(F, K, T, r, sigma):
        """
        Black call price in forward form:

            C = exp(-rT) [F N(d1) - K N(d2)]
        """
        K = np.asarray(K, dtype=float)

        if T <= 0:
            return np.maximum(F - K, 0.0)

        if sigma <= 0:
            return np.exp(-r * T) * np.maximum(F - K, 0.0)

        vol_sqrt_T = sigma * np.sqrt(T)

        d1 = (np.log(F / K) + 0.5 * sigma**2 * T) / vol_sqrt_T

        d2 = d1 - vol_sqrt_T

        return np.exp(-r * T) * (F * norm.cdf(d1) - K * norm.cdf(d2))

    @staticmethod
    def black_implied_vol(price, F, K, T, r):
        """
        Recover Black implied volatility from a European call price.
        """
        discount = np.exp(-r * T)

        intrinsic = discount * max(F - K, 0.0)
        upper_bound = discount * F

        if price <= intrinsic + 1e-12:
            return 0.0

        if price >= upper_bound:
            return np.nan

        def objective(sigma):
            return LognormalMixtureModel.black_call(F, K, T, r, sigma) - price

        try:
            return brentq(objective, 1e-8, 5.0, xtol=1e-12, rtol=1e-10)
        except ValueError:
            return np.nan

    # ============================================================
    # Parameterisation
    # ============================================================

    def _component_forwards(self, p, displacement):
        """
        Construct component forwards while enforcing

            p F1 + (1-p) F2 = F.
        """
        F1 = self.F * np.exp(displacement)
        F2 = (self.F - p * F1) / (1.0 - p)
        return F1, F2

    @staticmethod
    def _sigmoid(x):
        return 1.0 / (1.0 + np.exp(-x))

    @classmethod
    def _transform_parameters(cls, raw_params):
        """
        Map unconstrained optimisation parameters to valid model parameters.

        p        in (0, 1)
        sigma1   > 0
        sigma2   > 0
        """
        p_raw, displacement, sigma1_raw, sigma2_raw = raw_params

        p = cls._sigmoid(p_raw)
        sigma1 = np.exp(sigma1_raw)
        sigma2 = np.exp(sigma2_raw)

        return p, displacement, sigma1, sigma2

    # ============================================================
    # Pricing
    # ============================================================

    def price(self, K, p=None, displacement=None, sigma1=None, sigma2=None):
        """
        Price European calls under the two-component mixture.
        """
        p = self.p if p is None else p
        displacement = self.displacement if displacement is None else displacement
        sigma1 = self.sigma1 if sigma1 is None else sigma1
        sigma2 = self.sigma2 if sigma2 is None else sigma2

        if any(x is None for x in [p, displacement, sigma1, sigma2]):
            raise ValueError("Model parameters are not set. " "Calibrate the model first or provide parameters.")

        F1, F2 = self._component_forwards(p, displacement)

        if F1 <= 0 or F2 <= 0:
            raise ValueError("Invalid component forward.")

        C1 = self.black_call(F1, K, self.T, self.r, sigma1)

        C2 = self.black_call(F2, K, self.T, self.r, sigma2)

        return p * C1 + (1.0 - p) * C2

    def implied_vol(self, K):
        """
        Return model implied volatility for one or several strikes.
        """
        K = np.asarray(K, dtype=float)

        prices = np.atleast_1d(self.price(K))

        vols = np.array(
            [
                self.black_implied_vol(price, self.F, strike, self.T, self.r)
                for price, strike in zip(prices, np.atleast_1d(K))
            ]
        )

        if K.ndim == 0:
            return vols[0]

        return vols

    # ============================================================
    # Density
    # ============================================================

    @staticmethod
    def _lognormal_density(ST, F_component, T, sigma):
        """
        Lognormal terminal density satisfying

            E[S_T] = F_component.
        """
        ST = np.asarray(ST, dtype=float)

        density = np.zeros_like(ST)

        positive = ST > 0
        x = ST[positive]

        mu = np.log(F_component) - 0.5 * sigma**2 * T

        std = sigma * np.sqrt(T)

        density[positive] = 1.0 / (x * std * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * ((np.log(x) - mu) / std) ** 2)

        return density

    def density(self, ST):
        """
        Evaluate the calibrated terminal risk-neutral density.
        """
        if self.p is None:
            raise ValueError("Model has not been calibrated.")

        g1 = self._lognormal_density(ST, self.F1, self.T, self.sigma1)

        g2 = self._lognormal_density(ST, self.F2, self.T, self.sigma2)

        return self.p * g1 + (1.0 - self.p) * g2

    # ============================================================
    # Calibration
    # ============================================================

    def _residuals(self, raw_params, strikes, market_vols, weights):
        p, displacement, sigma1, sigma2 = self._transform_parameters(raw_params)

        F1, F2 = self._component_forwards(p, displacement)

        # Invalid forward configuration
        if F1 <= 0 or F2 <= 0:
            return np.ones_like(market_vols) * 1e3

        try:
            prices = self.price(strikes, p=p, displacement=displacement, sigma1=sigma1, sigma2=sigma2)
        except ValueError:
            return np.ones_like(market_vols) * 1e3

        model_vols = np.array(
            [self.black_implied_vol(price, self.F, K, self.T, self.r) for price, K in zip(prices, strikes)]
        )

        if np.any(~np.isfinite(model_vols)):
            return np.ones_like(market_vols) * 1e3

        residuals = model_vols - market_vols

        if weights is not None:
            residuals = np.sqrt(weights) * residuals

        return residuals

    def calibrate(self, strikes, market_vols, weights=None, initial_guess=None):
        """
        Calibrate the model to market implied volatilities.

        Parameters
        ----------
        strikes : array-like
            Market strikes.
        market_vols : array-like
            Market implied volatilities in decimal form.
        weights : array-like, optional
            Calibration weights.
        initial_guess : array-like, optional
            Initial guess in raw optimisation coordinates.

        Returns
        -------
        dict
            Calibrated parameters.
        """
        strikes = np.asarray(strikes, dtype=float)
        market_vols = np.asarray(market_vols, dtype=float)

        if weights is not None:
            weights = np.asarray(weights, dtype=float)

        if initial_guess is None:
            atm_vol = market_vols[np.argmin(np.abs(strikes - self.F))]

            initial_guess = np.array(
                [
                    0.0,  # p ~ 0.50
                    -0.03,  # displacement
                    np.log(0.75 * atm_vol),  # sigma1
                    np.log(1.50 * atm_vol),  # sigma2
                ]
            )

        result = least_squares(self._residuals, initial_guess, args=(strikes, market_vols, weights), max_nfev=10000)

        p, displacement, sigma1, sigma2 = self._transform_parameters(result.x)

        F1, F2 = self._component_forwards(p, displacement)

        self.p = p
        self.displacement = displacement
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.F1 = F1
        self.F2 = F2

        self.calibration_result = result

        return self.parameters

    # ============================================================
    # Diagnostics
    # ============================================================

    @property
    def parameters(self):
        """
        Return calibrated model parameters.
        """
        if self.p is None:
            return None

        return {
            "p": self.p,
            "1-p": 1.0 - self.p,
            "F1": self.F1,
            "F2": self.F2,
            "sigma1": self.sigma1,
            "sigma2": self.sigma2,
            "displacement": self.displacement,
        }

    @property
    def forward_check(self):
        """
        Verify the risk-neutral forward condition.
        """
        if self.p is None:
            return None

        return self.p * self.F1 + (1.0 - self.p) * self.F2

    def calibration_errors(self, strikes, market_vols):
        """
        Return model IVs and calibration errors.
        """
        strikes = np.asarray(strikes, dtype=float)
        market_vols = np.asarray(market_vols, dtype=float)

        model_vols = self.implied_vol(strikes)

        return {
            "strikes": strikes,
            "market_vols": market_vols,
            "model_vols": model_vols,
            "errors": model_vols - market_vols,
        }

    def __repr__(self):
        if self.p is None:
            return (
                f"LognormalMixtureModel("
                f"F={self.F:.4f}, "
                f"T={self.T:.4f}, "
                f"r={self.r:.4f}, "
                f"calibrated=False)"
            )

        return (
            f"LognormalMixtureModel("
            f"F={self.F:.4f}, "
            f"T={self.T:.4f}, "
            f"p={self.p:.4f}, "
            f"F1={self.F1:.4f}, "
            f"F2={self.F2:.4f}, "
            f"sigma1={self.sigma1:.4f}, "
            f"sigma2={self.sigma2:.4f})"
        )
