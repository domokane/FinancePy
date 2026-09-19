# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 19:29:25 2026

@author: Dominic
"""

import numpy as np
from .lognormal_mixture_model import LognormalMixtureModel


class LognormalMixtureSurface:
    """
    Calibrates one LognormalMixtureModel per maturity.

    Parameters
    ----------
    maturities : array-like
        Times to maturity.
    forwards : array-like
        Forward corresponding to each maturity.
    rates : array-like or float
        Continuously compounded rates.
    strikes : array-like
        Common strike grid across maturities.
    vol_surface : 2D array-like
        Implied volatility surface with shape

            (n_maturities, n_strikes)
    """

    def __init__(self, maturities, forwards, strikes, vol_surface, rates=0.0):
        self.maturities = np.asarray(maturities, dtype=float)

        self.forwards = np.asarray(forwards, dtype=float)

        self.strikes = np.asarray(strikes, dtype=float)

        self.vol_surface = np.asarray(vol_surface, dtype=float)

        if np.isscalar(rates):
            self.rates = np.full(len(self.maturities), rates, dtype=float)
        else:
            self.rates = np.asarray(rates, dtype=float)

        self._validate_inputs()

        self.models = []
        self.calibrated = False

    # ============================================================
    # Validation
    # ============================================================

    def _validate_inputs(self):

        n_t = len(self.maturities)
        n_k = len(self.strikes)

        if len(self.forwards) != n_t:
            raise ValueError("forwards must have one value per maturity.")

        if len(self.rates) != n_t:
            raise ValueError("rates must have one value per maturity.")

        if self.vol_surface.shape != (n_t, n_k):
            raise ValueError("vol_surface must have shape " "(n_maturities, n_strikes).")

        if np.any(self.maturities <= 0.0):
            raise ValueError("Maturities must be positive.")

        if np.any(self.forwards <= 0.0):
            raise ValueError("Forwards must be positive.")

        if np.any(self.vol_surface <= 0.0):
            raise ValueError("Implied volatilities must be positive.")

    # ============================================================
    # Calibration
    # ============================================================

    def calibrate(self, weights=None, warm_start=True):
        """
        Calibrate one mixture model for each maturity.

        Parameters
        ----------
        weights : None or array-like
            Optional calibration weights.

            Can be:
                - 1D array of length n_strikes
                - 2D array matching vol_surface

        warm_start : bool
            If True, use previous maturity's fitted parameters
            as the starting point for the next maturity.
        """

        self.models = []

        previous_raw_guess = None

        for i, T in enumerate(self.maturities):

            model = LognormalMixtureModel(F=self.forwards[i], T=T, r=self.rates[i])

            vols = self.vol_surface[i]

            if weights is None:
                slice_weights = None

            else:
                weights_array = np.asarray(weights, dtype=float)

                if weights_array.ndim == 1:
                    slice_weights = weights_array

                elif weights_array.ndim == 2:
                    slice_weights = weights_array[i]

                else:
                    raise ValueError("weights must be 1D or 2D.")

            if warm_start and previous_raw_guess is not None:

                model.calibrate(
                    strikes=self.strikes, market_vols=vols, weights=slice_weights, initial_guess=previous_raw_guess
                )

            else:

                model.calibrate(strikes=self.strikes, market_vols=vols, weights=slice_weights)

            self.models.append(model)

            if warm_start:
                previous_raw_guess = model.calibration_result.x.copy()

        self.calibrated = True

        return self.models

    # ============================================================
    # Surface reconstruction
    # ============================================================

    def fitted_vol_surface(self):
        """
        Return calibrated implied vols on the original strike grid.
        """

        self._check_calibrated()

        fitted = np.zeros_like(self.vol_surface)

        for i, model in enumerate(self.models):
            fitted[i] = model.implied_vol(self.strikes)

        return fitted

    # ============================================================
    # Pricing
    # ============================================================

    def call_price_surface(self):
        """
        Return call prices on the original strike grid.
        """

        self._check_calibrated()

        prices = np.zeros_like(self.vol_surface)

        for i, model in enumerate(self.models):
            prices[i] = model.price(self.strikes)

        return prices

    # ============================================================
    # Interpolation in maturity
    # ============================================================

    def implied_vol(self, K, T):
        """
        Interpolate implied volatility between calibrated maturities.

        Interpolation is done in total variance:

            w(K,T) = sigma(K,T)^2 T

        which is preferable to direct interpolation of volatility.
        """

        self._check_calibrated()

        K = float(K)
        T = float(T)

        if T < self.maturities[0]:
            raise ValueError("Requested maturity is below calibrated range.")

        if T > self.maturities[-1]:
            raise ValueError("Requested maturity is above calibrated range.")

        # Exact maturity
        idx = np.where(np.isclose(self.maturities, T))[0]

        if len(idx) > 0:
            return self.models[idx[0]].implied_vol(K)

        upper = np.searchsorted(self.maturities, T)

        lower = upper - 1

        T1 = self.maturities[lower]
        T2 = self.maturities[upper]

        sigma1 = self.models[lower].implied_vol(K)

        sigma2 = self.models[upper].implied_vol(K)

        w1 = sigma1**2 * T1
        w2 = sigma2**2 * T2

        alpha = (T - T1) / (T2 - T1)

        w = (1.0 - alpha) * w1 + alpha * w2

        return np.sqrt(w / T)

    # ============================================================
    # Diagnostics
    # ============================================================

    def calibration_errors(self):
        """
        Return model vol errors on the original grid.
        """

        fitted = self.fitted_vol_surface()

        return fitted - self.vol_surface

    def rmse_by_maturity(self):
        """
        RMSE of fitted implied vols by maturity.
        """

        errors = self.calibration_errors()

        return np.sqrt(np.mean(errors**2, axis=1))

    def max_abs_error_by_maturity(self):

        errors = self.calibration_errors()

        return np.max(np.abs(errors), axis=1)

    # ============================================================
    # Calendar arbitrage diagnostic
    # ============================================================

    def calendar_arbitrage_flags(self):
        """
        Simple diagnostic using total variance on the original strike grid.

        Returns True where total variance decreases between maturities.
        """

        fitted = self.fitted_vol_surface()

        total_variance = fitted**2 * self.maturities[:, None]

        differences = np.diff(total_variance, axis=0)

        return differences < -1e-12

    # ============================================================
    # Helpers
    # ============================================================

    def _check_calibrated(self):

        if not self.calibrated:
            raise ValueError("Surface has not been calibrated.")

    def __repr__(self):

        return (
            f"LognormalMixtureSurface("
            f"n_maturities={len(self.maturities)}, "
            f"n_strikes={len(self.strikes)}, "
            f"calibrated={self.calibrated})"
        )
