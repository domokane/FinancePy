########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
from scipy.interpolate import PchipInterpolator

from ...utils.error import FinError
from ...utils.math import test_monotonicity
from ...models.option_implied_dbn import option_implied_dbn


class EquityVolCurveNew:
    """Manage a volatility smile/skew at a single maturity.

    The smile is interpolated using total variance in log-forward-moneyness:

        F = S * exp((r - q) * T)
        k = log(K / F)
        w(k) = sigma(k)^2 * T

    PCHIP interpolation is used to reduce artificial oscillation and
    curvature relative to a standard cubic spline.
    """

    ####################################################################################

    def __init__(
        self,
        strikes: np.ndarray,
        volatilities: np.ndarray,
        s: float,
        t_exp: float,
        r: float,
        q: float,
    ) -> None:

        strikes = np.asarray(strikes, dtype=float)
        volatilities = np.asarray(volatilities, dtype=float)

        if strikes.ndim != 1 or volatilities.ndim != 1:
            raise FinError(
                "Strikes and volatilities must be one-dimensional."
            )

        if len(strikes) < 2:
            raise FinError(
                "At least two volatility points are required."
            )

        if len(strikes) != len(volatilities):
            raise FinError(
                "Strike and volatility vectors not same length."
            )

        if np.any(strikes <= 0.0):
            raise FinError(
                "Strikes must be positive."
            )

        if test_monotonicity(strikes) is False:
            raise FinError(
                "Strikes must be strictly monotonic."
            )

        if np.any(volatilities <= 0.0):
            raise FinError(
                "Volatilities must be positive."
            )

        if s <= 0.0:
            raise FinError(
                "Spot price must be positive."
            )

        if t_exp <= 0.0:
            raise FinError(
                "Time to expiry must be positive."
            )

        self._strikes = strikes
        self._volatilities = volatilities

        self._s = float(s)
        self._t_exp = float(t_exp)
        self._r = float(r)
        self._q = float(q)

        #
        # Forward price.
        #
        self._forward = (
            self._s
            * np.exp(
                (self._r - self._q)
                * self._t_exp
            )
        )

        #
        # Log-forward-moneyness:
        #
        #     k = log(K / F)
        #
        self._k = np.log(
            self._strikes / self._forward
        )

        #
        # Total variance:
        #
        #     w = sigma^2 T
        #
        self._w = (
            self._volatilities
            * self._volatilities
            * self._t_exp
        )

        #
        # Shape-preserving cubic interpolation.
        #
        # Extrapolation is disabled deliberately. The wings are handled
        # explicitly below.
        #
        self._f = PchipInterpolator(
            self._k,
            self._w,
            extrapolate=False,
        )

    ####################################################################################

    def _total_variance(
        self,
        strike: float | np.ndarray,
    ) -> float | np.ndarray:
        """Return interpolated total variance."""

        scalar_input = np.isscalar(strike)

        strikes = np.asarray(
            strike,
            dtype=float,
        )

        if np.any(strikes <= 0.0):
            raise FinError(
                "Strikes must be positive."
            )

        k = np.log(
            strikes / self._forward
        )

        w = np.empty_like(
            k,
            dtype=float,
        )

        k_min = self._k[0]
        k_max = self._k[-1]

        inside = (
            (k >= k_min)
            & (k <= k_max)
        )

        left = k < k_min
        right = k > k_max

        #
        # Interpolation inside the quoted range.
        #
        w[inside] = self._f(
            k[inside]
        )

        #
        # Conservative wing treatment.
        #
        # Flat total variance is preferable to uncontrolled cubic
        # extrapolation when the objective is a stable implied density.
        #
        w[left] = self._w[0]
        w[right] = self._w[-1]

        if np.any(w <= 0.0):
            raise FinError(
                "Interpolated total variance is non-positive."
            )

        if scalar_input:
            return float(w)

        return w

    ####################################################################################

    def volatility(
        self,
        strike: float | np.ndarray,
    ) -> float | np.ndarray:
        """Return interpolated volatility for the supplied strike."""

        w = self._total_variance(
            strike
        )

        return np.sqrt(
            w / self._t_exp
        )

    ####################################################################################

    def calculate_pdf(
        self,
        smin: float,
        smax: float,
        n_intervals: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the smile-implied risk-neutral distribution.

        Uses the spot, expiry, rates and dividend yield supplied when the
        volatility curve was constructed.
        """

        if smin <= 0.0:
            raise FinError(
                "Minimum strike must be positive."
            )

        if smax <= smin:
            raise FinError(
                "Maximum strike must be greater than minimum strike."
            )

        if n_intervals < 2:
            raise FinError(
                "Number of intervals must be at least 2."
            )

        strikes = np.linspace(
            smin,
            smax,
            n_intervals + 1,
        )

        sigmas = self.volatility(
            strikes
        )

        density_dk = option_implied_dbn(
            self._s,
            self._t_exp,
            self._r,
            self._q,
            strikes,
            sigmas,
        )

        tol = 1.0e-5

        pmin = np.min(
            density_dk
        )

        if pmin < -tol:
            print(
                "WARNING: "
                "Volatility curve implies negative probability density. "
                f"Minimum value = {pmin:.8g}"
            )

        density_dk = np.maximum(
            density_dk,
            0.0,
        )

        return strikes, density_dk

    ####################################################################################

    @property
    def forward(self) -> float:
        """Return the forward price used by the smile."""

        return self._forward############################################################