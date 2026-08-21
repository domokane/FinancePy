########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from enum import Enum

import numpy as np

from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import UnivariateSpline

from ...utils.error import FinError
from ...utils.math import test_monotonicity
from ...models.option_implied_dbn import option_implied_dbn

from ...models.black_scholes_analytic import fwd as bsa_fwd
from ...models.black_scholes_analytic import delta as bsa_delta
from ...models.black_scholes_analytic import vega as bsa_vega

########################################################################################


class EquityVolCurveInterpTypes(Enum):
    """Volatility interp scheme."""

    CUBIC_SPLINE = 0
    PCHIP = 1
    SMOOTHING_SPLINE = 2
    LOG_STRIKE_SMOOTHING_SPLINE = 3


########################################################################################


class EquityVolCurveExtrapTypes(Enum):
    """Volatility extrapolation scheme."""

    INTERPOLATOR = 0
    TANGENT = 1


########################################################################################


class EquityVolCurveDeltaTypes(Enum):
    """Smile dynamics used when calculating option delta."""

    STICKY_STRIKE = 0
    STICKY_MONEYNESS = 1


########################################################################################


class EquityVolCurve:
    """Manage a volatility smile/skew to a single maturity."""

    ####################################################################################

    def __init__(
        self,
        strikes: np.ndarray,
        volatilities: np.ndarray,
        s: float,
        t_exp: float,
        r: float,
        q: float,
        interp: EquityVolCurveInterpTypes = (EquityVolCurveInterpTypes.CUBIC_SPLINE),
        extrap: EquityVolCurveExtrapTypes = (EquityVolCurveExtrapTypes.TANGENT),
        smoothing: float = 1.0e-4,
        derivative_bump: float = 1.0e-4,
    ) -> None:

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        volatilities = np.asarray(
            volatilities,
            dtype=float,
        )

        if strikes.ndim != 1 or volatilities.ndim != 1:
            raise FinError("Strikes and volatilities must be one-dimensional.")

        if len(strikes) < 2:
            raise FinError("At least two volatility points are required.")

        if len(strikes) != len(volatilities):
            raise FinError("Strike and volatility vectors not same length.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if test_monotonicity(strikes) is False:
            raise FinError("Strikes must be strictly monotonic.")

        if np.any(np.diff(strikes) <= 0.0):
            raise FinError("Strikes must be strictly increasing.")

        if np.any(volatilities <= 0.0):
            raise FinError("Volatilities must be positive.")

        if s <= 0.0:
            raise FinError("Spot price must be positive.")

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        if smoothing < 0.0:
            raise FinError("Smoothing parameter must be non-negative.")

        if derivative_bump <= 0.0:
            raise FinError("Derivative bump must be positive.")

        if not np.isfinite(r):
            raise FinError("Interest rate must be finite.")

        if not np.isfinite(q):
            raise FinError("Dividend yield must be finite.")

        if not isinstance(
            interp,
            EquityVolCurveInterpTypes,
        ):
            raise FinError("Invalid volatility interp scheme.")

        if not isinstance(
            extrap,
            EquityVolCurveExtrapTypes,
        ):
            raise FinError("Invalid volatility extrapolation scheme.")

        self._strikes = strikes
        self._volatilities = volatilities

        self._s = float(s)
        self._t_exp = float(t_exp)
        self._r = float(r)
        self._q = float(q)

        self._interp = interp
        self._extrap = extrap
        self._smoothing = float(smoothing)
        self._derivative_bump = float(derivative_bump)

        self._forward = bsa_fwd(s, t_exp, r, q)

        if interp == EquityVolCurveInterpTypes.LOG_STRIKE_SMOOTHING_SPLINE:
            self._x = np.log(self._strikes)
            self._use_log_strike = True
        else:
            self._x = self._strikes.copy()
            self._use_log_strike = False

        self._f = self._build_interpolator()

    ####################################################################################

    def _build_interpolator(
        self,
    ):
        """Construct the selected volatility interpolator."""

        x = self._x
        vols = self._volatilities

        interp = self._interp

        if interp == EquityVolCurveInterpTypes.CUBIC_SPLINE:

            if len(x) < 3:
                raise FinError("Cubic spline requires at least three points.")

            return CubicSpline(x, vols, bc_type="natural", extrapolate=True)

        if interp == EquityVolCurveInterpTypes.PCHIP:

            return PchipInterpolator(
                x,
                vols,
                extrapolate=True,
            )

        if interp == EquityVolCurveInterpTypes.SMOOTHING_SPLINE:

            if len(x) < 4:
                raise FinError("Cubic smoothing spline requires at least four points.")

            scale = np.mean(vols) ** 2
            ss = self._smoothing * len(x) * scale

            return UnivariateSpline(
                x,
                vols,
                k=3,
                s=ss,
                ext=0,
            )

        if interp == EquityVolCurveInterpTypes.LOG_STRIKE_SMOOTHING_SPLINE:

            if len(x) < 4:
                raise FinError("Log-strike smoothing spline requires at least " "four points.")

            scale = np.mean(vols) ** 2
            ss = self._smoothing * len(x) * scale

            return UnivariateSpline(
                x,
                vols,
                k=3,
                s=ss,
                ext=0,
            )

        raise FinError("Unknown volatility interp scheme.")

    ####################################################################################

    def _coordinate(
        self,
        strike: np.ndarray,
    ) -> np.ndarray:
        """Convert strike into the interp coordinate."""

        if self._use_log_strike:
            return np.log(strike)

        return strike

    ####################################################################################

    def _endpoint_slopes(
        self,
    ) -> tuple[float, float]:
        """Return endpoint slopes with respect to interp coordinate."""

        x0 = self._x[0]
        x1 = self._x[-1]

        if self._interp in (
            EquityVolCurveInterpTypes.CUBIC_SPLINE,
            EquityVolCurveInterpTypes.PCHIP,
        ):
            left = self._f(
                x0,
                1,
            )
            right = self._f(
                x1,
                1,
            )
            return float(left), float(right)

        if self._interp in (
            EquityVolCurveInterpTypes.SMOOTHING_SPLINE,
            EquityVolCurveInterpTypes.LOG_STRIKE_SMOOTHING_SPLINE,
        ):

            derivative = self._f.derivative(1)
            left = derivative(x0)
            right = derivative(x1)
            return float(left), float(right)

        raise FinError("Unable to calculate endpoint slopes.")

    ####################################################################################

    def volatility(
        self,
        strike: float | np.ndarray,
    ) -> float | np.ndarray:
        """Return interpolated volatility for the supplied strike."""

        scalar = np.isscalar(strike)

        strikes = np.atleast_1d(
            np.asarray(
                strike,
                dtype=float,
            )
        )

        if np.any(strikes <= 0.0):
            raise FinError("Strike must be positive.")

        x = self._coordinate(strikes)

        xmin = self._x[0]
        xmax = self._x[-1]

        left = x < xmin
        right = x > xmax
        inside = ~(left | right)

        vol = np.empty_like(
            strikes,
            dtype=float,
        )

        if np.any(inside):

            vol[inside] = self._f(x[inside])

        if self._extrap == EquityVolCurveExtrapTypes.INTERPOLATOR:

            if np.any(left):

                vol[left] = self._f(x[left])

            if np.any(right):

                vol[right] = self._f(x[right])

        elif self._extrap == EquityVolCurveExtrapTypes.TANGENT:

            left_slope, right_slope = self._endpoint_slopes()

            if np.any(left):

                left_vol = float(self._f(xmin))

                vol[left] = left_vol + left_slope * (x[left] - xmin)

            if np.any(right):

                right_vol = float(self._f(xmax))

                vol[right] = right_vol + right_slope * (x[right] - xmax)

        else:

            raise FinError("Unknown volatility extrapolation scheme.")

        if np.any(~np.isfinite(vol)):
            raise FinError("Invalid interpolated volatility.")

        if np.any(vol <= 0.0):
            raise FinError("Non-positive volatility. Not permitted.")

        if scalar:
            return float(vol[0])

        return vol

    ####################################################################################

    def _volatility_derivative(
        self,
        strike: float | np.ndarray,
    ) -> float | np.ndarray:
        """Return d sigma / dK using a central strike bump.

        The finite-difference derivative uses ``volatility()`` directly,
        so it automatically respects the selected interp and
        extrapolation schemes.
        """

        scalar = np.isscalar(strike)
        strikes = np.atleast_1d(np.asarray(strike, dtype=float))

        if np.any(strikes <= 0.0):
            raise FinError("Strike must be positive.")

        dK = self._derivative_bump * strikes
        dK = np.minimum(dK, 0.5 * strikes)

        strike_up = strikes + dK
        strike_down = strikes - dK

        vol_up = self.volatility(strike_up)
        vol_down = self.volatility(strike_down)
        derivative = (vol_up - vol_down) / (2.0 * dK)

        if np.any(~np.isfinite(derivative)):
            raise FinError("Invalid volatility derivative.")

        if scalar:
            return float(derivative[0])

        return derivative

    ####################################################################################

    def delta_sticky_strike(
        self,
        strike: float | np.ndarray,
        option_type,
    ) -> float | np.ndarray:
        """Return delta assuming sticky-strike volatility.
        The result is the standard Black-Scholes delta evaluated
        using the volatility from this curve.
        """

        scalar = np.isscalar(strike)

        strikes = np.atleast_1d(
            np.asarray(
                strike,
                dtype=float,
            )
        )

        if np.any(strikes <= 0.0):
            raise FinError("Strike must be positive.")

        sigmas = self.volatility(strikes)

        deltas = bsa_delta(self._s, self._t_exp, strikes, self._r, self._q, sigmas, option_type.value)

        if scalar:
            return float(deltas[0])

        return deltas

    ####################################################################################

    def delta_sticky_moneyness(
        self,
        strike: float | np.ndarray,
        option_type,
    ) -> float | np.ndarray:
        """Return delta assuming sticky-moneyness volatility.
        The smile is held fixed as a function of K/S. So if sigma(K/S) then
            d sigma / dS = -(K/S) d sigma/dK
        and
            Delta = Delta_BS + Vega * d sigma/dS.
        """

        scalar = np.isscalar(strike)

        strikes = np.atleast_1d(
            np.asarray(
                strike,
                dtype=float,
            )
        )

        if np.any(strikes <= 0.0):
            raise FinError("Strike must be positive.")

        sigmas = self.volatility(strikes)
        dsigma_dk = self._volatility_derivative(strikes)

        delta_bs = bsa_delta(self._s, self._t_exp, strikes, self._r, self._q, sigmas, option_type.value)
        vega_bs = bsa_vega(self._s, self._t_exp, strikes, self._r, self._q, sigmas, option_type.value)

        dsigma_ds = -(strikes / self._s) * dsigma_dk
        delta = delta_bs + vega_bs * dsigma_ds

        if scalar:
            return float(delta[0])

        return delta

    ####################################################################################

    def delta(
        self,
        strike: float | np.ndarray,
        option_type,
        convention: EquityVolCurveDeltaTypes = (EquityVolCurveDeltaTypes.STICKY_STRIKE),
    ) -> float | np.ndarray:
        """Return delta using the selected smile-dynamics convention."""

        if not isinstance(
            convention,
            EquityVolCurveDeltaTypes,
        ):
            raise FinError("Invalid volatility delta convention.")

        if convention == EquityVolCurveDeltaTypes.STICKY_STRIKE:

            return self.delta_sticky_strike(
                strike,
                option_type,
            )

        if convention == EquityVolCurveDeltaTypes.STICKY_MONEYNESS:

            return self.delta_sticky_moneyness(
                strike,
                option_type,
            )

        raise FinError("Unknown volatility delta convention.")

    ####################################################################################

    def pdf(
        self,
        smin: float,
        smax: float,
        n_intervals: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the smile/skew-implied probability distribution.

        Uses the market state stored in this volatility curve.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Strike grid and probability masses f(K) * dK.
        """

        if smin <= 0.0:
            raise FinError("Minimum strike must be positive.")

        if smax <= smin:
            raise FinError("Maximum strike must be greater than minimum strike.")

        if n_intervals < 2:
            raise FinError("Number of intervals must be at least 2.")

        strikes = np.linspace(
            smin,
            smax,
            n_intervals + 1,
        )

        sigmas = self.volatility(strikes)

        density_dk = option_implied_dbn(
            self._s,
            self._t_exp,
            self._r,
            self._q,
            strikes,
            sigmas,
        )

        tol = 1.0e-5

        pmin = np.min(density_dk)

        if pmin < -tol:

            print("WARNING: " "Volatility curve implies negative probability density. " f"Minimum value = {pmin:.8g}")

        density_dk = np.maximum(
            density_dk,
            0.0,
        )

        return strikes, density_dk

    ####################################################################################
