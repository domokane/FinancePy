from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from ..utils.error import FinError


def normal_sabr_vol(forward, strike, expiry, alpha, rho, nu):
    """Hagan's beta=0 SABR normal (Bachelier) implied-volatility expansion.

    The beta=0 model permits negative forwards and strikes. Formula is
    alpha * z/x(z) * (1 + (2 - 3*rho**2)*nu**2*T/24), with
    z = nu*(forward-strike)/alpha.
    """
    if not np.isfinite([forward, strike, expiry, alpha, rho, nu]).all():
        raise FinError("Normal SABR parameters must be finite.")
    if expiry < 0.0 or alpha <= 0.0 or not -1.0 < rho < 1.0 or nu < 0.0:
        raise FinError("Require T>=0, alpha>0, -1<rho<1, nu>=0.")

    z = nu * (forward - strike) / alpha
    if abs(z) < 1e-7:
        ratio = 1.0 - rho * z / 2.0
    else:
        # asinh formulation of x(z); stable for negative forwards.
        root = np.sqrt(1.0 - rho * rho)
        x = np.arcsinh((z - rho) / root) - np.arcsinh(-rho / root)
        ratio = z / x
    correction = 1.0 + (2.0 - 3.0 * rho * rho) * nu * nu * expiry / 24.0
    return alpha * ratio * correction


class NormalSABR:
    """Normal-volatility SABR with beta fixed at zero."""

    def __init__(self, alpha: float, rho: float, nu: float):
        normal_sabr_vol(0.0, 0.0, 0.0, alpha, rho, nu)
        self.alpha, self.rho, self.nu = float(alpha), float(rho), float(nu)

    def normal_vol(self, forward: float, strike: float, expiry: float) -> float:
        return normal_sabr_vol(
            forward, strike, expiry, self.alpha, self.rho, self.nu
        )

    @classmethod
    def calibrate(cls, forward, strikes, normal_vols, expiry,
                  initial_guess=None):
        """Fit alpha, rho, nu to at least three observed normal volatilities."""
        k = np.asarray(strikes, dtype=float)
        vol = np.asarray(normal_vols, dtype=float)
        if k.ndim != 1 or vol.shape != k.shape or k.size < 3:
            raise FinError("At least three strike/normal-vol observations required.")
        if (not np.isfinite(k).all() or not np.isfinite(vol).all()
                or np.any(vol <= 0.0)):
            raise FinError("Strike/vol observations must be finite with vol > 0.")
        if not np.isfinite([forward, expiry]).all() or expiry < 0.0:
            raise FinError("Invalid forward or expiry.")

        if initial_guess is None:
            alpha0 = float(vol[np.argmin(np.abs(k - forward))])
            # nu is seeded from the observed smile curvature scale.
            initial_guess = (alpha0, 0.0, 0.3)
        guess = np.asarray(initial_guess, dtype=float)
        if guess.shape != (3,) or not np.isfinite(guess).all():
            raise FinError("initial_guess must contain alpha, rho, nu.")
        if guess[0] <= 0 or not -1 < guess[1] < 1 or guess[2] < 0:
            raise FinError("Invalid Normal SABR initial guess.")

        def residuals(params):
            return np.array([
                normal_sabr_vol(forward, strike, expiry, *params)
                for strike in k
            ]) - vol

        fit = least_squares(
            residuals, guess,
            bounds=([1e-12, -0.999, 0.0], [np.inf, 0.999, np.inf]),
            xtol=1e-12, ftol=1e-12, gtol=1e-12,
            max_nfev=4000,
        )
        if not fit.success:
            raise FinError("Normal SABR calibration did not converge.")
        return cls(*fit.x)

