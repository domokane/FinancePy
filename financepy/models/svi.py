import numpy as np

from scipy.optimize import least_squares
from ..utils.error import FinError

########################################################################################


class SVI:
    """
    Raw SVI implied-volatility smile.

    The total implied variance is

        w(k)
        =
        a
        + b [
            rho (k-m)
            + sqrt((k-m)^2 + sigma^2)
        ]

    where

        k = log(K/F).

    The model may either be constructed from known parameters or calibrated
    directly to an implied-volatility smile.
    """

    def __init__(
        self,
        a=None,
        b=None,
        rho=None,
        m=None,
        sigma=None,
    ) -> None:

        self._a = a
        self._b = b
        self._rho = rho
        self._m = m
        self._sigma = sigma

        if a is not None or b is not None or rho is not None or m is not None or sigma is not None:
            self._validate_parameters()

    ####################################################################################

    def calibrate(
        self,
        forward,
        strikes,
        implied_volatilities,
        t_exp,
    ):
        """
        Calibrate raw SVI parameters to a single implied-volatility smile.

        Calibration is performed in total implied variance.
        """

        if forward <= 0.0:
            raise FinError("Forward must be positive.")

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        implied_volatilities = np.asarray(
            implied_volatilities,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if implied_volatilities.ndim != 1:
            raise FinError("Implied volatilities must be one-dimensional.")

        if len(strikes) != len(implied_volatilities):
            raise FinError("Strikes and implied volatilities must have the same length.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if np.any(implied_volatilities <= 0.0):
            raise FinError("Implied volatilities must be positive.")

        k = np.log(strikes / forward)

        market_variance = implied_volatilities * implied_volatilities * t_exp

        ################################################################################
        # OBJECTIVE
        ################################################################################

        def residuals(params):

            a, b, rho, m, sigma = params

            x = k - m

            model_variance = a + b * (rho * x + np.sqrt(x * x + sigma * sigma))

            return model_variance - market_variance

        ################################################################################
        # INITIAL GUESS
        ################################################################################

        atm_index = np.argmin(np.abs(k))

        atm_variance = market_variance[atm_index]

        x0 = np.array(
            [
                0.5 * atm_variance,
                0.10,
                -0.40,
                0.00,
                0.20,
            ]
        )

        lower_bounds = np.array(
            [
                -1.00,
                1.0e-8,
                -0.999,
                -2.00,
                1.0e-6,
            ]
        )

        upper_bounds = np.array(
            [
                1.00,
                5.00,
                0.999,
                2.00,
                5.00,
            ]
        )

        result = least_squares(
            residuals,
            x0,
            bounds=(
                lower_bounds,
                upper_bounds,
            ),
            xtol=1.0e-12,
            ftol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=10000,
        )

        if not result.success:
            raise FinError("SVI calibration failed.")

        (
            self._a,
            self._b,
            self._rho,
            self._m,
            self._sigma,
        ) = result.x

        self._validate_parameters()

        return result.cost

    ####################################################################################

    def total_variance_from_log_moneyness(
        self,
        k,
    ):

        self._check_calibrated()

        x = k - self._m

        return self._a + self._b * (self._rho * x + np.sqrt(x * x + self._sigma * self._sigma))

    ####################################################################################

    def total_variance(
        self,
        forward,
        strike,
    ):

        if forward <= 0.0:
            raise FinError("Forward must be positive.")

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

        k = np.log(strike / forward)

        return self.total_variance_from_log_moneyness(k)

    ####################################################################################

    def implied_volatility(
        self,
        forward,
        strike,
        t_exp,
    ):

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        w = self.total_variance(
            forward,
            strike,
        )

        if w < 0.0:
            return np.nan

        return np.sqrt(w / t_exp)

    ####################################################################################

    def implied_volatility_curve(
        self,
        forward,
        strikes,
        t_exp,
    ):

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        vols = np.empty(
            len(strikes),
            dtype=float,
        )

        for i, strike in enumerate(strikes):

            vols[i] = self.implied_volatility(
                forward,
                strike,
                t_exp,
            )

        return vols

    ####################################################################################

    def _check_calibrated(
        self,
    ):

        if self._a is None or self._b is None or self._rho is None or self._m is None or self._sigma is None:
            raise FinError("SVI model has not been calibrated.")

    ####################################################################################

    def _validate_parameters(
        self,
    ):

        self._check_calibrated()

        if self._b <= 0.0:
            raise FinError("SVI b must be positive.")

        if self._rho <= -1.0 or self._rho >= 1.0:
            raise FinError("SVI rho must lie between -1 and 1.")

        if self._sigma <= 0.0:
            raise FinError("SVI sigma must be positive.")

    ####################################################################################

    def parameters(
        self,
    ):
        self._check_calibrated()

        return np.array(
            [
                self._a,
                self._b,
                self._rho,
                self._m,
                self._sigma,
            ],
            dtype=float,
        )

    ####################################################################################

    def __repr__(
        self,
    ):

        s = "OBJECT TYPE: SVI\n"
        s += "A: " + str(self._a) + "\n"
        s += "B: " + str(self._b) + "\n"
        s += "RHO: " + str(self._rho) + "\n"
        s += "M: " + str(self._m) + "\n"
        s += "SIGMA: " + str(self._sigma) + "\n"
        return s
