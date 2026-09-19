import numpy as np

from scipy.optimize import least_squares

from ..utils.error import FinError
from .implied_volatility_surface import ImpliedVolatilitySurface

########################################################################################


class SSVIPowerLawPhi:
    """
    Power-law SSVI phi function

        phi(theta) = eta * theta^(-gamma)
    """

    def __init__(
        self,
        eta,
        gamma,
    ) -> None:

        if eta <= 0.0:
            raise FinError("Eta must be positive.")

        if gamma < 0.0 or gamma > 1.0:
            raise FinError("Gamma must lie between 0 and 1.")

        self._eta = eta
        self._gamma = gamma

    ####################################################################################

    def value(
        self,
        theta,
    ):

        if theta <= 0.0:
            raise FinError("Theta must be positive.")

        return self._eta * theta ** (-self._gamma)

    ####################################################################################

    def __call__(
        self,
        theta,
    ):

        return self.value(theta)

    ####################################################################################

    def parameters(
        self,
    ):

        return np.array(
            [
                self._eta,
                self._gamma,
            ],
            dtype=float,
        )

    ####################################################################################

    def __repr__(
        self,
    ):

        s = "OBJECT TYPE: SSVI POWER LAW PHI\n"
        s += "ETA: " + str(self._eta) + "\n"
        s += "GAMMA: " + str(self._gamma) + "\n"

        return s


########################################################################################


class SSVISurface(ImpliedVolatilitySurface):
    """
    SSVI implied-volatility surface

        w(k, theta)
        =
        theta / 2
        [
            1
            + rho phi(theta) k
            + sqrt(
                (phi(theta) k + rho)^2
                + 1 - rho^2
            )
        ]

    using the power-law specification

        phi(theta)
        =
        eta theta^(-gamma).

    The surface can either be constructed from known parameters
    or calibrated directly to a market implied-volatility grid.
    """

    def __init__(
        self,
        expiries=None,
        atm_total_variances=None,
        rho=None,
        phi_function=None,
    ) -> None:

        self._expiries = None
        self._atm_total_variances = None
        self._rho = None
        self._phi_function = None

        supplied = [
            expiries is not None,
            atm_total_variances is not None,
            rho is not None,
            phi_function is not None,
        ]

        num_supplied = sum(supplied)

        if num_supplied != 0 and num_supplied != 4:
            raise FinError("Either supply all SSVI parameters or none.")

        if num_supplied == 4:

            self._expiries = np.asarray(
                expiries,
                dtype=float,
            )

            self._atm_total_variances = np.asarray(
                atm_total_variances,
                dtype=float,
            )

            self._rho = float(rho)

            self._phi_function = phi_function

            self._validate_parameters()

    ####################################################################################

    def calibrate(
        self,
        forwards,
        strikes,
        expiries,
        implied_volatilities,
    ):
        """
        Jointly calibrate

            theta(T_1), ..., theta(T_N), rho, eta, gamma

        to the complete implied-volatility grid.

        The theta term structure is parameterized so that

            0 < theta_1 < theta_2 < ... < theta_N.

        This guarantees monotonic ATM total variance.
        """

        forwards = np.asarray(
            forwards,
            dtype=float,
        )

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        implied_volatilities = np.asarray(
            implied_volatilities,
            dtype=float,
        )

        self._validate_market_data(
            forwards,
            strikes,
            expiries,
            implied_volatilities,
        )

        num_expiries = len(expiries)

        ################################################################################
        # INITIAL THETA ESTIMATES
        ################################################################################
        #
        # These are only starting values for the optimizer.
        #
        ################################################################################

        theta0 = np.empty(
            num_expiries,
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            k = np.log(strikes / forwards[i])

            market_variance = implied_volatilities[i] * implied_volatilities[i] * t_exp

            theta0[i] = np.interp(
                0.0,
                k,
                market_variance,
            )

        ################################################################################
        # ENSURE STRICTLY INCREASING INITIAL THETA
        ################################################################################

        minimum_increment = 1.0e-8

        theta0[0] = max(
            theta0[0],
            minimum_increment,
        )

        for i in range(
            1,
            num_expiries,
        ):

            theta0[i] = max(
                theta0[i],
                theta0[i - 1] + minimum_increment,
            )

        ################################################################################
        # PARAMETERIZE THETA USING LOG INCREMENTS
        ################################################################################
        #
        # y[0] gives theta_1:
        #
        #     theta_1 = exp(y[0])
        #
        # and for i > 0:
        #
        #     theta_i = theta_{i-1} + exp(y[i])
        #
        ################################################################################

        theta_parameters0 = np.empty(
            num_expiries,
            dtype=float,
        )

        theta_parameters0[0] = np.log(theta0[0])

        for i in range(
            1,
            num_expiries,
        ):

            theta_parameters0[i] = np.log(theta0[i] - theta0[i - 1])

        ################################################################################
        # FULL INITIAL VECTOR
        ################################################################################

        x0 = np.concatenate(
            (
                theta_parameters0,
                np.array(
                    [
                        -0.50,  # rho
                        1.00,  # eta
                        0.30,  # gamma
                    ],
                    dtype=float,
                ),
            )
        )

        ################################################################################
        # PARAMETER DECODER
        ################################################################################

        def unpack_parameters(
            params,
        ):

            theta_parameters = params[:num_expiries]

            rho = params[num_expiries]

            eta = params[num_expiries + 1]

            gamma = params[num_expiries + 2]

            theta = np.empty(
                num_expiries,
                dtype=float,
            )

            theta[0] = np.exp(theta_parameters[0])

            for i in range(
                1,
                num_expiries,
            ):

                theta[i] = theta[i - 1] + np.exp(theta_parameters[i])

            return (
                theta,
                rho,
                eta,
                gamma,
            )

        ################################################################################
        # OBJECTIVE FUNCTION
        ################################################################################

        def residuals(
            params,
        ):

            (
                theta,
                rho,
                eta,
                gamma,
            ) = unpack_parameters(params)

            errors = []

            for i, t_exp in enumerate(expiries):

                k = np.log(strikes / forwards[i])

                phi = eta * theta[i] ** (-gamma)

                x = phi * k

                model_variance = 0.5 * theta[i] * (1.0 + rho * x + np.sqrt((x + rho) ** 2 + 1.0 - rho * rho))

                model_volatility = np.sqrt(model_variance / t_exp)

                errors.extend(model_volatility - implied_volatilities[i])

            return np.asarray(
                errors,
                dtype=float,
            )

        def residuals_old(
            params,
        ):

            (
                theta,
                rho,
                eta,
                gamma,
            ) = unpack_parameters(params)

            errors = []

            for i, t_exp in enumerate(expiries):

                k = np.log(strikes / forwards[i])

                market_variance = implied_volatilities[i] * implied_volatilities[i] * t_exp

                phi = eta * theta[i] ** (-gamma)

                x = phi * k

                model_variance = 0.5 * theta[i] * (1.0 + rho * x + np.sqrt((x + rho) ** 2 + 1.0 - rho * rho))

                errors.extend(model_variance - market_variance)

            return np.asarray(
                errors,
                dtype=float,
            )

        ################################################################################
        # BOUNDS
        ################################################################################
        #
        # The theta parameters are unconstrained because exponentiation
        # guarantees positive increments.
        #
        ################################################################################

        lower_bounds = np.concatenate(
            (
                np.full(
                    num_expiries,
                    -np.inf,
                ),
                np.array(
                    [
                        -0.999,
                        1.0e-8,
                        0.0,
                    ]
                ),
            )
        )

        upper_bounds = np.concatenate(
            (
                np.full(
                    num_expiries,
                    np.inf,
                ),
                np.array(
                    [
                        0.999,
                        10.0,
                        1.0,
                    ]
                ),
            )
        )

        ################################################################################
        # CALIBRATION
        ################################################################################

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
            max_nfev=20000,
        )

        if not result.success:

            raise FinError("SSVI calibration failed.")

        (
            theta,
            rho,
            eta,
            gamma,
        ) = unpack_parameters(result.x)

        self._expiries = expiries.copy()

        self._atm_total_variances = theta

        self._rho = float(rho)

        self._phi_function = SSVIPowerLawPhi(
            eta,
            gamma,
        )

        self._validate_parameters()

        return result.cost

    ####################################################################################

    def theta(
        self,
        t_exp,
    ):
        """
        Return ATM total variance theta(T).

        Linear interpolation is performed in total variance.
        """

        self._check_calibrated()

        if t_exp < self._expiries[0]:

            raise FinError("Expiry is below SSVI surface.")

        if t_exp > self._expiries[-1]:

            raise FinError("Expiry is above SSVI surface.")

        return np.interp(
            t_exp,
            self._expiries,
            self._atm_total_variances,
        )

    ####################################################################################

    def total_variance_from_log_moneyness(
        self,
        k,
        t_exp,
    ):
        """
        Return SSVI total implied variance.
        """

        self._check_calibrated()

        theta = self.theta(t_exp)

        phi = self._phi_function(theta)

        x = phi * k

        return 0.5 * theta * (1.0 + self._rho * x + np.sqrt((x + self._rho) ** 2 + 1.0 - self._rho * self._rho))

    ####################################################################################

    def total_variance(
        self,
        forward,
        strike,
        t_exp,
    ):

        if forward <= 0.0:

            raise FinError("Forward must be positive.")

        if strike <= 0.0:

            raise FinError("Strike must be positive.")

        if t_exp <= 0.0:

            raise FinError("Expiry must be positive.")

        k = np.log(strike / forward)

        return self.total_variance_from_log_moneyness(
            k,
            t_exp,
        )

    ####################################################################################

    def implied_volatility(
        self,
        forward,
        strike,
        t_exp,
    ):

        w = self.total_variance(
            forward,
            strike,
            t_exp,
        )

        if not np.isfinite(w) or w < 0.0:

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

    def implied_volatility_surface(
        self,
        forwards,
        strikes,
        expiries,
    ):

        forwards = np.asarray(
            forwards,
            dtype=float,
        )

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        if len(forwards) != len(expiries):

            raise FinError("One forward is required for each expiry.")

        vols = np.empty(
            (
                len(expiries),
                len(strikes),
            ),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            vols[i, :] = self.implied_volatility_curve(
                forwards[i],
                strikes,
                t_exp,
            )

        return vols

    ####################################################################################

    def parameters(
        self,
    ):
        """
        Return global SSVI parameters

            rho, eta, gamma.
        """

        self._check_calibrated()

        eta, gamma = self._phi_function.parameters()

        return np.array(
            [
                self._rho,
                eta,
                gamma,
            ],
            dtype=float,
        )

    ####################################################################################

    def atm_total_variances(
        self,
    ):

        self._check_calibrated()

        return self._atm_total_variances.copy()

    ####################################################################################

    def expiries(
        self,
    ):

        self._check_calibrated()

        return self._expiries.copy()

    ####################################################################################

    def _validate_market_data(
        self,
        forwards,
        strikes,
        expiries,
        implied_volatilities,
    ):

        if forwards.ndim != 1:
            raise FinError("Forwards must be one-dimensional.")

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if len(forwards) != len(expiries):
            raise FinError("One forward is required for each expiry.")

        if np.any(forwards <= 0.0):
            raise FinError("Forwards must be positive.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if np.any(expiries <= 0.0):
            raise FinError("Expiries must be positive.")

        if np.any(np.diff(expiries) <= 0.0):
            raise FinError("Expiries must be strictly increasing.")

        if not np.all(np.isfinite(forwards)):
            raise FinError("Forwards must be finite.")

        if not np.all(np.isfinite(strikes)):
            raise FinError("Strikes must be finite.")

        if not np.all(np.isfinite(expiries)):
            raise FinError("Expiries must be finite.")

        expected_shape = (
            len(expiries),
            len(strikes),
        )

        if implied_volatilities.shape != expected_shape:

            raise FinError("Implied-volatility surface has incorrect shape.")

        if not np.all(np.isfinite(implied_volatilities)):
            raise FinError("Implied volatilities must be finite.")

        if np.any(implied_volatilities <= 0.0):
            raise FinError("Implied volatilities must be positive.")

    ####################################################################################

    def _validate_parameters(
        self,
    ):

        if self._rho <= -1.0 or self._rho >= 1.0:

            raise FinError("Rho must lie between -1 and 1.")

        if self._expiries.ndim != 1:

            raise FinError("Expiries must be one-dimensional.")

        if self._atm_total_variances.ndim != 1:

            raise FinError("ATM total variances must be one-dimensional.")

        if len(self._expiries) != len(self._atm_total_variances):

            raise FinError("One ATM variance is required for each expiry.")

        if np.any(self._expiries <= 0.0):

            raise FinError("Expiries must be positive.")

        if np.any(np.diff(self._expiries) <= 0.0):

            raise FinError("Expiries must be strictly increasing.")

        if np.any(self._atm_total_variances <= 0.0):

            raise FinError("ATM total variances must be positive.")

        if np.any(np.diff(self._atm_total_variances) <= 0.0):

            raise FinError("ATM total variances must be strictly increasing.")

    ####################################################################################

    def _check_calibrated(
        self,
    ):

        if (
            self._expiries is None
            or self._atm_total_variances is None
            or self._rho is None
            or self._phi_function is None
        ):

            raise FinError("SSVI surface has not been calibrated.")

    ####################################################################################

    def __repr__(
        self,
    ):

        s = "OBJECT TYPE: SSVI SURFACE\n"

        if self._rho is None:

            s += "STATUS: NOT CALIBRATED\n"

            return s

        rho, eta, gamma = self.parameters()

        s += "RHO: " + str(rho) + "\n"

        s += "ETA: " + str(eta) + "\n"

        s += "GAMMA: " + str(gamma) + "\n"

        s += "EXPIRIES: " + str(self._expiries) + "\n"

        s += "ATM TOTAL VARIANCES: " + str(self._atm_total_variances) + "\n"

        return s
