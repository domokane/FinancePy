import numpy as np

from ..utils.error import FinError
from .svi import SVI
from .implied_volatility_surface import ImpliedVolatilitySurface

########################################################################################


class SVISurface(ImpliedVolatilitySurface):
    """
    SVI implied-volatility surface.

    The surface consists of one raw-SVI smile for each market expiry.

    Each smile is parameterized in forward log-moneyness

        k = log(K / F_T)

    using total implied variance

        w(k,T) = sigma_imp(k,T)^2 T.

    Between calibrated expiries, total variance is interpolated linearly
    in time.
    """

    def __init__(
        self,
        expiries=None,
        svi_parameters=None,
    ) -> None:

        self._expiries = None
        self._smiles = None

        if expiries is not None or svi_parameters is not None:

            if expiries is None or svi_parameters is None:
                raise FinError("Expiries and SVI parameters must both be supplied.")

            self._set_parameters(
                expiries,
                svi_parameters,
            )

    ####################################################################################

    def calibrate(
        self,
        forwards,
        strikes,
        expiries,
        implied_volatilities,
    ):
        """
        Calibrate one raw-SVI smile to each expiry.

        Parameters
        ----------
        forwards : array_like
            Forward price for each expiry.

        strikes : array_like
            Strike grid.

        expiries : array_like
            Expiry times.

        implied_volatilities : array_like
            Market implied-volatility surface with shape

                (num_expiries, num_strikes).

        Returns
        -------
        calibration_errors : ndarray
            Calibration cost for each maturity slice.
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

        smiles = []
        calibration_errors = np.empty(
            len(expiries),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            smile = SVI()

            calibration_errors[i] = smile.calibrate(
                forwards[i],
                strikes,
                implied_volatilities[i],
                t_exp,
            )

            smiles.append(smile)

        self._expiries = expiries
        self._smiles = smiles

        return calibration_errors

    ####################################################################################

    def total_variance(
        self,
        forward,
        strike,
        t_exp,
    ):
        """
        Return total implied variance

            w(K,T) = sigma_imp(K,T)^2 T.

        At calibrated expiries the corresponding SVI smile is evaluated
        directly. Between expiries, total variance at fixed forward
        log-moneyness is interpolated linearly in time.
        """

        self._check_calibrated()

        if forward <= 0.0:
            raise FinError("Forward must be positive.")

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        if t_exp < self._expiries[0]:
            raise FinError("Expiry is below SVI surface.")

        if t_exp > self._expiries[-1]:
            raise FinError("Expiry is above SVI surface.")

        k = np.log(strike / forward)

        ################################################################################
        # EXACT MARKET EXPIRY
        ################################################################################

        index = np.searchsorted(
            self._expiries,
            t_exp,
        )

        if index < len(self._expiries) and np.isclose(
            self._expiries[index],
            t_exp,
        ):

            return self._smiles[index].total_variance_from_log_moneyness(k)

        ################################################################################
        # INTERPOLATE BETWEEN MARKET EXPIRIES
        ################################################################################

        upper = index
        lower = index - 1

        t1 = self._expiries[lower]

        t2 = self._expiries[upper]

        w1 = self._smiles[lower].total_variance_from_log_moneyness(k)

        w2 = self._smiles[upper].total_variance_from_log_moneyness(k)

        alpha = (t_exp - t1) / (t2 - t1)

        return (1.0 - alpha) * w1 + alpha * w2

    ####################################################################################

    def implied_volatility(
        self,
        forward,
        strike,
        t_exp,
    ):
        """
        Return Black implied volatility.
        """

        w = self.total_variance(
            forward,
            strike,
            t_exp,
        )

        if not np.isfinite(w):
            return np.nan

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
        """
        Return an implied-volatility smile for a fixed expiry.
        """

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

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
        """
        Return an implied-volatility surface.

        The returned array has shape

            (num_expiries, num_strikes).
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

        if forwards.ndim != 1:
            raise FinError("Forwards must be one-dimensional.")

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

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

            for j, strike in enumerate(strikes):

                vols[i, j] = self.implied_volatility(
                    forwards[i],
                    strike,
                    t_exp,
                )

        return vols

    ####################################################################################

    def parameters(
        self,
    ):
        """
        Return calibrated raw-SVI parameters.

        Each row contains

            (a, b, rho, m, sigma).
        """

        self._check_calibrated()

        params = np.empty(
            (
                len(self._smiles),
                5,
            ),
            dtype=float,
        )

        for i, smile in enumerate(self._smiles):

            params[i, 0] = smile._a
            params[i, 1] = smile._b
            params[i, 2] = smile._rho
            params[i, 3] = smile._m
            params[i, 4] = smile._sigma

        return params

    ####################################################################################

    def _set_parameters(
        self,
        expiries,
        svi_parameters,
    ):

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        svi_parameters = np.asarray(
            svi_parameters,
            dtype=float,
        )

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if np.any(expiries <= 0.0):
            raise FinError("Expiries must be positive.")

        if np.any(np.diff(expiries) <= 0.0):
            raise FinError("Expiries must be strictly increasing.")

        expected_shape = (
            len(expiries),
            5,
        )

        if svi_parameters.shape != expected_shape:
            raise FinError("SVI parameters must have shape " "(num_expiries, 5).")

        smiles = []

        for params in svi_parameters:

            smiles.append(
                SVI(
                    a=params[0],
                    b=params[1],
                    rho=params[2],
                    m=params[3],
                    sigma=params[4],
                )
            )

        self._expiries = expiries
        self._smiles = smiles

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

    def _check_calibrated(
        self,
    ):

        if self._expiries is None or self._smiles is None:
            raise FinError("SVI surface has not been calibrated.")

    ####################################################################################

    def __repr__(
        self,
    ):

        s = "OBJECT TYPE: SVI SURFACE\n"

        if self._expiries is None:

            s += "STATUS: NOT CALIBRATED\n"

        else:

            s += "NUM EXPIRIES: " + str(len(self._expiries)) + "\n"

        return s
