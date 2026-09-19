###############################################################################
# Copyright (C) 2026
#
# Merton Jump-Diffusion Model
#
# Risk-neutral dynamics:
#
#     dS_t / S_{t-}
#         = (r - q - lambda * kappa_j) dt
#           + sigma dW_t
#           + (J - 1) dN_t
#
#     log(J) ~ N(mu_j, delta_j^2)
#
#     kappa_j = E[J - 1]
#             = exp(mu_j + 0.5 * delta_j^2) - 1
#
###############################################################################

import numpy as np

from scipy.optimize import least_squares

from financepy.models.black_scholes_analytic import value
from financepy.models.black_scholes_analytic import vega
from financepy.models.black_scholes_analytic import implied_volatility
from financepy.utils.global_types import OptionTypes

###############################################################################
# Helpers
###############################################################################


def _as_option_type_value(option_type):
    """Convert OptionTypes enum or integer to FinancePy integer value."""

    if isinstance(option_type, OptionTypes):
        return option_type.value

    return int(option_type)


###############################################################################


def _expand_term_structure(x, num_expiries):
    """Allow r and q to be supplied as scalars or arrays by expiry."""

    if np.isscalar(x):
        return np.full(num_expiries, float(x))

    x = np.asarray(x, dtype=float)

    if len(x) != num_expiries:
        raise ValueError("Term structure must contain one value per expiry.")

    return x


###############################################################################


def _check_surface_inputs(
    strikes,
    expiries,
    market_vols,
):
    """Validate volatility-surface dimensions."""

    strikes = np.asarray(strikes, dtype=float)
    expiries = np.asarray(expiries, dtype=float)
    market_vols = np.asarray(market_vols, dtype=float)

    if strikes.ndim != 1:
        raise ValueError("strikes must be a one-dimensional array.")

    if expiries.ndim != 1:
        raise ValueError("expiries must be a one-dimensional array.")

    if market_vols.shape != (len(expiries), len(strikes)):
        raise ValueError("market_vols must have shape " "(len(expiries), len(strikes)).")

    if np.any(strikes <= 0.0):
        raise ValueError("All strikes must be positive.")

    if np.any(expiries <= 0.0):
        raise ValueError("All expiries must be positive.")

    if np.any(market_vols <= 0.0):
        raise ValueError("All market volatilities must be positive.")

    return strikes, expiries, market_vols


###############################################################################
# Merton Model
###############################################################################


class MertonJumpDiffusion:
    """
    Merton jump-diffusion model for European equity options.

    The risk-neutral stock-price process is

        dS/S_- =
            (r - q - lambda * kappa_j) dt
            + sigma dW
            + (J - 1) dN

    where

        log(J) ~ Normal(mu_j, delta_j^2)

    and

        kappa_j = E[J - 1]
                = exp(mu_j + 0.5 * delta_j^2) - 1.

    Parameters
    ----------
    sigma : float
        Diffusive volatility.

    jump_intensity : float
        Poisson jump intensity lambda, in jumps per year.

    jump_mean : float
        Mean log jump size mu_j.

    jump_volatility : float
        Standard deviation delta_j of log jump sizes.

    poisson_tolerance : float
        Remaining Poisson probability below which the pricing
        summation is terminated.

    max_jumps : int
        Maximum number of jump terms in the pricing summation.
    """

    ###########################################################################

    def __init__(
        self,
        sigma,
        jump_intensity,
        jump_mean,
        jump_volatility,
        poisson_tolerance=1.0e-12,
        max_jumps=200,
    ):

        self.sigma = float(sigma)
        self.jump_intensity = float(jump_intensity)
        self.jump_mean = float(jump_mean)
        self.jump_volatility = float(jump_volatility)

        self.poisson_tolerance = float(poisson_tolerance)
        self.max_jumps = int(max_jumps)

        self._validate_parameters()

    ###########################################################################

    def _validate_parameters(self):

        if self.sigma < 0.0:
            raise ValueError("sigma must be non-negative.")

        if self.jump_intensity < 0.0:
            raise ValueError("jump_intensity must be non-negative.")

        if self.jump_volatility < 0.0:
            raise ValueError("jump_volatility must be non-negative.")

        if self.poisson_tolerance <= 0.0:
            raise ValueError("poisson_tolerance must be positive.")

        if self.max_jumps < 1:
            raise ValueError("max_jumps must be at least one.")

    ###########################################################################

    @property
    def jump_compensator(self):
        """
        Expected proportional jump:

            E[J - 1]
            =
            exp(mu_j + 0.5 delta_j^2) - 1.
        """

        mu_j = self.jump_mean
        delta_j = self.jump_volatility

        return np.exp(mu_j + 0.5 * delta_j * delta_j) - 1.0

    ###########################################################################

    @property
    def expected_jump_multiplier(self):
        """Return E[J]."""

        return self.jump_compensator + 1.0

    ###########################################################################

    def value(
        self,
        stock_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_yield,
        option_type=OptionTypes.EUROPEAN_CALL,
    ):
        """
        Value a European option using Merton's Poisson-mixture solution.

        Conditional on N_T=n jumps, log(S_T) is Gaussian. Therefore
        the option value is a Poisson-weighted sum of Black-Scholes
        values.
        """

        s = float(stock_price)
        t = float(time_to_expiry)
        k = float(strike_price)
        r = float(risk_free_rate)
        q = float(dividend_yield)

        opt_type_value = _as_option_type_value(option_type)

        if t <= 0.0:

            if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
                return max(s - k, 0.0)

            if opt_type_value == OptionTypes.EUROPEAN_PUT.value:
                return max(k - s, 0.0)

            raise ValueError("Only European calls and puts are supported.")

        sigma = self.sigma
        lam = self.jump_intensity
        mu_j = self.jump_mean
        delta_j = self.jump_volatility

        kappa_j = self.jump_compensator

        lambda_t = lam * t

        # P(N_T = 0)
        poisson_prob = np.exp(-lambda_t)

        cumulative_prob = 0.0
        option_value = 0.0

        for n in range(self.max_jumps + 1):

            # --------------------------------------------------------------
            # Conditional total variance
            #
            # sigma_n^2 T
            # =
            # sigma^2 T + n delta_j^2
            # --------------------------------------------------------------

            variance_n = sigma * sigma + n * delta_j * delta_j / t

            sigma_n = np.sqrt(max(variance_n, 1.0e-16))

            # --------------------------------------------------------------
            # Conditional forward
            #
            # E[S_T | N_T=n]
            #
            # = S exp[
            #     (r-q-lambda*kappa_j)T
            #     + n(mu_j + delta_j^2/2)
            #   ]
            #
            # bs_value uses
            #
            #     F = S exp[(r-q_n)T]
            #
            # so define q_n such that the forward is correct while
            # retaining the actual discount rate r.
            # --------------------------------------------------------------

            q_n = q + lam * kappa_j - n * (mu_j + 0.5 * delta_j * delta_j) / t

            value_n = value(
                s,
                t,
                k,
                r,
                q_n,
                sigma_n,
                opt_type_value,
            )

            option_value += poisson_prob * value_n

            cumulative_prob += poisson_prob

            # Stop once the remaining Poisson mass is negligible.
            if n >= lambda_t and 1.0 - cumulative_prob < self.poisson_tolerance:
                break

            # Recursive Poisson probability:
            #
            # p_(n+1) = p_n lambda T / (n+1)
            #
            poisson_prob *= lambda_t / (n + 1.0)

        return float(option_value)

    ###########################################################################

    def implied_volatility(
        self,
        stock_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_yield,
        option_type=OptionTypes.EUROPEAN_CALL,
    ):
        """Return Black-Scholes implied volatility of the Merton price."""

        opt_type_value = _as_option_type_value(option_type)

        price = self.value(
            stock_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_yield,
            option_type,
        )

        vol = implied_volatility(
            stock_price,
            time_to_expiry,
            strike_price,
            risk_free_rate,
            dividend_yield,
            price,
            opt_type_value,
        )

        return float(vol)

    ###########################################################################

    def volatility_smile(
        self,
        stock_price,
        time_to_expiry,
        strikes,
        risk_free_rate,
        dividend_yield,
        option_type=OptionTypes.EUROPEAN_CALL,
    ):
        """Generate the Merton Black-Scholes implied-volatility smile."""

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        vols = np.empty(len(strikes))

        for i, strike in enumerate(strikes):

            vols[i] = self.implied_volatility(
                stock_price,
                time_to_expiry,
                strike,
                risk_free_rate,
                dividend_yield,
                option_type,
            )

        return vols

    ###########################################################################

    def volatility_surface(
        self,
        stock_price,
        strikes,
        expiries,
        risk_free_rates,
        dividend_yields,
        option_type=OptionTypes.EUROPEAN_CALL,
    ):
        """
        Generate a matrix of Black-Scholes implied volatilities.

        Returns
        -------
        vols : ndarray
            Shape = (number of expiries, number of strikes).
        """

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        num_expiries = len(expiries)
        num_strikes = len(strikes)

        rates = _expand_term_structure(
            risk_free_rates,
            num_expiries,
        )

        dividends = _expand_term_structure(
            dividend_yields,
            num_expiries,
        )

        vols = np.empty((num_expiries, num_strikes))

        for i, t in enumerate(expiries):

            for j, k in enumerate(strikes):

                vols[i, j] = self.implied_volatility(
                    stock_price,
                    t,
                    k,
                    rates[i],
                    dividends[i],
                    option_type,
                )

        return vols

    ###########################################################################

    @classmethod
    def calibrate(
        cls,
        stock_price,
        strikes,
        expiries,
        market_vols,
        risk_free_rates,
        dividend_yields,
        option_type=OptionTypes.EUROPEAN_CALL,
        initial_guess=None,
        lower_bounds=None,
        upper_bounds=None,
        weights=None,
        calibration_type="VOL",
        max_nfev=3000,
    ):
        """
        Calibrate one global Merton parameter set to an entire
        implied-volatility surface.

        The calibrated parameters are

            sigma
            jump_intensity lambda
            jump_mean mu_j
            jump_volatility delta_j

        Parameters
        ----------
        market_vols : ndarray
            Matrix of Black-Scholes implied volatilities with shape

                (len(expiries), len(strikes))

        calibration_type : str
            "VOL"
                Minimise errors directly in Black-Scholes implied
                volatility.

            "VEGA"
                Minimise Black-Scholes vega-scaled price errors.

        weights : ndarray or None
            Optional calibration weights with the same dimensions
            as market_vols.

        Returns
        -------
        MertonCalibrationResult
        """

        strikes, expiries, market_vols = _check_surface_inputs(
            strikes,
            expiries,
            market_vols,
        )

        s = float(stock_price)

        num_expiries = len(expiries)
        num_strikes = len(strikes)

        rates = _expand_term_structure(
            risk_free_rates,
            num_expiries,
        )

        dividends = _expand_term_structure(
            dividend_yields,
            num_expiries,
        )

        opt_type_value = _as_option_type_value(option_type)

        # --------------------------------------------------------------
        # Calibration weights
        # --------------------------------------------------------------

        if weights is None:

            weights = np.ones_like(market_vols)

        else:

            weights = np.asarray(
                weights,
                dtype=float,
            )

            if weights.shape != market_vols.shape:
                raise ValueError("weights must have the same shape " "as market_vols.")

        # --------------------------------------------------------------
        # Initial parameters
        # --------------------------------------------------------------

        if initial_guess is None:

            # Use the volatility nearest spot as a rough starting
            # value for the diffusion volatility.

            atm_idx = int(np.argmin(np.abs(strikes - s)))

            initial_sigma = float(np.min(market_vols[:, atm_idx]))

            initial_guess = np.array(
                [
                    initial_sigma,
                    0.50,
                    -0.10,
                    0.20,
                ],
                dtype=float,
            )

        else:

            initial_guess = np.asarray(
                initial_guess,
                dtype=float,
            )

        # --------------------------------------------------------------
        # Parameter bounds
        # --------------------------------------------------------------

        if lower_bounds is None:

            lower_bounds = np.array(
                [
                    0.001,  # sigma
                    0.000,  # lambda
                    -2.000,  # mu_j
                    0.001,  # delta_j
                ],
                dtype=float,
            )

        if upper_bounds is None:

            upper_bounds = np.array(
                [
                    2.00,  # sigma
                    20.00,  # lambda
                    1.00,  # mu_j
                    2.00,  # delta_j
                ],
                dtype=float,
            )

        lower_bounds = np.asarray(
            lower_bounds,
            dtype=float,
        )

        upper_bounds = np.asarray(
            upper_bounds,
            dtype=float,
        )

        # --------------------------------------------------------------
        # Generate market prices and market vegas once.
        #
        # This is required for the VEGA calibration objective and
        # is useful for diagnostics even when calibrating in VOL.
        # --------------------------------------------------------------

        market_prices = np.empty_like(market_vols)

        market_vegas = np.empty_like(market_vols)

        for i in range(num_expiries):

            t = expiries[i]
            r = rates[i]
            q = dividends[i]

            for j in range(num_strikes):

                k = strikes[j]
                vol = market_vols[i, j]

                market_prices[i, j] = value(
                    s,
                    t,
                    k,
                    r,
                    q,
                    vol,
                    opt_type_value,
                )

                market_vegas[i, j] = vega(
                    s,
                    t,
                    k,
                    r,
                    q,
                    vol,
                    opt_type_value,
                )

        # Avoid division by essentially zero vega.
        vega_floor = max(
            1.0e-8,
            1.0e-6 * s,
        )

        market_vegas = np.maximum(
            market_vegas,
            vega_floor,
        )

        calibration_type = calibration_type.upper()

        # --------------------------------------------------------------
        # Objective
        # --------------------------------------------------------------

        def residuals(params):

            sigma = params[0]
            lam = params[1]
            mu_j = params[2]
            delta_j = params[3]

            model = cls(
                sigma=sigma,
                jump_intensity=lam,
                jump_mean=mu_j,
                jump_volatility=delta_j,
            )

            errors = np.empty(num_expiries * num_strikes)

            n = 0

            for i in range(num_expiries):

                t = expiries[i]
                r = rates[i]
                q = dividends[i]

                for j in range(num_strikes):

                    k = strikes[j]

                    model_price = model.value(
                        s,
                        t,
                        k,
                        r,
                        q,
                        option_type,
                    )

                    if calibration_type == "VOL":

                        model_vol = implied_volatility(
                            s,
                            t,
                            k,
                            r,
                            q,
                            model_price,
                            opt_type_value,
                        )

                        if model_vol is None or not np.isfinite(model_vol):
                            error = 1.0

                        else:

                            error = model_vol - market_vols[i, j]

                    elif calibration_type == "VEGA":

                        error = (model_price - market_prices[i, j]) / market_vegas[i, j]

                    else:

                        raise ValueError("calibration_type must be " "'VOL' or 'VEGA'.")

                    errors[n] = weights[i, j] * error

                    n += 1

            return errors

        # --------------------------------------------------------------
        # Optimisation
        # --------------------------------------------------------------

        result = least_squares(
            residuals,
            x0=initial_guess,
            bounds=(
                lower_bounds,
                upper_bounds,
            ),
            xtol=1.0e-10,
            ftol=1.0e-10,
            gtol=1.0e-10,
            max_nfev=max_nfev,
        )

        sigma = result.x[0]
        lam = result.x[1]
        mu_j = result.x[2]
        delta_j = result.x[3]

        model = cls(
            sigma=sigma,
            jump_intensity=lam,
            jump_mean=mu_j,
            jump_volatility=delta_j,
        )

        # --------------------------------------------------------------
        # Generate exact fitted implied-volatility surface.
        # --------------------------------------------------------------

        fitted_vols = model.volatility_surface(
            s,
            strikes,
            expiries,
            rates,
            dividends,
            option_type,
        )

        vol_errors = fitted_vols - market_vols

        rmse = np.sqrt(np.mean(vol_errors * vol_errors))

        mae = np.mean(np.abs(vol_errors))

        max_abs_error = np.max(np.abs(vol_errors))

        return MertonCalibrationResult(
            model=model,
            sigma=sigma,
            jump_intensity=lam,
            jump_mean=mu_j,
            jump_volatility=delta_j,
            fitted_vols=fitted_vols,
            market_vols=market_vols.copy(),
            vol_errors=vol_errors,
            rmse=rmse,
            mae=mae,
            max_abs_error=max_abs_error,
            optimizer_result=result,
        )

    ###########################################################################

    def __repr__(self):

        return (
            "MertonJumpDiffusion("
            f"sigma={self.sigma:.6f}, "
            f"jump_intensity={self.jump_intensity:.6f}, "
            f"jump_mean={self.jump_mean:.6f}, "
            f"jump_volatility={self.jump_volatility:.6f})"
        )


###############################################################################
# Calibration Result
###############################################################################


class MertonCalibrationResult:
    """Container for a Merton volatility-surface calibration."""

    ###########################################################################

    def __init__(
        self,
        model,
        sigma,
        jump_intensity,
        jump_mean,
        jump_volatility,
        fitted_vols,
        market_vols,
        vol_errors,
        rmse,
        mae,
        max_abs_error,
        optimizer_result,
    ):

        self.model = model

        self.sigma = sigma
        self.jump_intensity = jump_intensity
        self.jump_mean = jump_mean
        self.jump_volatility = jump_volatility

        self.fitted_vols = fitted_vols
        self.market_vols = market_vols
        self.vol_errors = vol_errors

        self.rmse = rmse
        self.mae = mae
        self.max_abs_error = max_abs_error

        self.optimizer_result = optimizer_result

    ###########################################################################

    @property
    def success(self):
        return self.optimizer_result.success

    ###########################################################################

    @property
    def jump_compensator(self):

        return np.exp(self.jump_mean + 0.5 * self.jump_volatility * self.jump_volatility) - 1.0

    ###########################################################################

    @property
    def expected_jump_size(self):
        """Expected proportional jump E[J-1]."""

        return self.jump_compensator

    ###########################################################################

    def __repr__(self):

        return (
            "MertonCalibrationResult(\n"
            f"  sigma              = {self.sigma:.8f}\n"
            f"  jump_intensity     = {self.jump_intensity:.8f}\n"
            f"  jump_mean          = {self.jump_mean:.8f}\n"
            f"  jump_volatility    = {self.jump_volatility:.8f}\n"
            f"  expected_jump_size = {self.expected_jump_size:.8f}\n"
            f"  vol_rmse           = {self.rmse:.8f}\n"
            f"  vol_mae            = {self.mae:.8f}\n"
            f"  max_abs_vol_error  = {self.max_abs_error:.8f}\n"
            f"  success            = {self.success}\n"
            ")"
        )


###############################################################################
# Example
###############################################################################

if __name__ == "__main__":

    stock_price = 100.0
    risk_free_rate = 0.05
    dividend_yield = 0.02

    strikes = np.array(
        [
            70.0,
            80.0,
            90.0,
            100.0,
            110.0,
            120.0,
            130.0,
        ]
    )

    expiries = np.array(
        [
            0.25,
            0.50,
            1.00,
            2.00,
            5.00,
        ]
    )

    market_vols = np.array(
        [
            [0.300, 0.270, 0.240, 0.210, 0.195, 0.187, 0.185],
            [0.285, 0.260, 0.232, 0.205, 0.192, 0.185, 0.184],
            [0.270, 0.250, 0.225, 0.200, 0.190, 0.185, 0.185],
            [0.250, 0.235, 0.215, 0.198, 0.191, 0.189, 0.190],
            [0.230, 0.220, 0.208, 0.197, 0.193, 0.191, 0.192],
        ]
    )

    calibration = MertonJumpDiffusion.calibrate(
        stock_price=stock_price,
        strikes=strikes,
        expiries=expiries,
        market_vols=market_vols,
        risk_free_rates=risk_free_rate,
        dividend_yields=dividend_yield,
        option_type=OptionTypes.EUROPEAN_CALL,
        calibration_type="VOL",
    )

    print(calibration)

    print("\nMarket vols:")
    print(calibration.market_vols)

    print("\nFitted vols:")
    print(calibration.fitted_vols)

    print("\nVol errors:")
    print(calibration.vol_errors)
