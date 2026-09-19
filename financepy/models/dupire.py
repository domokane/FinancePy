# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from scipy.interpolate import RectBivariateSpline

from ..models.black_scholes_analytic import european_value
from ..utils.error import FinError
from ..utils.global_types import OptionTypes

########################################################################################


class Dupire:
    """
    Dupire local-volatility model.

    The model recovers local volatility from a smooth surface of European
    call-option prices.

    For constant continuously compounded interest rate r and dividend
    yield q,

        sigma_loc^2(K,T)

                 dC/dT + (r-q) K dC/dK + q C
        = ------------------------------------------------
                        0.5 K^2 d2C/dK2

    The input surface may be supplied either as:

        call_prices

    or as:

        implied_volatilities

    If implied volatilities are supplied, stock_price must also be supplied.
    The implied-volatility surface is converted internally into European
    call prices.

    The surface arrays must have shape

        (num_expiries, num_strikes).
    """

    def __init__(
        self,
        expiries,
        strikes,
        interest_rate,
        dividend_yield,
        call_prices=None,
        implied_volatilities=None,
        stock_price=None,
    ) -> None:

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        self._validate_grid(
            expiries,
            strikes,
        )

        if call_prices is None and implied_volatilities is None:
            raise FinError("Call prices or implied volatilities must be supplied.")

        if call_prices is not None and implied_volatilities is not None:
            raise FinError("Supply call prices or implied volatilities, not both.")

        if call_prices is not None:

            call_prices = np.asarray(
                call_prices,
                dtype=float,
            )

        else:

            if stock_price is None:
                raise FinError("Stock price is required when using implied volatilities.")

            if stock_price <= 0.0:
                raise FinError("Stock price must be positive.")

            implied_volatilities = np.asarray(
                implied_volatilities,
                dtype=float,
            )

            call_prices = self._vols_to_call_prices(
                stock_price,
                expiries,
                strikes,
                implied_volatilities,
                interest_rate,
                dividend_yield,
            )

        self._validate_call_prices(
            expiries,
            strikes,
            call_prices,
        )

        self._expiries = expiries
        self._strikes = strikes
        self._call_prices = call_prices

        self._interest_rate = interest_rate
        self._dividend_yield = dividend_yield

        self._call_spline = RectBivariateSpline(
            expiries,
            strikes,
            call_prices,
            kx=3,
            ky=3,
            s=0.0,
        )

    ####################################################################################

    def call_value(
        self,
        strike,
        t_exp,
    ):
        """
        Return the interpolated European call-option price.
        """

        self._validate_point(
            strike,
            t_exp,
        )

        return float(
            self._call_spline(
                t_exp,
                strike,
                dx=0,
                dy=0,
            )[0, 0]
        )

    ####################################################################################

    def local_variance(
        self,
        strike,
        t_exp,
    ):
        """
        Return Dupire local variance at state level S = strike and
        time t = t_exp.
        """

        self._validate_point(
            strike,
            t_exp,
        )

        (
            call_value,
            dcdt,
            dcdk,
            d2cdk2,
        ) = self._call_derivatives(
            strike,
            t_exp,
        )

        r = self._interest_rate
        q = self._dividend_yield

        numerator = dcdt + (r - q) * strike * dcdk + q * call_value

        denominator = 0.5 * strike * strike * d2cdk2

        if not np.isfinite(numerator):
            return np.nan

        if not np.isfinite(denominator):
            return np.nan

        # Positive call-price convexity in strike is required.
        if denominator <= 0.0:
            return np.nan

        local_variance = numerator / denominator

        if not np.isfinite(local_variance):
            return np.nan

        if local_variance <= 0.0:
            return np.nan

        return local_variance

    ####################################################################################

    def local_volatility(
        self,
        strike,
        t_exp,
    ):
        """
        Return Dupire local volatility at state level S = strike and
        time t = t_exp.
        """

        local_variance = self.local_variance(
            strike,
            t_exp,
        )

        if not np.isfinite(local_variance):
            return np.nan

        return np.sqrt(local_variance)

    ####################################################################################

    def local_volatility_curve(
        self,
        strikes,
        t_exp,
    ):
        """
        Return local volatility across strike for a fixed expiry.
        """

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        local_vols = np.empty(
            len(strikes),
            dtype=float,
        )

        for i, strike in enumerate(strikes):

            local_vols[i] = self.local_volatility(
                strike,
                t_exp,
            )

        return local_vols

    ####################################################################################

    def local_volatility_surface(
        self,
        strikes,
        expiries,
    ):
        """
        Return the local-volatility surface.

        The returned array has shape

            (num_expiries, num_strikes).
        """

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if np.any(expiries <= 0.0):
            raise FinError("Expiries must be positive.")

        local_vols = np.empty(
            (
                len(expiries),
                len(strikes),
            ),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            for j, strike in enumerate(strikes):

                local_vols[i, j] = self.local_volatility(
                    strike,
                    t_exp,
                )

        return local_vols

    ####################################################################################

    def _call_derivatives(
        self,
        strike,
        t_exp,
    ):
        """
        Return

            C,
            dC/dT,
            dC/dK,
            d2C/dK2.
        """

        call_value = float(
            self._call_spline(
                t_exp,
                strike,
                dx=0,
                dy=0,
            )[0, 0]
        )

        dcdt = float(
            self._call_spline(
                t_exp,
                strike,
                dx=1,
                dy=0,
            )[0, 0]
        )

        dcdk = float(
            self._call_spline(
                t_exp,
                strike,
                dx=0,
                dy=1,
            )[0, 0]
        )

        d2cdk2 = float(
            self._call_spline(
                t_exp,
                strike,
                dx=0,
                dy=2,
            )[0, 0]
        )

        return (
            call_value,
            dcdt,
            dcdk,
            d2cdk2,
        )

    ####################################################################################

    def _vols_to_call_prices(
        self,
        stock_price,
        expiries,
        strikes,
        implied_volatilities,
        interest_rate,
        dividend_yield,
    ):
        """
        Convert a Black-Scholes implied-volatility surface into a
        European call-price surface.
        """

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

        call_prices = np.empty(
            expected_shape,
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            for j, strike in enumerate(strikes):

                call_prices[i, j] = european_value(
                    stock_price,
                    t_exp,
                    strike,
                    interest_rate,
                    dividend_yield,
                    implied_volatilities[i, j],
                    OptionTypes.EUROPEAN_CALL.value,
                )

        return call_prices

    ####################################################################################

    def _validate_grid(
        self,
        expiries,
        strikes,
    ):

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if len(expiries) < 4:
            raise FinError("At least four expiries are required.")

        if len(strikes) < 4:
            raise FinError("At least four strikes are required.")

        if np.any(expiries <= 0.0):
            raise FinError("Expiries must be positive.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if np.any(np.diff(expiries) <= 0.0):
            raise FinError("Expiries must be strictly increasing.")

        if np.any(np.diff(strikes) <= 0.0):
            raise FinError("Strikes must be strictly increasing.")

    ####################################################################################

    def _validate_call_prices(
        self,
        expiries,
        strikes,
        call_prices,
    ):

        expected_shape = (
            len(expiries),
            len(strikes),
        )

        if call_prices.shape != expected_shape:
            raise FinError("Call-price surface has incorrect shape.")

        if not np.all(np.isfinite(call_prices)):
            raise FinError("Call prices must be finite.")

        if np.any(call_prices < 0.0):
            raise FinError("Call prices cannot be negative.")

    ####################################################################################

    def _validate_point(
        self,
        strike,
        t_exp,
    ):

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        if strike < self._strikes[0]:
            raise FinError("Strike is below Dupire surface.")

        if strike > self._strikes[-1]:
            raise FinError("Strike is above Dupire surface.")

        if t_exp < self._expiries[0]:
            raise FinError("Expiry is below Dupire surface.")

        if t_exp > self._expiries[-1]:
            raise FinError("Expiry is above Dupire surface.")

    ####################################################################################

    def __repr__(self):

        s = "OBJECT TYPE: Dupire\n"

        s += "NUM EXPIRIES: " + str(len(self._expiries)) + "\n"

        s += "NUM STRIKES: " + str(len(self._strikes)) + "\n"

        s += "INTEREST RATE: " + str(self._interest_rate) + "\n"

        s += "DIVIDEND YIELD: " + str(self._dividend_yield) + "\n"

        return s
