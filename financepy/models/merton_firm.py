# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union
import numpy as np

from ..utils.error import FinError
from ..utils.helpers import label_to_string, check_argument_types
from ..utils.math import normcdf_vect

# TODO: Redesign this class


########################################################################################
class MertonFirm:
    """
    Merton structural firm-value model.

    Inputs may be scalars or NumPy arrays and are broadcast using NumPy rules.

    Parameters
    ----------
    asset_value : float or np.ndarray
        Market value of the firm's assets, A(t).

    debt_face_value : float or np.ndarray
        Face value of zero-coupon debt, F.

    time_to_maturity : float or np.ndarray
        Time to debt maturity, tau = T - t, in years.

    risk_free_rate : float or np.ndarray
        Continuously compounded risk-free rate.

    asset_drift : float or np.ndarray
        Physical expected return on the firm's assets.

    asset_volatility : float or np.ndarray
        Volatility of the firm's asset value.
    """

    def __init__(
        self,
        asset_value: Union[float, np.ndarray],
        debt_face_value: Union[float, np.ndarray],
        time_to_maturity: Union[float, np.ndarray],
        risk_free_rate: Union[float, np.ndarray],
        asset_drift: Union[float, np.ndarray],
        asset_volatility: Union[float, np.ndarray],
    ) -> None:

        check_argument_types(MertonFirm.__init__, locals())

        self._asset_value = np.asarray(asset_value, dtype=float)
        self._debt_face_value = np.asarray(debt_face_value, dtype=float)
        self._time_to_maturity = np.asarray(time_to_maturity, dtype=float)
        self._risk_free_rate = np.asarray(risk_free_rate, dtype=float)
        self._asset_drift = np.asarray(asset_drift, dtype=float)
        self._asset_volatility = np.asarray(asset_volatility, dtype=float)

        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """Validate model inputs."""

        if np.any(self._asset_value <= 0.0):
            raise FinError("Asset value must be positive.")

        if np.any(self._debt_face_value <= 0.0):
            raise FinError("Debt face value must be positive.")

        if np.any(self._time_to_maturity <= 0.0):
            raise FinError("Time to maturity must be positive.")

        if np.any(self._asset_volatility <= 0.0):
            raise FinError("Asset volatility must be positive.")

        try:
            np.broadcast_arrays(
                self._asset_value,
                self._debt_face_value,
                self._time_to_maturity,
                self._risk_free_rate,
                self._asset_drift,
                self._asset_volatility,
            )
        except ValueError as exc:
            raise FinError("Model inputs are not broadcast-compatible.") from exc

    def _d1_d2(self):
        """Return the Black-Scholes d1 and d2 terms."""

        a = self._asset_value
        f = self._debt_face_value
        tau = self._time_to_maturity
        r = self._risk_free_rate
        sigma_a = self._asset_volatility

        sigma_root_tau = sigma_a * np.sqrt(tau)

        d1 = (np.log(a / f) + (r + 0.5 * sigma_a**2) * tau) / sigma_root_tau

        d2 = d1 - sigma_root_tau

        return d1, d2

    def asset_value(self) -> np.ndarray:
        """Return the firm's asset value."""

        return self._asset_value

    def debt_face_value(self) -> np.ndarray:
        """Return the debt face value."""

        return self._debt_face_value

    def time_to_maturity(self) -> np.ndarray:
        """Return time to maturity."""

        return self._time_to_maturity

    def risk_free_rate(self) -> np.ndarray:
        """Return the risk-free rate."""

        return self._risk_free_rate

    def asset_drift(self) -> np.ndarray:
        """Return the physical asset drift."""

        return self._asset_drift

    def asset_volatility(self) -> np.ndarray:
        """Return the asset volatility."""

        return self._asset_volatility

    def asset_to_debt_ratio(self) -> np.ndarray:
        """Return A(t) / F."""

        return self._asset_value / self._debt_face_value

    def equity_value(self) -> np.ndarray:
        """Return the market value of equity."""

        d1, d2 = self._d1_d2()

        a = self._asset_value
        f = self._debt_face_value
        tau = self._time_to_maturity
        r = self._risk_free_rate

        return a * normcdf_vect(d1) - f * np.exp(-r * tau) * normcdf_vect(d2)

    def debt_value(self) -> np.ndarray:
        """Return the market value of risky zero-coupon debt."""

        d1, d2 = self._d1_d2()

        a = self._asset_value
        f = self._debt_face_value
        tau = self._time_to_maturity
        r = self._risk_free_rate

        return a * normcdf_vect(-d1) + f * np.exp(-r * tau) * normcdf_vect(d2)

    def equity_volatility(self) -> np.ndarray:
        """Return the equity volatility implied by the Merton model."""

        d1, _ = self._d1_d2()

        a = self._asset_value
        sigma_a = self._asset_volatility
        e = self.equity_value()

        return (a / e) * normcdf_vect(d1) * sigma_a

    def risky_yield(self) -> np.ndarray:
        """Return the continuously compounded risky debt yield."""

        d = self.debt_value()
        f = self._debt_face_value
        tau = self._time_to_maturity

        return -(1.0 / tau) * np.log(d / f)

    def credit_spread(self) -> np.ndarray:
        """Return the continuously compounded credit spread."""

        return self.risky_yield() - self._risk_free_rate

    def distance_to_default(self) -> np.ndarray:
        """
        Return the physical-measure distance to default.

        This uses the physical asset drift rather than the risk-free rate.
        """

        a = self._asset_value
        f = self._debt_face_value
        tau = self._time_to_maturity
        mu = self._asset_drift
        sigma_a = self._asset_volatility

        return (np.log(a / f) + (mu - 0.5 * sigma_a**2) * tau) / (sigma_a * np.sqrt(tau))

    def physical_default_probability(self) -> np.ndarray:
        """Return the physical probability of default by maturity."""

        dd = self.distance_to_default()

        return normcdf_vect(-dd)

    def risk_neutral_default_probability(self) -> np.ndarray:
        """Return the risk-neutral probability of default by maturity."""

        _, d2 = self._d1_d2()

        return normcdf_vect(-d2)

    def __repr__(self) -> str:

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("ASSET VALUE", self._asset_value)
        s += label_to_string("DEBT FACE VALUE", self._debt_face_value)
        s += label_to_string("TIME TO MATURITY", self._time_to_maturity)
        s += label_to_string("RISK-FREE RATE", self._risk_free_rate)
        s += label_to_string("ASSET DRIFT", self._asset_drift)
        s += label_to_string("ASSET VOLATILITY", self._asset_volatility)

        return s
