# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union, Any

import numpy as np
from scipy import optimize

from ..utils.error import FinError
from ..utils.helpers import check_argument_types, label_to_string
from ..utils.math import normcdf
from .merton_firm import MertonFirm


def _merton_equations(
    x: np.ndarray,
    equity_value: float,
    equity_volatility: float,
    bond_face: float,
    years_to_maturity: float,
    risk_free_rate: float,
) -> np.ndarray:
    """
    Equations used to infer asset value and asset volatility from
    observed equity value and equity volatility.
    """

    asset_value, asset_volatility = x

    sigma = asset_volatility
    a = asset_value
    e = equity_value
    r = risk_free_rate
    t = years_to_maturity
    f = bond_face

    if a <= 0.0 or sigma <= 0.0:
        return np.array([1.0e10, 1.0e10])

    sigma_root_t = sigma * np.sqrt(t)

    d1 = (np.log(a / f) + (r + 0.5 * sigma**2) * t) / sigma_root_t
    d2 = d1 - sigma_root_t

    model_equity_value = a * normcdf(d1) - f * np.exp(-r * t) * normcdf(d2)
    model_equity_volatility = a / e * normcdf(d1) * sigma

    return np.array(
        [
            model_equity_value - equity_value,
            model_equity_volatility - equity_volatility,
        ]
    )


class MertonFirmMkt(MertonFirm):
    """
    Market implementation of the Merton firm-value model.

    The observable inputs are equity value and equity volatility. The firm's
    asset value and asset volatility are inferred by solving the Merton equity
    value and equity volatility equations simultaneously.

    Parameters may be scalars or NumPy arrays. NumPy broadcasting rules are
    applied across inputs.
    """

    def __init__(
        self,
        equity_value: Union[float, np.ndarray],
        bond_face: Union[float, np.ndarray],
        years_to_maturity: Union[float, np.ndarray],
        risk_free_rate: Union[float, np.ndarray],
        asset_growth_rate: Union[float, np.ndarray],
        equity_volatility: Union[float, np.ndarray],
    ) -> None:

        check_argument_types(self.__init__, locals())

        equity_value = np.asarray(equity_value, dtype=float)
        bond_face = np.asarray(bond_face, dtype=float)
        years_to_maturity = np.asarray(years_to_maturity, dtype=float)
        risk_free_rate = np.asarray(risk_free_rate, dtype=float)
        asset_growth_rate = np.asarray(asset_growth_rate, dtype=float)
        equity_volatility = np.asarray(equity_volatility, dtype=float)

        self._validate_market_inputs(
            equity_value,
            bond_face,
            years_to_maturity,
            equity_volatility,
        )

        try:
            (
                self._e,
                self._l,
                self._t,
                self._r,
                self._mu,
                self._ve,
            ) = np.broadcast_arrays(
                equity_value,
                bond_face,
                years_to_maturity,
                risk_free_rate,
                asset_growth_rate,
                equity_volatility,
            )
        except ValueError as exc:
            raise FinError("MertonFirmMkt inputs are not broadcast-compatible.") from exc

        asset_value, asset_volatility = self._solve_for_asset_value_and_volatility()

        # Preserve observed market quantities before initialising parent.
        market_equity_value = self._e.copy()
        market_equity_volatility = self._ve.copy()

        # Initialise the parent class with the inferred asset quantities.
        super().__init__(
            asset_value=asset_value,
            debt_face_value=self._l,
            time_to_maturity=self._t,
            risk_free_rate=self._r,
            asset_drift=self._mu,
            asset_volatility=asset_volatility,
        )

        # Store the observed market quantities separately.
        self._market_equity_value = market_equity_value
        self._market_equity_volatility = market_equity_volatility

    @staticmethod
    def _validate_market_inputs(
        equity_value: np.ndarray,
        bond_face: np.ndarray,
        years_to_maturity: np.ndarray,
        equity_volatility: np.ndarray,
    ) -> None:
        """Validate observable market inputs."""

        if np.any(equity_value <= 0.0):
            raise FinError("Equity value must be positive.")

        if np.any(bond_face <= 0.0):
            raise FinError("Bond face value must be positive.")

        if np.any(years_to_maturity <= 0.0):
            raise FinError("Years to maturity must be positive.")

        if np.any(equity_volatility <= 0.0):
            raise FinError("Equity volatility must be positive.")

    def _solve_for_asset_value_and_volatility(
        self,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Infer asset value and asset volatility point by point."""

        shape = self._e.shape

        asset_values = np.empty(shape, dtype=float)
        asset_volatilities = np.empty(shape, dtype=float)

        iterator = np.ndindex(shape)

        for idx in iterator:

            e = float(self._e[idx])
            ve = float(self._ve[idx])
            l = float(self._l[idx])
            t = float(self._t[idx])
            r = float(self._r[idx])

            # Natural initial approximation:
            # assets ~= equity + present value of debt.
            asset_value_0 = e + l * np.exp(-r * t)

            # Approximate asset volatility using the equity-to-asset ratio.
            asset_volatility_0 = ve * e / asset_value_0

            x0 = np.array(
                [
                    asset_value_0,
                    asset_volatility_0,
                ]
            )

            result = optimize.root(
                _merton_equations,
                x0,
                args=(e, ve, l, t, r),
            )

            if not result.success:
                raise FinError("Unable to solve for Merton asset value and volatility: " f"{result.message}")

            asset_value = result.x[0]
            asset_volatility = result.x[1]

            if asset_value <= 0.0:
                raise FinError("Solved Merton asset value is not positive.")

            if asset_volatility <= 0.0:
                raise FinError("Solved Merton asset volatility is not positive.")

            asset_values[idx] = asset_value
            asset_volatilities[idx] = asset_volatility

        return asset_values, asset_volatilities

    def market_equity_value(self) -> np.ndarray:
        """Return the observed market equity value."""

        return self._market_equity_value

    def market_equity_volatility(self) -> np.ndarray:
        """Return the observed market equity volatility."""

        return self._market_equity_volatility

    def __repr__(self) -> str:

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string(
            "EQUITY VALUE",
            self._market_equity_value,
        )
        s += label_to_string(
            "BOND FACE",
            self._l,
        )
        s += label_to_string(
            "YEARS TO MATURITY",
            self._t,
        )
        s += label_to_string(
            "RISK FREE RATE",
            self._r,
        )
        s += label_to_string(
            "ASSET GROWTH",
            self._mu,
        )
        s += label_to_string(
            "EQUITY VOLATILITY",
            self._market_equity_volatility,
        )
        s += label_to_string(
            "IMPLIED ASSET VALUE",
            self._a,
        )
        s += label_to_string(
            "IMPLIED ASSET VOLATILITY",
            self._va,
        )

        return s
