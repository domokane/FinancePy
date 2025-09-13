# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union
import numpy as np

from ..utils.helpers import label_to_string, check_argument_types
from ..utils.math import normcdf, normcdf_vect

# TODO: Redesign this class

########################################################################################


from typing import Any

class MertonFirm:
    """Implementation of the Merton Firm Value Model according to the original
    formulation by Merton with the inputs being the asset value of the firm,
    the liabilities (bond face), the time to maturity in years, the risk-free
    rate, the asset growth rate and the asset value volatility."""

    ####################################################################################

    def __init__(
        self,
        asset_value: Union[float, list, np.ndarray],
        bond_face: Union[float, list, np.ndarray],
        years_to_maturity: Union[float, list, np.ndarray],
        risk_free_rate: Union[float, list, np.ndarray],
        asset_growth_rate: Union[float, list, np.ndarray],
        asset_volatility: Union[float, list, np.ndarray],
    ) -> None:
        """Create an object that holds all of the model parameters. These
        parameters may be vectorised."""

        check_argument_types(self.__init__, locals())

        self._a = np.array(asset_value)
        self._l = np.array(bond_face)
        self._t = np.array(years_to_maturity)
        self._r = np.array(risk_free_rate)
        self._mu = np.array(asset_growth_rate)
        self._va = np.array(asset_volatility)
        self._d = self.debt_value()
        self._e = self.equity_value()
        self._ve = self.equity_vol()

    ####################################################################################

    def leverage(self) -> np.ndarray:
        """Calculate the leverage."""

        lvg = self._a / self._l
        return lvg

    ####################################################################################

    def asset_value(self) -> np.ndarray:
        """Calculate the asset value."""

        return self._a

    ####################################################################################

    def debt_face_value(self) -> np.ndarray:
        """Calculate the asset value."""

        return self._l

    ####################################################################################

    def asset_vol(self) -> np.ndarray:
        """Return the asset volatility."""

        return self._va

    ####################################################################################

    def equity_vol(self) -> np.ndarray:
        """Calculate the equity volatility."""

        e = self.equity_value()

        lvg = self._a / self._l
        sigma_root_t = self._va * np.sqrt(self._t)

        d1 = np.log(lvg) + (self._r + 0.5 * self._va**2) * self._t
        d1 = d1 / sigma_root_t
        evol = (self._a / e) * normcdf_vect(d1) * self._va
        return evol

    ####################################################################################

    def equity_value(self) -> np.ndarray:
        """Calculate the equity value."""

        lvg = self._a / self._l
        sigma_root_t = self._va * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._va**2) * self._t
        d1 = d1 / sigma_root_t
        d2 = d1 - sigma_root_t
        evalue = self._a * normcdf_vect(d1) - self._l * np.exp(
            -self._r * self._t
        ) * normcdf_vect(d2)
        return evalue

    ####################################################################################

    def debt_value(self) -> np.ndarray:
        """Calculate the debt value"""

        lvg = self._a / self._l
        sigma_root_t = self._va * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._va**2) * self._t
        d1 = d1 / sigma_root_t
        d2 = d1 - sigma_root_t
        dvalue = self._a * normcdf_vect(-d1) + self._l * np.exp(
            -self._r * self._t
        ) * normcdf_vect(d2)
        return dvalue

    ####################################################################################

    def credit_spread(self) -> np.ndarray:
        """Calculate the credit spread from the debt value."""

        dvalue = self.debt_value()
        spd = -(1.0 / self._t) * np.log(dvalue / self._l) - self._r
        return spd

    ####################################################################################

    def prob_default(self) -> np.ndarray:
        """Calculate the default probability. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate."""

        lvg = self._a / self._l
        dd = np.log(lvg) + (self._mu - (self._va**2) / 2.0) * self._t
        dd = dd / self._va / np.sqrt(self._t)
        pd = 1.0 - normcdf_vect(dd)
        return pd

    ####################################################################################

    def dist_default(self) -> np.ndarray:
        """Calculate the distance to default. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate."""

        lvg = self._a / self._l
        dd = np.log(lvg) + (self._mu - (self._va**2) / 2.0) * self._t
        dd = dd / self._va / np.sqrt(self._t)
        return dd

    ####################################################################################

    def __repr__(self) -> str:

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ASSET VALUE", self._a)
        s += label_to_string("BOND FACE", self._l)
        s += label_to_string("YEARS TO MATURITY", self._t)
        s += label_to_string("ASSET GROWTH", self._mu)
        s += label_to_string("ASSET VOLATILITY", self._va)
        return s
