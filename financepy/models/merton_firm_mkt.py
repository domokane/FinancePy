# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union, Any

import numpy as np
from scipy import optimize

from ..utils.math import normcdf
from ..utils.helpers import label_to_string, check_argument_types
from ..utils.error import FinError
from .merton_firm import MertonFirm

########################################################################################


def _fobj(x: Any, *args: Any) -> float:
    """Find value of asset value and vol that fit equity value and vol"""

    a, v_a = x

    e = args[0]
    v_e = args[1]
    l = args[2]
    t = args[3]
    r = args[4]

    lvg = a / l
    sigma_root_t = v_a * np.sqrt(t)
    d1 = np.log(lvg) + (r + 0.5 * v_a**2) * t
    d1 = d1 / sigma_root_t
    d2 = d1 - sigma_root_t

    v_e_lhs = (a / e) * normcdf(d1) * v_a
    e_lhs = a * normcdf(d1) - l * np.exp(-r * t) * normcdf(d2)
    obj = (e - e_lhs) ** 2 + (v_e - v_e_lhs) ** 2

    return obj


########################################################################################


class MertonFirmMkt(MertonFirm):
    """
    Market Extension of the Merton Firm Model according to the original
    formulation by Merton with the inputs being the equity value of the firm,
    the liabilities (bond face), the time to maturity in years, the risk-free
    rate, the asset growth rate and the equity volatility. The asset value and
    asset volatility are computed internally by solving two non-linear
    simultaneous equations.
    """

    ####################################################################################

    def __init__(
        self,
        equity_value: Union[float, list, np.ndarray],
        bond_face: Union[float, list, np.ndarray],
        years_to_maturity: Union[float, list, np.ndarray],
        risk_free_rate: Union[float, list, np.ndarray],
        asset_growth_rate: Union[float, list, np.ndarray],
        equity_volatility: Union[float, list, np.ndarray],
    ):
        """Create an object that holds all of the model parameters. These
        parameters may be vectorised."""

        check_argument_types(self.__init__, locals())

        if isinstance(equity_value, float):
            equity_value = [equity_value]

        if isinstance(bond_face, float):
            bond_face = [bond_face]

        if isinstance(years_to_maturity, float):
            years_to_maturity = [years_to_maturity]

        if isinstance(risk_free_rate, float):
            risk_free_rate = [risk_free_rate]

        if isinstance(asset_growth_rate, float):
            asset_growth_rate = [asset_growth_rate]

        if isinstance(equity_volatility, float):
            equity_volatility = [equity_volatility]

        self._e = np.array(equity_value)
        self._l = np.array(bond_face)
        self._t = np.array(years_to_maturity)
        self._r = np.array(risk_free_rate)
        self._mu = np.array(asset_growth_rate)
        self._ve = np.array(equity_volatility)

        nmax = max(
            len(self._e),
            len(self._l),
            len(self._t),
            len(self._r),
            len(self._mu),
            len(self._ve),
        )

        if len(self._e) != nmax and len(self._e) > 1:
            raise FinError("Len e must be 1 or maximum length of arrays")

        if len(self._l) != nmax and len(self._l) > 1:
            raise FinError("Len l must be 1 or maximum length of arrays")

        if len(self._t) != nmax and len(self._t) > 1:
            raise FinError("Len T must be 1 or maximum length of arrays")

        if len(self._r) != nmax and len(self._r) > 1:
            raise FinError("Len r must be 1 or maximum length of arrays")

        if len(self._mu) != nmax and len(self._mu) > 1:
            raise FinError("Len mu must be 1 or maximum length of arrays")

        if len(self._ve) != nmax and len(self._ve) > 1:
            raise FinError("Len mu must be 1 or maximum length of arrays")

        self._nmax = nmax
        self._solve_for_asset_value_and_vol()
        self._d = self.debt_value()

    ####################################################################################

    def _solve_for_asset_value_and_vol(self) -> None:

        self._a = []
        self._va = []

        for i in range(0, self._nmax):

            argtuple = ()

            if len(self._e) == self._nmax:
                argtuple += (self._e[i],)
            else:
                argtuple += (self._e[0],)

            if len(self._ve) == self._nmax:
                argtuple += (self._ve[i],)
            else:
                argtuple += (self._ve[0],)

            if len(self._l) == self._nmax:
                argtuple += (self._l[i],)
            else:
                argtuple += (self._l[0],)

            if len(self._t) == self._nmax:
                argtuple += (self._t[i],)
            else:
                argtuple += (self._t[0],)

            if len(self._r) == self._nmax:
                argtuple += (self._r[i],)
            else:
                argtuple += (self._r[0],)

            # I initialise asset value and vol to equity value and vol
            x0 = np.array([argtuple[0], argtuple[1]])

            result = optimize.minimize(_fobj, x0, args=argtuple, tol=1e-9)

            self._a.append(result.x[0])
            self._va.append(result.x[1])

        self._a = np.array(self._a)
        self._va = np.array(self._va)

    ####################################################################################

    def __repr__(self) -> str:

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EQUITY VALUE", self._e)
        s += label_to_string("BOND FACE", self._l)
        s += label_to_string("YEARS TO MATURITY", self._t)
        s += label_to_string("ASSET GROWTH", self._mu)
        s += label_to_string("EQUITY VOLATILITY", self._ve)
        return s
