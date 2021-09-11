##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..utils.math import N

from scipy import optimize

from ..utils.helpers import label_to_string, check_argument_types
from ..utils.error import FinError
from .merton_firm import MertonFirm

###############################################################################


def _fobj(x, *args):
    """ Find value of asset value and vol that fit equity value and vol """

    A, vA = x

    E = args[0]
    vE = args[1]
    L = args[2]
    t = args[3]
    r = args[4]

    lvg = A / L
    sigmaRootT = vA * np.sqrt(t)
    d1 = np.log(lvg) + (r + 0.5 * vA ** 2) * t
    d1 = d1 / sigmaRootT
    d2 = d1 - sigmaRootT

    vE_LHS = (A / E) * N(d1) * vA
    E_LHS = A * N(d1) - L * np.exp(-r * t) * N(d2)
    obj = (E - E_LHS)**2 + (vE - vE_LHS)**2

    return obj

###############################################################################


class MertonFirmMkt(MertonFirm):
    """ Market Extension of the Merton Firm Model according to the original
    formulation by Merton with the inputs being the equity value of the firm,
    the liabilities (bond face), the time to maturity in years, the risk-free
    rate, the asset growth rate and the equity volatility. The asset value and
    asset volatility are computed internally by solving two non-linear
    simultaneous equations. """

    def __init__(self,
                 equity_value: (float, list, np.ndarray),
                 bondFace: (float, list, np.ndarray),
                 timeToMaturity: (float, list, np.ndarray),
                 risk_free_rate: (float, list, np.ndarray),
                 assetGrowthRate: (float, list, np.ndarray),
                 equity_volatility: (float, list, np.ndarray)):
        """ Create an object that holds all of the model parameters. These
        parameters may be vectorised. """

        check_argument_types(self.__init__, locals())

        if isinstance(equity_value, float):
            equity_value = [equity_value]

        if isinstance(bondFace, float):
            bondFace = [bondFace]

        if isinstance(timeToMaturity, float):
            timeToMaturity = [timeToMaturity]

        if isinstance(risk_free_rate, float):
            risk_free_rate = [risk_free_rate]

        if isinstance(assetGrowthRate, float):
            assetGrowthRate = [assetGrowthRate]

        if isinstance(equity_volatility, float):
            equity_volatility = [equity_volatility]

        self._E = np.array(equity_value)
        self._L = np.array(bondFace)
        self._t = np.array(timeToMaturity)
        self._r = np.array(risk_free_rate)
        self._mu = np.array(assetGrowthRate)
        self._vE = np.array(equity_volatility)

        nmax = max(len(self._E),
                   len(self._L),
                   len(self._t),
                   len(self._r),
                   len(self._mu),
                   len(self._vE))

        if len(self._E) != nmax and len(self._E) > 1:
            raise FinError("Len E must be 1 or maximum length of arrays")

        if len(self._L) != nmax and len(self._L) > 1:
            raise FinError("Len L must be 1 or maximum length of arrays")

        if len(self._t) != nmax and len(self._t) > 1:
            raise FinError("Len T must be 1 or maximum length of arrays")

        if len(self._r) != nmax and len(self._r) > 1:
            raise FinError("Len r must be 1 or maximum length of arrays")

        if len(self._mu) != nmax and len(self._mu) > 1:
            raise FinError("Len mu must be 1 or maximum length of arrays")

        if len(self._vE) != nmax and len(self._vE) > 1:
            raise FinError("Len mu must be 1 or maximum length of arrays")

        self._nmax = nmax
        self._solve_for_asset_value_and_vol()
        self._D = self.debt_value()

###############################################################################

    def _solve_for_asset_value_and_vol(self):

        self._A = []
        self._vA = []

        for i in range(0, self._nmax):

            argtuple = ()

            if len(self._E) == self._nmax:
                argtuple += (self._E[i],)
            else:
                argtuple += (self._E[0],)

            if len(self._vE) == self._nmax:
                argtuple += (self._vE[i],)
            else:
                argtuple += (self._vE[0],)

            if len(self._L) == self._nmax:
                argtuple += (self._L[i],)
            else:
                argtuple += (self._L[0],)

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

            self._A.append(result.x[0])
            self._vA.append(result.x[1])

        self._A = np.array(self._A)
        self._vA = np.array(self._vA)


###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EQUITY VALUE", self._E)
        s += label_to_string("BOND FACE", self._L)
        s += label_to_string("YEARS TO MATURITY", self._t)
        s += label_to_string("ASSET GROWTH", self._mu)
        s += label_to_string("EQUITY VOLATILITY", self._vE)
        return s

###############################################################################
