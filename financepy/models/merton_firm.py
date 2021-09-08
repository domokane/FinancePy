##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ..utils.helpers import label_to_string, check_argument_types
import numpy as np

from scipy.stats import norm
N = norm.cdf


# TODO: Redesign this class

###############################################################################


class MertonFirm():
    """ Implementation of the Merton Firm Value Model according to the original
    formulation by Merton with the inputs being the asset value of the firm,
    the liabilities (bond face), the time to maturity in years, the risk-free
    rate, the asset growth rate and the asset value volatility. """

    def __init__(self,
                 assetValue: (float, list, np.ndarray),
                 bondFace: (float, list, np.ndarray),
                 timeToMaturity: (float, list, np.ndarray),
                 risk_free_rate: (float, list, np.ndarray),
                 assetGrowthRate: (float, list, np.ndarray),
                 assetVolatility: (float, list, np.ndarray)):
        """ Create an object that holds all of the model parameters. These
        parameters may be vectorised. """

        check_argument_types(self.__init__, locals())

        self._A = np.array(assetValue)
        self._L = np.array(bondFace)
        self._t = np.array(timeToMaturity)
        self._r = np.array(risk_free_rate)
        self._mu = np.array(assetGrowthRate)
        self._vA = np.array(assetVolatility)
        self._D = self.debt_value()
        self._E = self.equity_value()
        self._vE = self.equity_vol()

###############################################################################

    def leverage(self):
        """ Calculate the leverage. """

        lvg = self._A / self._L
        return lvg

###############################################################################

    def asset_value(self):
        """ Calculate the asset value. """

        return self._A

###############################################################################

    def debt_face_value(self):
        """ Calculate the asset value. """

        return self._L

###############################################################################

    def equity_vol(self):
        """ Calculate the equity volatility. """

        E = self.equity_value()

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)

        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT
        evol = (self._A / E) * N(d1) * self._vA
        return evol

###############################################################################

    def equity_value(self):
        """ Calculate the equity value. """

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT
        d2 = d1 - sigmaRootT
        evalue = self._A * N(d1) - self._L * np.exp(-self._r * self._t) * N(d2)
        return evalue

###############################################################################

    def debt_value(self):
        """ Calculate the debt value """

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT
        d2 = d1 - sigmaRootT
        dvalue = self._A * N(-d1) + self._L * \
            np.exp(-self._r * self._t) * N(d2)
        return dvalue

###############################################################################

    def credit_spread(self):
        """ Calculate the credit spread from the debt value. """

        dvalue = self.debt_value()
        spd = -(1.0 / self._t) * np.log(dvalue / self._L) - self._r
        return spd

###############################################################################

    def prob_default(self):
        """ Calculate the default probability. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate. """

        lvg = self._A / self._L
        dd = np.log(lvg) + (self._mu - (self._vA**2)/2.0) * self._t
        dd = dd / self._vA / np.sqrt(self._t)
        pd = 1.0 - N(dd)
        return pd

###############################################################################

    def dist_default(self):
        """ Calculate the distance to default. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate. """

        lvg = self._A / self._L
        dd = np.log(lvg) + (self._mu - (self._vA**2)/2.0) * self._t
        dd = dd / self._vA / np.sqrt(self._t)
        return dd

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ASSET VALUE", self._A)
        s += label_to_string("BOND FACE", self._L)
        s += label_to_string("YEARS TO MATURITY", self._t)
        s += label_to_string("ASSET GROWTH", self._mu)
        s += label_to_string("ASSET VOLATILITY", self._vA)
        return s

###############################################################################
