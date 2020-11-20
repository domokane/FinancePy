##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from scipy.stats import norm
N = norm.cdf

from ..finutils.FinHelperFunctions import labelToString, checkArgumentTypes

# TODO: Redesign this class

###############################################################################


class FinModelMertonCredit():
    ''' Implementation of the Merton Credit Model according to the original
    formulation by Merton with the inputs being the asset value of the firm, 
    the liabilities (bond face), the time to maturity in years, the risk-free 
    rate, the asset growth rate and the asset value volatility. '''

    def __init__(self,
                 assetValue: (float, list, np.ndarray),
                 bondFace: (float, list,np.ndarray),
                 timeToMaturity: (float, list, np.ndarray),
                 riskFreeRate: (float, list, np.ndarray),
                 assetGrowthRate: (float, list, np.ndarray),
                 assetVolatility: (float, list, np.ndarray)):
        ''' Create an object that holds all of the model parameters. These
        parameters may be vectorised. '''

        checkArgumentTypes(self.__init__, locals())

        self._A = np.array(assetValue)
        self._L = np.array(bondFace)
        self._t = np.array(timeToMaturity)
        self._r = np.array(riskFreeRate)
        self._mu = np.array(assetGrowthRate)
        self._vA = np.array(assetVolatility)
        self._D = self.debtValue()
        self._E = self.equityValue()
        self._vE = self.equityVol()

###############################################################################

    def leverage(self):
        ''' Calculate the leverage. '''

        lvg = self._A / self._L
        return lvg

###############################################################################

    def assetValue(self):
        ''' Calculate the asset value. '''

        return self._A

###############################################################################

    def debtFacevalue(self):
        ''' Calculate the asset value. '''

        return self._L

###############################################################################

    def equityVol(self):
        ''' Calculate the equity volatility. '''

        E = self.equityValue()

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)

        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT        
        evol = (self._A / E) * N(d1) * self._vA        
        return evol

###############################################################################

    def equityValue(self):
        ''' Calculate the equity value. '''

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT
        d2 = d1 - sigmaRootT
        evalue = self._A * N(d1) - self._L * np.exp(-self._r * self._t) * N(d2)
        return evalue

###############################################################################

    def debtValue(self):
        ''' Calculate the debt value '''

        lvg = self._A / self._L
        sigmaRootT = self._vA * np.sqrt(self._t)
        d1 = np.log(lvg) + (self._r + 0.5 * self._vA ** 2) * self._t
        d1 = d1 / sigmaRootT
        d2 = d1 - sigmaRootT
        dvalue = self._A * N(-d1) + self._L * np.exp(-self._r * self._t) * N(d2)
        return dvalue

###############################################################################

    def creditSpread(self):
        ''' Calculate the credit spread from the debt value. '''

        dvalue = self.debtValue()
        spd = -(1.0 / self._t) * np.log(dvalue / self._L) - self._r
        return spd

###############################################################################

    def probDefault(self):
        ''' Calculate the default probability. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate. '''
        
        lvg = self._A / self._L
        dd = np.log(lvg) + (self._mu - (self._vA**2)/2.0) * self._t
        dd = dd / self._vA / np.sqrt(self._t)
        pd = 1.0 - N(dd)
        return pd

###############################################################################

    def distDefault(self):
        ''' Calculate the distance to default. This is not risk-neutral so it
        uses the real world drift rather than the risk-free rate. '''
        
        lvg = self._A / self._L
        dd = np.log(lvg) + (self._mu - (self._vA**2)/2.0) * self._t
        dd = dd / self._vA / np.sqrt(self._t)
        return dd

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ASSET VALUE", self._A)
        s += labelToString("BOND FACE", self._L)
        s += labelToString("YEARS TO MATURITY", self._t)
        s += labelToString("ASSET GROWTH", self._mu)
        s += labelToString("ASSET VOLATILITY", self._vA)
        return s

###############################################################################
