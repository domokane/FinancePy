##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from scipy.stats import norm
N = norm.cdf

from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

# TODO: Redesign this class

###############################################################################


class FinMertonCreditModel():

    def __init__(self,
                 assetValue: float,
                 bondFace: float,
                 timeToMaturity: float,
                 riskFreeRate: float,
                 assetGrowthRate: float,
                 volatility: float):

        checkArgumentTypes(self.__init__, locals())

        self._assetValue = assetValue
        self._bondFace = bondFace
        self._timeToMaturity = timeToMaturity
        self._assetGrowthRate = assetGrowthRate
        self._volatility = volatility

###############################################################################

    def leverage(self):
        lvg = self._assetValue / self._bondFace
        return lvg

###############################################################################

    def equityValue(self):

        lvg = self._assetValue / self._bondFace

        d1 = np.log(lvg) +\
            (self._riskFreeRate + 0.5 * self._volatility ** 2) * self._timeToMaturity

        d1 = d1 / (self._volatility * np.sqrt(self._timeToMaturity))
        d2 = d1 - self._volatility * np.sqrt(self._timeToMaturity)

        evalue = self._assetValue * N(d1) - self._bondFace * \
            np.exp(-self._riskFreeRate * self._timeToMaturity) * N(d2)

        return evalue

###############################################################################

    def debtValue(self):

        lvg = self._assetValue / self._bondFace
        d1 = np.log(lvg) +\
            (self._riskFreeRate + 0.5 * self._volatility ** 2) * self._timeToMaturity
        d1 = d1 / (self._volatility * np.sqrt(self._timeToMaturity))
        d2 = d1 - self._volatility * np.sqrt(self._timeToMaturity)

        dvalue = self._assetValue * N(-d1) + self._bondFace * \
            np.exp(-self._riskFreeRate * self._timeToMaturity) * N(d2)

        return dvalue

###############################################################################

    def creditSpread(self):

        dvalue = self.debtValue()
        spd = -(1.0 / self._timeToMaturity) * np.log(dvalue / self._bondFace)
        - self._riskFreeRate

###############################################################################

    def probDefault(self):

        lvg = self._assetValue / self._bondFace
        dd = np.log(lvg)
        dd += (self._assetGrowthRate - (self._volatility**2)/2.0) * self._timeToMaturity
        dd = dd / self._volatility / np.sqrt(self._timeToMaturity)
        pd = 1.0 - N(dd)
        return pd

###############################################################################
