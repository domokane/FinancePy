##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from scipy.stats import norm
N = norm.cdf

from ..finutils.FinHelperFunctions import labelToString, checkArgumentTypes

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
        self._riskFreeRate = riskFreeRate
        self._assetGrowthRate = assetGrowthRate
        self._volatility = volatility

###############################################################################

    def leverage(self):
        ''' Calculate the leverage. '''
        lvg = self._assetValue / self._bondFace
        return lvg

###############################################################################

    def equityValue(self):
        ''' Calculate the equity value. '''

        lvg = self._assetValue / self._bondFace

        d1 = np.log(lvg) +\
            (self._riskFreeRate + 0.5 * self._volatility ** 2) \
            * self._timeToMaturity

        d1 = d1 / (self._volatility * np.sqrt(self._timeToMaturity))
        d2 = d1 - self._volatility * np.sqrt(self._timeToMaturity)

        evalue = self._assetValue * N(d1) - self._bondFace * \
            np.exp(-self._riskFreeRate * self._timeToMaturity) * N(d2)

        return evalue

###############################################################################

    def debtValue(self):
        ''' Calculate the debt value '''
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
        ''' Calculate the credit spread '''
        dvalue = self.debtValue()
        spd = -(1.0 / self._timeToMaturity) * np.log(dvalue / self._bondFace)
        - self._riskFreeRate

        return spd

###############################################################################

    def probDefault(self):
        ''' Calculate the default probability. '''
        lvg = self._assetValue / self._bondFace
        dd = np.log(lvg)
        dd += (self._assetGrowthRate - (self._volatility**2)/2.0) * \
            self._timeToMaturity
        dd = dd / self._volatility / np.sqrt(self._timeToMaturity)
        pd = 1.0 - N(dd)
        return pd

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ASSET VALUE", self._assetValue)
        s += labelToString("BOND FACE", self._bondFace)
        s += labelToString("YEARS TO MATURITY", self._timeToMaturity)
        s += labelToString("ASSET GROWTH", self._assetGrowthRate)
        s += labelToString("VOLATILITY", self._volatility)
        return s

###############################################################################
