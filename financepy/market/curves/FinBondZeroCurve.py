# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from ...finutils.FinMath import scale
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinInterpolate import FinInterpMethods
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################


def f(df, *args):
    curve = args[0]
    valueDate = args[1]
    bond = args[2]
    marketCleanPrice = args[3]
    numPoints = len(curve._times)
    curve._values[numPoints - 1] = df
    bondDiscountPrice = bond.cleanPriceFromDiscountCurve(valueDate, curve)
    objFn = bondDiscountPrice - marketCleanPrice
    return objFn

###############################################################################


class FinBondZeroCurve(FinDiscountCurve):
    ''' Class to do fitting of the bond zero rate curve. '''

    def __init__(self, settlementDate, bonds, cleanPrices,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Fit the curve to a set of bond yields using the type of curve 
        specified. Weights can be provided to deal with illiquid bonds. '''

        if len(bonds) != len(cleanPrices):
            raise ValueError("Num bonds does not equal number of prices.")

        self._settlementDate = settlementDate
        self._curveDate = settlementDate
        self._bonds = bonds
        self._cleanPrices = np.array(cleanPrices)
        self._discountCurve = None
        self._interpMethod = interpMethod

        tmats = []
        for bond in self._bonds:
            tmat = (bond._maturityDate-self._settlementDate)/gDaysInYear
            tmats.append(tmat)
        self._yearsToMaturity = np.array(tmats)

        self.bootstrapZeroRates()

###############################################################################

    def bootstrapZeroRates(self):

        self._times = np.array([0.0])
        self._values = np.array([1.0])

        df = 1.0

        for i in range(0, len(self._bonds)):
            bond = self._bonds[i]
            maturityDate = bond._maturityDate
            cleanPrice = self._cleanPrices[i]
            tmat = (maturityDate - self._settlementDate) / gDaysInYear
            argtuple = (self, self._settlementDate, bond, cleanPrice)
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, df)

            optimize.newton(f, x0=df, fprime=None, args=argtuple,
                            tol=1e-8, maxiter=100, fprime2=None)

##########################################################################

    def display(self, title):
        ''' Display yield curve. '''

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel('Time to Maturity (years)')
        plt.ylabel('Zero Rate (%)')

        tmax = np.max(self._yearsToMaturity)
        t = np.linspace(0.0, int(tmax+0.5), 100)

        zeroRate = self.zeroRate(t)
        zeroRate = scale(zeroRate, 100.0)
        plt.plot(t, zeroRate, label="Zero Rate Bootstrap")
        plt.legend(loc='lower right')
        plt.ylim((min(zeroRate)-0.3, max(zeroRate)*1.1))
        plt.grid(True)
