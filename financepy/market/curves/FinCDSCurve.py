# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from math import log, ceil

import scipy.optimize as optimize
import numpy as np
from numba import njit, float64

from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinInterpolate import interpolate, FinInterpMethods

from .FinDiscountCurve import FinDiscountCurve

###############################################################################


@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def uniformToDefaultTime(u, t, v):

    if u == 0.0:
        return 99999.0

    if u == 1.0:
        return 0.0

    numPoints = len(v)

    index = 0

    for i in range(1, numPoints):
        if u <= v[i - 1] and u > v[i]:
            index = i
            break

    if index == numPoints + 1:
        t1 = t[numPoints - 1]
        q1 = v[numPoints - 1]
        t2 = t[numPoints]
        q2 = v[numPoints]
        lam = log(q1 / q2) / (t2 - t1)
        tau = t2 - log(u / q2) / lam
    else:
        t1 = t[index - 1]
        q1 = v[index - 1]
        t2 = t[index]
        q2 = v[index]
        tau = (t1 * log(q2 / u) + t2 * log(u / q1)) / log(q2 / q1)

    return tau

###############################################################################


def f(q, *args):
    self = args[0]
    valueDate = args[1]
    cds = args[2]
    numPoints = len(self._times)
    self._values[numPoints - 1] = q
    # This is important - we calibrate a curve that makes the clean PV of the
    # CDS equal to zero and so we select the second element of the value tuple
    objFn = cds.value(valueDate, self)[1]
    return objFn

###############################################################################


class FinCDSCurve(FinDiscountCurve):
    ''' Generate a survival probability curve implied by the value of CDS
    contracts given a Libor curve and an assumed recovery rate. And interpolation
    scheme for the survival probabilities is also required. '''

    def __init__(self,
                 curveDate,
                 cdsContracts,
                 liborCurve,
                 recoveryRate=0.40,
                 useCache=False,
                 interpolationMethod=FinInterpMethods.FLAT_FORWARDS):

        self._curveDate = curveDate
        self._cdsContracts = cdsContracts
        self._recoveryRate = recoveryRate
        self._liborCurve = liborCurve
        self._interpolationMethod = interpolationMethod
        self._builtOK = False

        self._times = []
        self._values = []

        if cdsContracts is not None:
            if len(cdsContracts) > 0:
                self.validate(cdsContracts)
                self.buildCurve()

###############################################################################

    def validate(self, cdsContracts):
        ''' Ensure that contracts are in increasinbg maturity. '''

        if len(cdsContracts) == 0:
            raise ValueError("No CDS contracts have been supplied.")

        maturityDate = cdsContracts[0]._maturityDate

        for cds in cdsContracts[1:]:
            if cds._maturityDate <= maturityDate:
                raise ValueError("CDS contracts not in increasing maturity.")

            maturityDate = cds._maturityDate

###############################################################################

    def survivalProbability(self, dt):
        ''' Extract the survival proibability to date dt. '''

        if isinstance(dt, FinDate):
            t = (dt - self._curveDate) / gDaysInYear
        else:
            t = dt

        if t < 0.0:
            print(t, dt, self._curveDate)
            raise ValueError("Survival Date before curve anchor date")

        if t == 0.0:
            return 1.0

#        if self._useCache == True and t < self._cachedTimeLimit:
#            iDay = int(round(t * gDaysInYear))
#            q = self._cachedDailySurvivalProbs[iDay]
#        else:
        q = interpolate(t, self._times, self._values,
                        self._interpolationMethod.value)

        return q

###############################################################################

    def df(self, t):
        ''' Extract the discount factor from the underlying Liubor curve. '''

        if isinstance(t, FinDate):
            t = (t - self._curveDate) / gDaysInYear

        return self._liborCurve.df(t)

##########################################################################

    def buildCurve(self):

        numTimes = len(self._cdsContracts)

        # we size the vectors to include time zero

        self._times = np.array([])
        self._values = np.array([])

        self._times = np.append(self._times, 0.0)
        self._values = np.append(self._values, 1.0)

        valuationDate = self._curveDate

        for i in range(0, numTimes):

            maturityDate = self._cdsContracts[i]._maturityDate

            argtuple = (self, valuationDate, self._cdsContracts[i])
            tmat = (maturityDate - valuationDate) / gDaysInYear
            q = self._values[i]

            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, q)

            optimize.newton(f, x0=q, fprime=None, args=argtuple,
                            tol=1e-7, maxiter=50, fprime2=None)

##########################################################################

    def buildCachedIssuerCurve(self):

        if self._useCache:
            self._UseCache = False

        timeCutoff = self._cachedTimeLimit
        dt = 1.0 / gDaysInYear

        if timeCutoff <= 0:
            raise Exception("Cutoff must be greater than zero")

        self.m_CachedTimeLimit = timeCutoff

        numDaysToTimeLimit = int(ceil(timeCutoff / dt) + 1)  # check +1

        self._cachedDailySurvivalProbs = np.zeros(numDaysToTimeLimit)

        t = 0.0

        for i in range(0, numDaysToTimeLimit):

            z = self.survivalProbability(self, t)
            self._cachedDailySurvivalProbs[i] = z
            t += dt

        self._useCache = True

###############################################################################

    def dump(self):

        numPoints = len(self._times)

        for i in range(0, numPoints):
            print(self._times[i], self._values[i])

###############################################################################
