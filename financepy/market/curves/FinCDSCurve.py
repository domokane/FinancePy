##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import log
import numpy as np
import scipy.optimize as optimize
from numba import njit, float64

from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...market.curves.FinInterpolate import uinterpolate, FinInterpMethods
from ...finutils.FinHelperFunctions import inputTime, inputFrequency
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinHelperFunctions import labelToString

###############################################################################


@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def uniformToDefaultTime(u, t, v):
    ''' Fast mapping of a uniform random variable to a default time given a
    survival probability curve. '''

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
    ''' Function that returns zero when the survival probability that gives a
    zero value of the CDS has been determined. '''

    self = args[0]
    valueDate = args[1]
    cds = args[2]
    numPoints = len(self._times)
    self._values[numPoints - 1] = q
    # This is important - we calibrate a curve that makes the clean PV of the
    # CDS equal to zero and so we select the second element of the value tuple
    objFn = cds.value(valueDate, self)['clean_pv']
    return objFn

###############################################################################


class FinCDSCurve():
    ''' Generate a survival probability curve implied by the value of CDS
    contracts given a Libor curve and an assumed recovery rate. A scheme for
    the interpolation of the survival probabilities is also required. '''

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
        else:
            pass # In some cases we allow None to be passed

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

    def survProb(self, dt):
        ''' Extract the survival probability to date dt. This function
        supports vectorisation. '''

        if isinstance(dt, FinDate):
            t = (dt - self._curveDate) / gDaysInYear
        elif isinstance(dt, list):
            t = np.array(dt)
        else:
            t = dt

        if np.any(t < 0.0):
            raise ValueError("Survival Date before curve anchor date")

        if isinstance(t, np.ndarray):
            n = len(t)
            qs = np.zeros(n)
            for i in range(0, n):
                qs[i] = uinterpolate(t[i], self._times, self._values,
                                     self._interpolationMethod.value)
            return qs
        else:
            q = uinterpolate(t, self._times, self._values,
                             self._interpolationMethod.value)
            return q

        return 999

###############################################################################

    def df(self, dt):
        ''' Extract the discount factor from the underlying Libor curve. This
        function supports vectorisation. '''

        if isinstance(dt, FinDate):
            t = (dt - self._curveDate) / gDaysInYear
        elif isinstance(dt, list):
            t = np.array(dt)
        else:
            t = dt

        return self._liborCurve.df(t)

###############################################################################

    def buildCurve(self):
        ''' Construct the CDS survival curve from a set of CDS contracts '''

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

##############################################################################

    def fwd(self, dt):
        ''' Calculate the instantaneous forward rate at the forward date dt
        using the numerical derivative. '''

        t = inputTime(dt, self)
        epsilon = 1e-8
        df1 = self.df(t) * self.survProb(t)
        df2 = self.df(t+epsilon) * self.survProb(t+epsilon)
        fwd = np.log(df1/df2)/dt
        return fwd

##############################################################################

    def fwdRate(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according between dates date1 and date2
        according to the specified day count convention. '''

        if date1 < self._curveDate:
            raise ValueError("Date1 before curve value date.")

        if date2 < date1:
            raise ValueError("Date2 must not be before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1, date2)
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / yearFrac
        return fwd

##############################################################################

    def zeroRate(self, dt, compoundingFreq=-1):
        ''' Calculate the zero rate to date dt in the chosen compounding
        frequency where -1 is continuous is the default. '''

        t = inputTime(dt, self)
        f = inputFrequency(compoundingFreq)
        df = self.df(t)
        q = self.survProb(t)
        dfq = df * q

        if f == 0:  # Simple interest
            zeroRate = (1.0/dfq-1.0)/t
        if f == -1:  # Continuous
            zeroRate = -np.log(dfq) / t
        else:
            zeroRate = (dfq**(-1.0/t) - 1) * f
        return zeroRate

##############################################################################

    def __repr__(self):
        ''' Print out the details of the survival probability curve. '''
        numPoints = len(self._times)
        s = labelToString("TIME","SURVIVAL_PROBABILITY")
        for i in range(0, numPoints):
            s += labelToString(self._times[i], self._values[i])

###############################################################################
