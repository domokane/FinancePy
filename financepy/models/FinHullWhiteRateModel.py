# -*- coding: utf-8 -*-

import numpy as np
from numba import jit, njit
from math import ceil, sqrt, exp, log
from ..market.curves.FinInterpolate import interpolate, FinInterpMethods
from ..finutils.FinMath import N
from ..finutils.FinError import FinError
from ..finutils.FinGlobalVariables import gDaysInYear

###############################################################################


def printTree(array):
    n1, n2 = array.shape
    for j in range(0, n2):
        for i in range(0, n1):
            x = array[i, j]
            if x != 99999999:
                print("%7.4f" % array[i, n2-j-1], end="")
            else:
                print("%7s" % '-', end="")
        print("")

###############################################################################


@njit(fastmath=True, cache=True)
def buildTreeFast(a, sigma, treeTimes, numTimeSteps, discountFactors):

    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)
    dR = sigma * sqrt(3.0 * dt)
    jmax = ceil(0.1835/(a * dt))
    jmin = - jmax
    N = jmax

    pu = np.zeros(shape=(2*jmax+1))
    pm = np.zeros(shape=(2*jmax+1))
    pd = np.zeros(shape=(2*jmax+1))

    # The short rate goes out one step extra to have the final short rate
    rt = np.zeros(shape=(numTimeSteps+2, 2*jmax+1))

    # probabilities start at time 0 and go out to one step before T
    # Branching is simple trinomial out to time step m=1 after which
    # the top node and bottom node connect internally to two lower nodes
    # and two upper nodes respectively. The probabilities only depend on j

    for j in range(-jmax, jmax+1):
        ajdt = a*j*dt
        jN = j + N
        if j == jmax:
            pu[jN] = 7.0/6.0 + 0.50*(ajdt*ajdt - 3.0*ajdt)
            pm[jN] = -1.0/3.0 - ajdt*ajdt + 2.0*ajdt
            pd[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt - ajdt)
        elif j == -jmax:
            pu[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt + ajdt)
            pm[jN] = -1.0/3.0 - ajdt*ajdt - 2.0*ajdt
            pd[jN] = 7.0/6.0 + 0.50*(ajdt*ajdt + 3.0*ajdt)
        else:
            pu[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt - ajdt)
            pm[jN] = 2.0/3.0 - ajdt*ajdt
            pd[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt + ajdt)

    # Arrow-Debreu array
    Q = np.zeros(shape=(numTimeSteps+2, 2*N+1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(numTimeSteps+1)

    # Time zero is trivial for the Arrow-Debreu price
    Q[0, N] = 1.0

    # Big loop over time steps
    for m in range(0, numTimeSteps + 1):

        nm = min(m, jmax)
        sumQZ = 0.0
        for j in range(-nm, nm+1):
            rdt = j*dR*dt
            sumQZ += Q[m, j+N] * exp(-rdt)
        alpha[m] = log(sumQZ/discountFactors[m+1]) / dt

        for j in range(-nm, nm+1):
            jN = j + N
            rt[m, jN] = alpha[m] + j*dR

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):
            jN = j + N
            rdt = rt[m, jN] * dt
            z = exp(-rdt)

            if j == jmax:
                Q[m+1, jN] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-2] += Q[m, jN] * pd[jN] * z
            elif j == jmin:
                Q[m+1, jN] += Q[m, jN] * pd[jN] * z
                Q[m+1, jN+1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN+2] += Q[m, jN] * pu[jN] * z
            else:
                Q[m+1, jN+1] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pd[jN] * z

    return Q, rt, dt

##########################################################################


class FinHullWhiteRateModel():

    def __init__(self, a, sigma):
        ''' Constructs the Hull-White rate model. The speed of mean reversion
        a and volatility are passed in. '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        self._a = a
        self._sigma = sigma

        self._Q = None
        self._r = None
        self._treeTimes = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discountCurve = None

###############################################################################

    def P(self, t, T, Rt, delta, pt, ptd, pT):
        ''' Forward discount factor as seen at time t for payment at time T
        where R is the delta-period short rate at time t '''

        sigma = self._sigma
        a = self._a

        BtT = (1.0 - exp(-self._a*(T-t)))/self._a
        BtDelta = (1.0 - exp(-self._a * delta))/self._a

        term1 = log(pT/pt) - (BtT/BtDelta) * log(ptd/pt)
        term2 = (sigma**2)*(1.0-exp(-2.0*a*t)) * BtT * (BtT - BtDelta)/(4.0*a)

        logAhat = term1 - term2
        BhattT = (BtT/BtDelta) * delta
        p = exp(logAhat - BhattT * Rt)
        return p

###############################################################################

    def optionOnZeroCouponBond_Anal(self, settlementDate,
                                    expiryDate, maturityDate,
                                    strikePrice, face, discountCurve):
        ''' Price an option on a zero coupon bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date. '''

        texp = (expiryDate - settlementDate) / gDaysInYear
        tmat = (maturityDate - settlementDate) / gDaysInYear

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        ptexp = discountCurve.df(texp)
        ptmat = discountCurve.df(tmat)

        sigma = self._sigma
        a = self._a

        sigmap = (sigma/a) * (1.0 - exp(-a*(tmat-texp)))
        sigmap *= sqrt((1.0-exp(-2.0*a*texp))/2.0/a)
        h = log((face * ptmat)/(strikePrice * ptexp)) / sigmap + sigmap/2.0

        callPrice = face * ptmat * N(h) - strikePrice * ptexp * N(h-sigmap)
        putPrice = strikePrice * ptexp * N(-h+sigmap) - face * ptmat * N(-h)

        return callPrice, putPrice

###############################################################################

    def optionOnZeroCouponBond_Tree(self, settlementDate,
                                    expiryDate, maturityDate,
                                    strikePrice, face):

        ''' Price an option on a zero coupon bond using a HW trinomial
        tree. The discount curve was already supplied to the tree build. '''

        texp = (expiryDate - settlementDate) / gDaysInYear
        tmat = (maturityDate - settlementDate) / gDaysInYear

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        if self._treeTimes is None:
            raise FinError("Tree has not been constructed.")

        if self._treeTimes[-1] < texp:
            raise FinError("Tree expiry must be >= option expiry date.")

        dt = self._dt
        tdelta = texp + dt
        ptexp = self._discountCurve.df(texp)
        ptdelta = self._discountCurve.df(tdelta)
        ptmat = self._discountCurve.df(tmat)

        numTimeSteps, numNodes = self._Q.shape
        expiryStep = int(texp/dt+0.50)

        callPrice = 0.0
        putPrice = 0.0

        for k in range(0, numNodes):
            q = self._Q[expiryStep, k]
            rt = self._rt[expiryStep, k]
            zcb = self.P(texp, tmat, rt, dt, ptexp, ptdelta, ptmat)
            putPayoff = max(strikePrice - zcb * face, 0.0)
            callPayoff = max(zcb * face - strikePrice, 0.0)
            putPrice += q * putPayoff
            callPrice += q * callPayoff

        return (callPrice, putPrice)

###############################################################################

    def americanOption_Tree(self, settlementDate, expiryDate, maturityDate,
                            strikePrice, face):
        ''' Right now this is a European put option on a zero!! '''

        texp = (expiryDate - settlementDate) / gDaysInYear
        tmat = (maturityDate - settlementDate) / gDaysInYear

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        numTimeSteps, numNodes = self._Q.shape
        dt = self._dt
        expiryStep = int(texp/dt + 0.50)
        print("Exp Step:", expiryStep, dt, texp, numTimeSteps)

        values = np.zeros(shape=(numTimeSteps, numNodes))
        for k in range(0, numNodes):
            rt = self._rt[expiryStep, k]
            zcb = self.P(texp, tmat, rt, dt)
            putPayoff = max(strikePrice-zcb*face, 0.0)
            values[expiryStep, k] = putPayoff

        N = self._jmax
        for m in range(expiryStep-1, -1, -1):
            nm = min(m, self._jmax)
            for k in range(-nm, nm+1):
                kN = k + N
                rt = self._rt[m, kN]
                df = exp(-rt*dt)
                pu = self._pu[kN]
                pm = self._pm[kN]
                pd = self._pd[kN]

                if k == self._jmax:
                    vu = values[m+1, kN]
                    vm = values[m+1, kN-1]
                    vd = values[m+1, kN-2]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    values[m, kN] = v
                elif k == -self._jmax:
                    vu = values[m+1, kN+2]
                    vm = values[m+1, kN+1]
                    vd = values[m+1, kN]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    values[m, kN] = v
                else:
                    vu = values[m+1, kN+1]
                    vm = values[m+1, kN]
                    vd = values[m+1, kN-1]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    values[m, kN] = v

        self._values = values
        return values[0, N]

###############################################################################

    def df_Tree(self, tmat):
        ''' Discount factor as seen from now to time tmat as long as the time
        is on the tree grid. '''

        if tmat == 0.0:
            return 1.0

        numTimeSteps, numNodes = self._Q.shape
        fn1 = tmat/self._dt
        fn2 = float(int(tmat/self._dt))
        if abs(fn1 - fn2) > 1e-6:
            raise FinError("Time not on tree time grid")

        timeStep = int(tmat / self._dt) + 1

        p = 0.0
        for i in range(0, numNodes):
            ad = self._Q[timeStep, i]
            p += ad
        zeroRate = -log(p)/tmat
        return p, zeroRate

###############################################################################

    def buildTree(self, startDate, endDate, numTimeSteps, discountCurve):

        maturity = (endDate - startDate) / gDaysInYear
        treeMaturity = maturity * (numTimeSteps+1)/numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, numTimeSteps + 2)

        discountFactors = np.zeros(shape=(numTimeSteps+2))
        discountFactors[0] = 1.0

        for i in range(1, numTimeSteps+2):
            t = treeTimes[i]
            discountFactors[i] = discountCurve.df(t)

        self._Q, self._rt, self._dt = buildTreeFast(self._a, self._sigma,
                                                    treeTimes, numTimeSteps,
                                                    discountFactors)
        self._treeTimes = treeTimes
        self._discountCurve = discountCurve

###############################################################################
