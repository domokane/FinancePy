# -*- coding: utf-8 -*-

import numpy as np
from numba import njit, float64, int64
from math import ceil, sqrt, exp

from ..finutils.FinError import FinError
from ..market.curves.FinInterpolate import FinInterpMethods, uinterpolate

###############################################################################
# TODO : Calculate accrued in bond option according to accrual convention
# TODO : Convergence is unstable - investigate how to improve it
###############################################################################

# (c) Dominic O'Kane - December-2019
# (Fergal O'Kane - search root function - 16-12-2019

###############################################################################

@njit(float64(float64,int64,float64[:],float64,float64,float64,int64),
	  fastmath=True, cache=True)
def f(alpha, nm, Q, P, dX, dt, N):

    sumQZ = 0.0
    for j in range(-nm, nm+1):
         x = alpha + j*dX
         rdt = exp(x)*dt
         sumQZ += Q[j+N] * exp(-rdt)

    objFn = sumQZ - P
    return objFn

###############################################################################

@njit(float64(float64,int64,float64[:],float64,float64,float64,int64),
	  fastmath=True, cache=True)
def fprime(alpha, nm, Q, P, dX, dt, N):

    sumQZdZ = 0.0
    for j in range(-nm, nm+1):
         x = alpha + j*dX
         rdt = exp(x)*dt
         sumQZdZ += Q[j+N] * exp(-rdt) * exp(x)

    deriv = -sumQZdZ*dt
    return deriv

###############################################################################
# This is the secant method which is not used as I computed the derivative of
# objective function with respect to the drift term

@njit(float64(float64,int64,float64[:],float64,float64,float64,int64),
	  fastmath=True, cache=True)
def searchRoot(x0, nm, Q, P, dX, dt, N):

    max_iter = 50
    max_error = 1e-7

    x1 = x0 * 1.0001
    f0 = f(x0, nm, Q, P, dX, dt, N)
    f1 = f(x1, nm, Q, P, dX, dt, N)

    for i in range(max_iter):
        x = x1 - f1 * (x1-x0)/(f1-f0)

        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, nm, Q, P, dX, dt, N)

        if (abs(f1) <= max_error):
            return x1

    return 999.0

###############################################################################
# This is Newton Raphson which is faster than the secant measure as it has the
# analytical derivative  that is easy to calculate
@njit(float64(float64,int64,float64[:],float64,float64,float64,int64),
	  fastmath=True, cache=True)
def searchRootDeriv(x0, nm, Q, P, dX, dt, N):

    max_iter = 50
    max_error = 1e-7

    for i in range(max_iter):

        fval = f(x0, nm, Q, P, dX, dt, N)

        if abs(fval) <= max_error:
            return x0

        fderiv = fprime(x0, nm, Q, P, dX, dt, N)
        step = fval/fderiv
        x0 = x0 - step

    return 999.0

###############################################################################

@njit(fastmath=True, cache = True)
def bondOptionFast(texp, tmat, strikePrice,  face, couponTimes, couponFlows,
                   americanExercise, _dfTimes, _dfValues,
                   _treeTimes, _Q, _pu, _pm, _pd, _rt, _dt, _a):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

#        print(texp, tmat, strikePrice, face)
#        print(couponTimes)
#        print(couponFlows)

        interp = FinInterpMethods.FLAT_FORWARDS.value

        #######################################################################

        numTimeSteps, numNodes = _Q.shape
        jmax = ceil(0.1835/(_a * _dt))
        N = jmax
        expiryStep = int(texp/_dt + 0.50)
        maturityStep = int(tmat/_dt + 0.50)

        #######################################################################

        treeFlows = np.zeros(numTimeSteps)
        numCoupons = len(couponTimes)

        for i in range(0, numCoupons):
            tcpn = couponTimes[i]
            n = int(round(tcpn/_dt, 0))
            ttree = _treeTimes[n]
            df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
            df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
            treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

        accrued = np.zeros(numTimeSteps)
        for m in range(0, expiryStep+1):
            treeTime = _treeTimes[m]

            for nextCpn in range(1, numCoupons):
                prevTime = couponTimes[nextCpn-1]
                nextTime = couponTimes[nextCpn]
                if treeTime > prevTime and treeTime < nextTime:
                    accdPeriod = treeTime - prevTime
                    period = (nextTime - prevTime)
                    accd = (accdPeriod / period) * couponFlows[nextCpn] * face
                    accrued[m] = accd
                    break

        #######################################################################

        callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
        putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
        bondValues = np.zeros(shape=(numTimeSteps, numNodes))

        # Start with the value of the bond at maturity
        for k in range(0, numNodes):
            bondValues[maturityStep, k] = (1.0 + treeFlows[maturityStep]) * face

        N = jmax

        for m in range(maturityStep-1, expiryStep-1, -1):

            nm = min(m, jmax)
            flow = treeFlows[m] * face

            for k in range(-nm, nm+1):
                kN = k + N
                df = exp(-_rt[m,kN]*_dt)
                pu = _pu[kN]
                pm = _pm[kN]
                pd = _pd[kN]

                if k == jmax:
                    vu = bondValues[m+1, kN]
                    vm = bondValues[m+1, kN-1]
                    vd = bondValues[m+1, kN-2]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v
                elif k == -jmax:
                    vu = bondValues[m+1, kN+2]
                    vm = bondValues[m+1, kN+1]
                    vd = bondValues[m+1, kN]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v
                else:
                    vu = bondValues[m+1, kN+1]
                    vm = bondValues[m+1, kN]
                    vd = bondValues[m+1, kN-1]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v

                bondValues[m, kN] += flow

        # Now consider exercise of the option on the expiry date
        # Start with the value of the bond at maturity and overwrite values
        nm = min(expiryStep, jmax)
        for k in range(-nm, nm+1):
            kN = k + N
            cleanPrice = bondValues[expiryStep, kN] - accrued[expiryStep]
            callOptionValues[expiryStep, kN] = max(cleanPrice-strikePrice, 0.0)
            putOptionValues[expiryStep, kN] = max(strikePrice-cleanPrice, 0.0)

        # Now step back to today considering early exercise
        for m in range(expiryStep-1, -1, -1):
            nm = min(m, jmax)
            flow = treeFlows[m] * face

            for k in range(-nm, nm+1):
                kN = k + N
                df = exp(-_rt[m,kN]*_dt)
                pu = _pu[kN]
                pm = _pm[kN]
                pd = _pd[kN]

                if k == jmax:
                    vu = bondValues[m+1, kN]
                    vm = bondValues[m+1, kN-1]
                    vd = bondValues[m+1, kN-2]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v
                elif k == jmax:
                    vu = bondValues[m+1, kN+2]
                    vm = bondValues[m+1, kN+1]
                    vd = bondValues[m+1, kN]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v
                else:
                    vu = bondValues[m+1, kN+1]
                    vm = bondValues[m+1, kN]
                    vd = bondValues[m+1, kN-1]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    bondValues[m, kN] = v

                bondValues[m, kN] += flow
                vcall = 0.0
                vput = 0.0

                if k == jmax:
                    vu = callOptionValues[m+1, kN]
                    vm = callOptionValues[m+1, kN-1]
                    vd = callOptionValues[m+1, kN-2]
                    vcall = (pu*vu + pm*vm + pd*vd) * df
                elif k == jmax:
                    vu = callOptionValues[m+1, kN+2]
                    vm = callOptionValues[m+1, kN+1]
                    vd = callOptionValues[m+1, kN]
                    vcall = (pu*vu + pm*vm + pd*vd) * df
                else:
                    vu = callOptionValues[m+1, kN+1]
                    vm = callOptionValues[m+1, kN]
                    vd = callOptionValues[m+1, kN-1]
                    vcall = (pu*vu + pm*vm + pd*vd) * df

                callOptionValues[m, kN] = vcall

                if k == jmax:
                    vu = putOptionValues[m+1, kN]
                    vm = putOptionValues[m+1, kN-1]
                    vd = putOptionValues[m+1, kN-2]
                    vput = (pu*vu + pm*vm + pd*vd) * df
                elif k == jmax:
                    vu = putOptionValues[m+1, kN+2]
                    vm = putOptionValues[m+1, kN+1]
                    vd = putOptionValues[m+1, kN]
                    vput = (pu*vu + pm*vm + pd*vd) * df
                else:
                    vu = putOptionValues[m+1, kN+1]
                    vm = putOptionValues[m+1, kN]
                    vd = putOptionValues[m+1, kN-1]
                    vput = (pu*vu + pm*vm + pd*vd) * df

                putOptionValues[m, kN] = vput

                if americanExercise is True:
                    cleanPrice = bondValues[m, kN] - accrued[m]
                    exercise = max(cleanPrice - strikePrice,0)

                    hold = callOptionValues[m, kN]
                    callOptionValues[m, kN] = max(exercise, hold)

                    hold = putOptionValues[m, kN]
                    putOptionValues[m, kN] = max(exercise, hold)

        return (callOptionValues[0, N], putOptionValues[0,N])

###############################################################################
# Unable to jit this totally due to scipty optimise
@njit(fastmath=True)
def buildTreeFast(a, sigma, treeTimes, numTimeSteps, discountFactors):

    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)
    dX = sigma * sqrt(3.0 * dt)
    jmax = ceil(0.1835/(a * dt))
    jmin = - jmax
    N = jmax

    pu = np.zeros(shape=(2*jmax+1))
    pm = np.zeros(shape=(2*jmax+1))
    pd = np.zeros(shape=(2*jmax+1))

    # The short rate goes out one step extra to have the final short rate
    # This is the BK model so x = log(r)
    X = np.zeros(shape=(numTimeSteps+2, 2*jmax+1))
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
    x0 = 3.0

    # Big loop over time steps
    for m in range(0, numTimeSteps + 1):

        nm = min(m, jmax)

        # Need to do drift adjustment which is non-linear and so requires a
        # root search algorithm - we choose to use Newton-Raphson
#        argtuple = (m, nm, Q[m], discountFactors[m+1], dX, dt, N)

        alpha[m] = searchRootDeriv(x0, nm, Q[m], discountFactors[m+1], dX, dt, N)
        x0 = alpha[m]

        for j in range(-nm, nm+1):
            jN = j + N
            X[m, jN] = alpha[m] + j*dX
            rt[m, jN] = exp(X[m, jN])

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):
            jN = j + N
            rdt = exp(X[m, jN]) * dt
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

    return (Q, pu, pm, pd, rt, dt)

##########################################################################


class FinBlackKarasinskiRateModel():

    def __init__(self, a, sigma):
        ''' Constructs the Black Karasinski rate model. The speed of mean
        reversion a and volatility are passed in. The short rate process
        is given by d(log(r)) = (theta(t) - a*log(r)) * dt  + sigma * dW '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        self._a = a
        self._sigma = sigma
        self._Q = None
        self._rt = None
        self._treeTimes = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discountCurve = None

###############################################################################

    def bondOption(self, texp, strikePrice,
                   face, couponTimes, couponFlows, americanExercise):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

        tmat = couponTimes[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        callValue, putValue \
            = bondOptionFast(texp, tmat, strikePrice, face, couponTimes,
                             couponFlows, americanExercise,
                             self._dfTimes, self._dfValues,
                             self._treeTimes, self._Q,
                             self._pu, self._pm, self._pd, self._rt,
                             self._dt, self._a)

        N = ceil(0.1835/(self._a * self._dt))

        return (callValue, putValue)

###############################################################################

    def buildTree(self, tmat, numTimeSteps, dfTimes, dfValues):

        interp = FinInterpMethods.FLAT_FORWARDS.value

        treeMaturity = tmat * (numTimeSteps+1)/numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, numTimeSteps + 2)

        dfTree = np.zeros(shape=(numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = uinterpolate(t, dfTimes, dfValues, interp)

        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = buildTreeFast(self._a, self._sigma,
                            treeTimes, numTimeSteps, dfTree)

        self._treeTimes = treeTimes
        self._dfTimes = dfTimes
        self._dfValues = dfValues

###############################################################################
