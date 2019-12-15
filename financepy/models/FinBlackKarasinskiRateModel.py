# -*- coding: utf-8 -*-

import numpy as np
from numba import njit, jit
from math import ceil, sqrt, exp, log
from scipy import optimize

from ..finutils.FinMath import N
from ..finutils.FinError import FinError
from ..finutils.FinGlobalVariables import gDaysInYear
from ..market.curves.FinInterpolate import FinInterpMethods, uinterpolate

###############################################################################
# TODO : Calculate accrued in bond option according to accrual convention
# TODO : Numba - this is not easy right now as scipy.optimise does not work but
#        could result in a significant speed up.
# TODO : Convergence is unstable - investigate how to improve it
###############################################################################

def f(alpha, *args):

    m = args[0]
    nm = args[1]
    Q = args[2]
    P = args[3]
    dX = args[4]
    dt = args[5]
    N = args[6]

    sumQZ = 0.0
    for j in range(-nm, nm+1):
         x = alpha + j*dX
         rdt = exp(x)*dt
         sumQZ += Q[m, j+N] * exp(-rdt)

    objFn = sumQZ - P
    return objFn

###############################################################################
# Unable to jit this totally due to scipty optimise
#@jit
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
    x0 = -3.0

    # Big loop over time steps
    for m in range(0, numTimeSteps + 1):

        nm = min(m, jmax)

        # Need to do drift adjustment which is non-linear and so requires a
        # root search algorithm - we choose to use Newton-Raphson
        argtuple = (m, nm, Q, discountFactors[m+1], dX, dt, N)

        alpha[m] = optimize.newton(f, x0=x0, fprime=None, args=argtuple,
                                   tol=1e-6, maxiter=50, fprime2=None)

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
        self._r = None
        self._treeTimes = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discountCurve = None

###############################################################################

    def bondOption(self, settlementDate, expiryDate, strikePrice,
                   face, bond, americanExercise):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

        interp = FinInterpMethods.FLAT_FORWARDS.value

        texp = (expiryDate - settlementDate) / gDaysInYear
        tmat = (bond._maturityDate - settlementDate) / gDaysInYear

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        dfTimes = self._discountCurve._times
        dfValues = self._discountCurve._values

        #######################################################################

        numTimeSteps, numNodes = self._Q.shape
        dt = self._dt
        jmax = ceil(0.1835/(self._a * dt))
        N = jmax
        expiryStep = int(texp/dt + 0.50)
        maturityStep = int(tmat/dt + 0.50)

        #######################################################################

        bond.calculateFlowDates(settlementDate)
        couponTimes = [0.0]
        couponFlows = [0.0]
        cpn = bond._coupon/bond._frequency
        for flowDate in bond._flowDates[1:]:
            flowTime = (flowDate - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)
        numCoupons = len(couponTimes)
        couponTimes = np.array(couponTimes)
        couponFlows = np.array(couponFlows)

        if np.any(couponTimes < 0.0):
            raise FinError("No coupon times can be before the value date.")

        if np.any(couponTimes > tmat):
            raise FinError("No coupon times can be after the maturity date.")

        treeFlows = np.zeros(numTimeSteps)

        for i in range(0, numCoupons):
            flowTime = couponTimes[i]
            n = int(round(flowTime/dt, 0))
            treeTime = self._treeTimes[n]
            df_flow = uinterpolate(flowTime, dfTimes, dfValues, interp)
            df_tree = uinterpolate(treeTime, dfTimes, dfValues, interp)
            treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

        accrued = np.zeros(numTimeSteps)
        for m in range(0, expiryStep+1):
            treeTime = self._treeTimes[m]

            for nextCpn in range(1, numCoupons):
                prevTime = couponTimes[nextCpn-1]
                nextTime = couponTimes[nextCpn]
                if treeTime > prevTime and treeTime < nextTime:
                    accdPeriod = treeTime - prevTime
                    period = (nextTime - prevTime)
                    accd = accdPeriod * cpn * face / period
                    accrued[m] = accd
                    break

        #######################################################################

        optionValues = np.zeros(shape=(numTimeSteps, numNodes))
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
                rt = self._rt[m, kN]
                df = exp(-rt*dt)
                pu = self._pu[kN]
                pm = self._pm[kN]
                pd = self._pd[kN]

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
            optionValues[expiryStep, kN] = max(cleanPrice - strikePrice, 0.0)

        # Now step back to today considering early exercise
        for m in range(expiryStep-1, -1, -1):
            nm = min(m, jmax)
            flow = treeFlows[m] * face

            for k in range(-nm, nm+1):
                kN = k + N
                rt = self._rt[m, kN]
                df = exp(-rt*dt)
                pu = self._pu[kN]
                pm = self._pm[kN]
                pd = self._pd[kN]

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

                if k == jmax:
                    vu = optionValues[m+1, kN]
                    vm = optionValues[m+1, kN-1]
                    vd = optionValues[m+1, kN-2]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    optionValues[m, kN] = v
                elif k == jmax:
                    vu = optionValues[m+1, kN+2]
                    vm = optionValues[m+1, kN+1]
                    vd = optionValues[m+1, kN]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    optionValues[m, kN] = v
                else:
                    vu = optionValues[m+1, kN+1]
                    vm = optionValues[m+1, kN]
                    vd = optionValues[m+1, kN-1]
                    v = (pu*vu + pm*vm + pd*vd) * df
                    optionValues[m, kN] = v

                if americanExercise is True:
                    cleanPrice = bondValues[m, kN] - accrued[m]
                    exercise = max(cleanPrice - strikePrice,0)
                    hold = optionValues[m, kN]
                    optionValues[m, kN] = max(exercise, hold)

        self._bondValues = bondValues
        self._optionValues = optionValues
        return optionValues[0, N], bondValues[0,N]

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

        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = buildTreeFast(self._a, self._sigma,
                            treeTimes, numTimeSteps, discountFactors)

        self._treeTimes = treeTimes
        self._discountCurve = discountCurve

###############################################################################
