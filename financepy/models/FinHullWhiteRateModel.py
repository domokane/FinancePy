# -*- coding: utf-8 -*-

import numpy as np
from numba import njit
from math import ceil, sqrt, exp, log
from ..finutils.FinMath import N
from ..finutils.FinError import FinError
from ..market.curves.FinInterpolate import FinInterpMethods, uinterpolate

from scipy import optimize

interp = FinInterpMethods.FLAT_FORWARDS.value

###############################################################################

@njit(fastmath=True, cache=True)
def accruedInterpolator(tset, couponTimes, couponAmounts):
    ''' Fast calulation of accrued interest using an Actual/Actual type of
    convention. This does not calculate according to other conventions. '''

    numCoupons = len(couponTimes)
    for i in range(1, numCoupons):
        if couponTimes[i-1] >= tset:
            denom = couponTimes[i]-couponTimes[i-1]
            accdFrac = (tset-couponTimes[i-1])/ denom
            accdCpn = accdFrac * couponAmounts[i]
            return accdCpn

###############################################################################

@njit(fastmath=True, cache=True)
def P_Fast(t, T, Rt, delta, pt, ptd, pT, _sigma, _a):
    ''' Forward discount factor as seen at some time t which may be in the
    future for payment at time T where Rt is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor. '''

    BtT = (1.0 - exp(-_a*(T-t)))/_a
    BtDelta = (1.0 - exp(-_a * delta))/_a

    term1 = log(pT/pt) - (BtT/BtDelta) * log(ptd/pt)
    term2 = (_sigma**2)*(1.0-exp(-2.0*_a*t)) * BtT * (BtT - BtDelta)/(4.0*_a)

    logAhat = term1 - term2
    BhattT = (BtT/BtDelta) * delta
    p = exp(logAhat - BhattT * Rt)
    return p

###############################################################################

@njit(fastmath=True, cache=True)
def buildTree_Fast(a, sigma, treeTimes, numTimeSteps, discountFactors):
    ''' Fast tree construction using Numba. '''
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

    return (Q, pu, pm, pd, rt, dt)

##########################################################################

@njit(fastmath=True, cache=True)
def americanBondOption_Tree_Fast(texp, strikePrice, face,
                    couponTimes, couponAmounts, americanExercise,
                    _sigma, _a, _Q, _pu, _pm, _pd, _rt, _dt, _treeTimes,
                    _dfTimes, _dfValues):
    ''' Value an option on a bond with coupons that can have European or
    American exercise. Some minor issues to do with handling coupons on
    the option expiry date need to be solved. Also this function should be
    moved out of the class so it can be sped up using NUMBA. '''

    numTimeSteps, numNodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * dt))
    N = jmax
    expiryStep = int(texp/dt + 0.50)

    #######################################################################

    if np.any(couponTimes < 0.0):
        raise FinError("No coupon times can be before the value date.")

    treeFlows = np.zeros(numTimeSteps)

    numCoupons = len(couponTimes)
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        if tcpn <= texp:
            n = int(round(tcpn/dt, 0))
            ttree = _treeTimes[n]
            df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
            df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
            treeFlows[n] += couponAmounts[i] * 1.0 * df_flow / df_tree

    accrued = np.zeros(numTimeSteps)
    for m in range(0, expiryStep+1):
        treeTime = _treeTimes[m]

        for nextCpn in range(1, numCoupons):
            prevTime = couponTimes[nextCpn-1]
            nextTime = couponTimes[nextCpn]
            if treeTime > prevTime and treeTime < nextTime:
                accdPeriod = treeTime - prevTime
                period = (nextTime - prevTime)
                accd = accdPeriod * couponAmounts[nextCpn] * face / period
                accrued[m] = accd
                break

    #######################################################################

    callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    bondValues = np.zeros(shape=(numTimeSteps, numNodes))

    ptexp = uinterpolate(texp, _dfTimes, _dfValues, interp)
    ptdelta = uinterpolate(texp+dt, _dfTimes, _dfValues, interp)

    flow = treeFlows[expiryStep] * face
    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + N
        rt = _rt[expiryStep, kN]
        bondPrice = 0.0
        for i in range(0, numCoupons):
            tflow = couponTimes[i]
            if tflow > _treeTimes[expiryStep]: # must be >
                ptflow = uinterpolate(tflow, _dfTimes, _dfValues, interp)
                zcb = P_Fast(texp, tflow, rt, dt, ptexp, ptdelta, ptflow,
                             _sigma, _a)
                bondPrice += couponAmounts[i] * face * zcb

        bondPrice += face * zcb
        bondValues[expiryStep, kN] = bondPrice + flow

    # Now consider exercise of the option on the expiry date
    # Start with the value of the bond at maturity and overwrite values
    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + N
        cleanPrice = bondValues[expiryStep, kN] - accrued[expiryStep]
        callOptionValues[expiryStep, kN] = max(cleanPrice - strikePrice, 0.0)
        putOptionValues[expiryStep, kN] = max(strikePrice - cleanPrice, 0.0)

    # Now step back to today considering early exercise
    for m in range(expiryStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face

        for k in range(-nm, nm+1):
            kN = k + N
            rt = _rt[m, kN]
            df = exp(-rt*dt)
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

    return callOptionValues[0, N], putOptionValues[0,N]

##########################################################################

def fwdFullBondPrice(rt, *args):
    ''' Price a coupon bearing bond on the option expiry date and return
    the difference from a strike price. This is used in a root search to
    find the future expiry time short rate that makes the bond price equal
    to the option strike price. It is a key step in the Jamshidian bond
    decomposition approach. The strike is a clean price. '''

    self = args[0]
    texp = args[1]
    cpnTimes = args[2]
    cpnAmounts = args[3]
    dfTimes = args[4]
    dfValues = args[5]
    strikePrice = args[6]
    face = args[7]

    dt = 0.001
    tdelta = texp + dt
    ptexp = uinterpolate(texp, dfTimes, dfValues, interp)
    ptdelta = uinterpolate(tdelta, dfTimes, dfValues, interp)

    numFlows = len(cpnTimes)
    pv = 0.0

    for i in range(1, numFlows):
        tcpn = cpnTimes[i]
        cpn = cpnAmounts[i]

        if tcpn >= texp:
            ptcpn = uinterpolate(tcpn, dfTimes, dfValues, interp)
            zcb = P_Fast(texp, tcpn, rt, dt, ptexp, ptdelta, ptcpn,
                         self._sigma, self._a)
            pv = pv + zcb * cpn

    if tcpn >= texp:
        pv = pv + zcb

    accd = accruedInterpolator(texp, cpnTimes, cpnAmounts)

    pv_clean = pv - accd
    obj = face * pv_clean - strikePrice
    return obj

##########################################################################

class FinHullWhiteRateModel():

    def __init__(self, a, sigma):
        ''' Constructs the Hull-White rate model. The speed of mean reversion
        a and volatility are passed in. The short rate process is given by
        dr = (theta(t) - ar) * dt  + sigma * dW '''

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
        self._treeBuilt = False

###############################################################################

    def optionOnZeroCouponBond(self, texp, tmat, strikePrice, face,
                               dfTimes, dfValues):
        ''' Price an option on a zero coupon bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date. '''

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        ptexp = uinterpolate(texp, dfTimes, dfValues, interp)
        ptmat = uinterpolate(tmat, dfTimes, dfValues, interp)

        sigma = self._sigma
        a = self._a

        sigmap = (sigma/a) * (1.0 - exp(-a*(tmat-texp)))
        sigmap *= sqrt((1.0-exp(-2.0*a*texp))/2.0/a)
        sigmap = abs(sigmap) + 1e-10

        h = log((face*ptmat)/(strikePrice * ptexp)) / sigmap + sigmap/2.0
        callPrice = face * ptmat * N(h) - strikePrice * ptexp * N(h-sigmap)
        putPrice = strikePrice * ptexp * N(-h+sigmap) - face * ptmat * N(-h)

        return callPrice, putPrice

###############################################################################

    def europeanBondOption_Jamshidian(self, texp, strikePrice, face, cpnTimes,
                              cpnAmounts, dfTimes, dfValues):
        ''' Valuation of a European bond option using the Jamshidian deconstruction
        of the bond into a strip of zero coupon bonds with the short rate that
        would make the bond option be at the money forward. '''

        argtuple = (self, texp, cpnTimes, cpnAmounts,
                    dfTimes, dfValues, strikePrice, face)

        x0 = 0.05
        rstar = optimize.newton(fwdFullBondPrice, x0=x0, fprime=None,
                                args=argtuple, tol=1e-8, maxiter=50,
                                fprime2=None)

        # Now we price a series of zero coupon bonds using this short rate
        dt = 0.001
        numCoupons = len(cpnTimes)

        ptexp = uinterpolate(texp, dfTimes, dfValues, interp)
        ptdelta = uinterpolate(texp+dt, dfTimes, dfValues, interp)

        callValue = 0.0
        putValue = 0.0

        for i in range(0, numCoupons):
            tcpn = cpnTimes[i]
            cpn = cpnAmounts[i]

            if tcpn >= texp:

                ptcpn = uinterpolate(tcpn, dfTimes, dfValues, interp)

                strike = P_Fast(texp, tcpn, rstar, dt, ptexp, ptdelta, ptcpn,
                                self._sigma, self._a)

                call, put = self.optionOnZeroCouponBond(texp, tcpn, strike,
                                                        1.0, dfTimes, dfValues)

                callValue += call * cpn * face
                putValue += put * cpn * face

        callValue += call * face
        putValue += put * face

        return (callValue, putValue)

###############################################################################

    def europeanBondOption_Tree(self, texp, strikePrice, face, cpnTimes,
                                cpnAmounts):
        ''' Price an option on a coupon-paying bond using tree to generate
        short rates at the expiry date and then to analytical solution of
        zero coupon bond in HW model to calculate the corresponding bond price.
        User provides bond object and option details. '''

        dt = self._dt
        tdelta = texp + dt

        ptexp = uinterpolate(texp, self._dfTimes, self._dfValues, interp)
        ptdelta = uinterpolate(tdelta, self._dfTimes, self._dfValues, interp)

        numTimeSteps, numNodes = self._Q.shape
        expiryStep = int(texp/dt+0.50)

        callPrice = 0.0
        putPrice = 0.0
        numCoupons = len(cpnTimes)

        for k in range(0, numNodes):
            q = self._Q[expiryStep, k]
            rt = self._rt[expiryStep, k]

            pv = 0.0
            for i in range(0, numCoupons):
                tcpn = cpnTimes[i]
                cpn = cpnAmounts[i]

                if tcpn >= texp:

                    ptcpn = uinterpolate(tcpn, self._dfTimes, self._dfValues,
                                         interp)

                    zcb = P_Fast(texp, tcpn, rt, dt, ptexp, ptdelta, ptcpn,
                                 self._sigma, self._a)

                    pv += cpn * zcb

            pv += zcb

            putPayoff = max(strikePrice - pv * face, 0.0)
            callPayoff = max(pv * face - strikePrice, 0.0)
            putPrice += q * putPayoff
            callPrice += q * callPayoff

        return (callPrice, putPrice)

###############################################################################

    def optionOnZeroCouponBond_Tree(self, texp, tmat, strikePrice, face):
        ''' Price an option on a zero coupon bond using a HW trinomial
        tree. The discount curve was already supplied to the tree build. '''

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

        ptexp = uinterpolate(texp, self._dfTimes, self._dfValues, interp)
        ptdelta = uinterpolate(tdelta, self._dfTimes, self._dfValues, interp)
        ptmat = uinterpolate(tmat, self._dfTimes, self._dfValues, interp)

        numTimeSteps, numNodes = self._Q.shape
        expiryStep = int(texp/dt+0.50)

        callPrice = 0.0
        putPrice = 0.0

        for k in range(0, numNodes):

            q = self._Q[expiryStep, k]

            rt = self._rt[expiryStep, k]

            zcb = P_fast(texp, tmat, rt, dt, ptexp, ptdelta, ptmat,
                         self._sigma, self._a)

            putPayoff = max(strikePrice - zcb * face, 0.0)
            callPayoff = max(zcb * face - strikePrice, 0.0)
            putPrice += q * putPayoff
            callPrice += q * callPayoff

        return (callPrice, putPrice)

###############################################################################

    def americanBondOption_Tree(self, texp, strikePrice, face,
                        couponTimes, couponAmounts, americanExercise):
        ''' Value an option on a bond with coupons that can have European or
        American exercise. Some minor issues to do with handling coupons on
        the option expiry date need to be solved. Also this function should be
        moved out of the class so it can be sped up using NUMBA. '''

        return americanBondOption_Tree_Fast(texp, strikePrice, face,
                                            couponTimes, couponAmounts,
                                            americanExercise,
                                            self._sigma, self._a,
                                            self._Q,
                                            self._pu, self._pm, self._pd,
                                            self._rt, self._dt, self._treeTimes,
                                            self._dfTimes, self._dfValues)

###############################################################################
# THIS IS THE VERSION THAT DOES NOT USE NUMBA

    def americanBondOption_Tree_OLD(self, texp, strikePrice, face,
                        couponTimes, couponAmounts, americanExercise):
        ''' Value an option on a bond with coupons that can have European or
        American exercise. Some minor issues to do with handling coupons on
        the option expiry date need to be solved. Also this function should be
        moved out of the class so it can be sped up using NUMBA. '''

        numTimeSteps, numNodes = self._Q.shape
        dt = self._dt
        jmax = ceil(0.1835/(self._a * dt))
        N = jmax
        expiryStep = int(texp/dt + 0.50)

        #######################################################################

        treeFlows = np.zeros(numTimeSteps)
        numCoupons = len(couponTimes)

        for i in range(0, numCoupons):
            tcpn = couponTimes[i]
            if tcpn <= texp:
                n = int(round(tcpn/dt, 0))
                ttree = self._treeTimes[n]
                df_flow = uinterpolate(tcpn, self._dfTimes, self._dfValues, interp)
                df_tree = uinterpolate(ttree, self._dfTimes, self._dfValues, interp)
                treeFlows[n] += couponAmounts[i] * 1.0 * df_flow / df_tree

        accrued = np.zeros(numTimeSteps)
        for m in range(0, expiryStep+1):
            treeTime = self._treeTimes[m]

            for nextCpn in range(1, numCoupons):
                prevTime = couponTimes[nextCpn-1]
                nextTime = couponTimes[nextCpn]
                if treeTime > prevTime and treeTime < nextTime:
                    accdPeriod = treeTime - prevTime
                    period = (nextTime - prevTime)
                    accd = accdPeriod * couponAmounts[nextCpn] * face / period
                    accrued[m] = accd
                    break

        #######################################################################

        callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
        putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
        bondValues = np.zeros(shape=(numTimeSteps, numNodes))

        ptexp = uinterpolate(texp, self._dfTimes, self._dfValues, interp)
        ptdelta = uinterpolate(texp+dt, self._dfTimes, self._dfValues, interp)

        flow = treeFlows[expiryStep] * face
        nm = min(expiryStep, jmax)
        for k in range(-nm, nm+1):
            kN = k + N
            rt = self._rt[expiryStep, kN]
            bondPrice = 0.0
            for i in range(0, numCoupons):
                tflow = couponTimes[i]
                if tflow > self._treeTimes[expiryStep]: # must be >
                    ptflow = uinterpolate(tflow, self._dfTimes, self._dfValues,
                                          interp)
                    zcb = self.P(texp, tflow, rt, dt, ptexp, ptdelta, ptflow)
                    bondPrice += couponAmounts[i] * face * zcb

            bondPrice += face * zcb
            bondValues[expiryStep, kN] = bondPrice + flow

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

        self._bondValues = bondValues
        self._callOptionValues = callOptionValues
        self._putOptionValues = putOptionValues

        return (callOptionValues[0, N], putOptionValues[0,N])

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

    def buildTree(self, treeMat, numTimeSteps, dfTimes, dfValues):

        treeMaturity = treeMat * (numTimeSteps+1)/numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, numTimeSteps + 2)

        dfTree = np.zeros(shape=(numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = uinterpolate(t, dfTimes, dfValues, interp)

        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = buildTree_Fast(self._a, self._sigma,
                           treeTimes, numTimeSteps, dfTree)

        self._treeTimes = treeTimes
        self._dfTimes = dfTimes
        self._dfValues = dfValues

###############################################################################
