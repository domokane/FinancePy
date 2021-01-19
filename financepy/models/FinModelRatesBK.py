##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64, int64
from math import ceil

from ..finutils.FinError import FinError
from ..finutils.FinMath import accruedInterpolator
from ..market.curves.FinInterpolator import FinInterpTypes, _uinterpolate
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinGlobalTypes import FinExerciseTypes
from ..finutils.FinGlobalVariables import gSmall

interp = FinInterpTypes.FLAT_FWD_RATES.value

###############################################################################
# TODO : Calculate accrued in bond option according to accrual convention
# TODO : Convergence is unstable - investigate how to improve it
# TODO : Write a fallback for gradient based alpha using bisection
# TODO : Fix treatment of accrued interest on the option expiry date.
###############################################################################

# (c) Dominic O'Kane - December-2019
# Fergal O'Kane - search root function - 16-12-2019

###############################################################################

def optionExerciseTypesToInt(optionExerciseType):

    if optionExerciseType == FinExerciseTypes.EUROPEAN:
        return 1
    if optionExerciseType == FinExerciseTypes.BERMUDAN:
        return 2
    if optionExerciseType == FinExerciseTypes.AMERICAN:
        return 3
    else:
        raise FinError("Unknown option exercise type.")

###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def f(alpha, nm, Q, P, dX, dt, N):

    sumQZ = 0.0
    for j in range(-nm, nm+1):
        x = alpha + j*dX
        rdt = np.exp(x)*dt
        sumQZ += Q[j+N] * np.exp(-rdt)

    objFn = sumQZ - P
    return objFn

###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def fprime(alpha, nm, Q, P, dX, dt, N):

    sumQZdZ = 0.0
    for j in range(-nm, nm+1):
        x = alpha + j*dX
        rdt = np.exp(x)*dt
        sumQZdZ += Q[j+N] * np.exp(-rdt) * np.exp(x)

    deriv = -sumQZdZ*dt
    return deriv

###############################################################################
# This is the secant method which is not used as I computed the derivative of
# objective function with respect to the drift term
###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def searchRoot(x0, nm, Q, P, dX, dt, N):

    #    print("Searching for root", x0)
    max_iter = 50
    max_error = 1e-8

    x1 = x0 * 1.0001
    f0 = f(x0, nm, Q, P, dX, dt, N)
    f1 = f(x1, nm, Q, P, dX, dt, N)

    for _ in range(0, max_iter):

        df = f1 - f0

        if df == 0.0:
            raise FinError("Search for alpha fails due to zero derivative")

        x = x1 - f1 * (x1 - x0) / df
        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, nm, Q, P, dX, dt, N)

        if (abs(f1) <= max_error):
            return x1

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################
# This is Newton Raphson which is faster than the secant measure as it has the
# analytical derivative  that is easy to calculate.
###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def searchRootDeriv(x0, nm, Q, P, dX, dt, N):

    max_iter = 50
    max_error = 1e-8

    for _ in range(0, max_iter):

        fval = f(x0, nm, Q, P, dX, dt, N)

        if abs(fval) <= max_error:
            return x0

        fderiv = fprime(x0, nm, Q, P, dX, dt, N)

        if abs(fderiv) == 0.0:
            print(x0, fval, fderiv)
            raise FinError("Function derivative is zero.")

        step = fval/fderiv
        x0 = x0 - step

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################


@njit(fastmath=True, cache=True)
def bermudanSwaption_Tree_Fast(texp, tmat,
                               strikePrice, faceAmount,
                               couponTimes, couponFlows,
                               exerciseTypeInt,
                               _dfTimes, _dfValues,
                               _treeTimes, _Q,
                               _pu, _pm, _pd,
                               _rt, _dt, _a):
    ''' Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time. '''

    numTimeSteps, numNodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################
    # I shove the floating rate value into the grid so it is handled in terms
    # of PV - it will not sit on a grid date so needs to be PV adjusted.
    # This is the value of the floating leg - this is a ONE CURVE approach.
    ###########################################################################

    fixedLegFlows = np.zeros(numTimeSteps)
    floatLegValues = np.ones(numTimeSteps) # Initialise it with ones CHANGE
    numCoupons = len(couponTimes)

    # swap fixed leg flows go all the way out to the swap maturity date
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _treeTimes[n]
        df_flow = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
        fixedLegFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree
        floatLegValues[n] = strikePrice # * df_flow / df_tree

    ############################## REMOVE START ###############################

    if 1==0:
        fixedpv = 0.0
        for n in range(0, numTimeSteps-1):
            ttree = _treeTimes[n]
            df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
            flow = fixedLegFlows[n]
            pvflow = flow * df_tree
            fixedpv += pvflow
            print(">>", n, ttree, df_tree, flow, fixedpv)
        fixedpv += df_tree
        df_tree = _uinterpolate(texp, _dfTimes, _dfValues, interp)
        floatpv = df_tree
        swaptionpv = (fixedpv/df_tree - 1.0) * df_tree
        print("PV:", fixedpv, floatpv, swaptionpv)
    
        fixedpv = 0.0
        for n in range(0, numCoupons):
            tcpn = couponTimes[n]
            df = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
            flow = couponFlows[n]
            pvflow = flow * df
            fixedpv += pvflow
            print("++", n, tcpn, df, flow, fixedpv)
        fixedpv += df
        df_tree = _uinterpolate(texp, _dfTimes, _dfValues, interp)
        floatpv = df_tree
        swaptionpv = (fixedpv/df_tree - 1.0) * df_tree
        print("PV:", fixedpv, floatpv, swaptionpv)
        
    ########################### REMOVE END ####################################

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mappedTimes = np.array([0.0])
    mappedAmounts = np.array([0.0])
    for n in range(1, len(_treeTimes)):

        accdAtExpiry = 0.0
        if _treeTimes[n-1] < texp and _treeTimes[n] >= texp:
            mappedTimes = np.append(mappedTimes, texp)
            mappedAmounts = np.append(mappedAmounts, accdAtExpiry)

        if fixedLegFlows[n] > 0.0:
            mappedTimes = np.append(mappedTimes, _treeTimes[n])
            mappedAmounts = np.append(mappedAmounts, fixedLegFlows[n])

    ###########################################################################

    accrued = np.zeros(numTimeSteps)
    for m in range(0, maturityStep+1):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)
        accrued[m] *= faceAmount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if fixedLegFlows[m] > gSmall:
            accrued[m] = fixedLegFlows[m] * faceAmount

    #######################################################################

    # The value of the swap at each time and node. Principal is exchanged.
    fixedLegValues = np.zeros(shape=(numTimeSteps, numNodes))
    # The value of the option to enter into a payer swap
    payValues = np.zeros(shape=(numTimeSteps, numNodes))
    # The value of the option to enter into a receiver swap
    recValues = np.zeros(shape=(numTimeSteps, numNodes))

    # Start with the value of the fixed leg at maturity
    for k in range(0, numNodes):
        flow = 1.0 + fixedLegFlows[maturityStep]
        fixedLegValues[maturityStep, k] = flow * faceAmount

    N = jmax
    
    # Now step back to today considering early exercise on coupon dates
    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = fixedLegFlows[m] * faceAmount

        for k in range(-nm, nm+1):
            kN = k + N
            rt = _rt[m, kN]
            df = np.exp(-rt * _dt)
            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = fixedLegValues[m+1, kN]
                vm = fixedLegValues[m+1, kN-1]
                vd = fixedLegValues[m+1, kN-2]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixedLegValues[m, kN] = v
            elif k == -jmax:
                vu = fixedLegValues[m+1, kN+2]
                vm = fixedLegValues[m+1, kN+1]
                vd = fixedLegValues[m+1, kN]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixedLegValues[m, kN] = v
            else:
                vu = fixedLegValues[m+1, kN+1]
                vm = fixedLegValues[m+1, kN]
                vd = fixedLegValues[m+1, kN-1]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixedLegValues[m, kN] = v

            fixedLegValues[m, kN] += flow
            vpay = 0.0
            vrec = 0.0

            if k == jmax:
                vu = payValues[m+1, kN]
                vm = payValues[m+1, kN-1]
                vd = payValues[m+1, kN-2]
                vpay = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = payValues[m+1, kN+2]
                vm = payValues[m+1, kN+1]
                vd = payValues[m+1, kN]
                vpay = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = payValues[m+1, kN+1]
                vm = payValues[m+1, kN]
                vd = payValues[m+1, kN-1]
                vpay = (pu*vu + pm*vm + pd*vd) * df

            payValues[m, kN] = vpay

            if k == jmax:
                vu = recValues[m+1, kN]
                vm = recValues[m+1, kN-1]
                vd = recValues[m+1, kN-2]
                vrec = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = recValues[m+1, kN+2]
                vm = recValues[m+1, kN+1]
                vd = recValues[m+1, kN]
                vrec = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = recValues[m+1, kN+1]
                vm = recValues[m+1, kN]
                vd = recValues[m+1, kN-1]
                vrec = (pu*vu + pm*vm + pd*vd) * df

            recValues[m, kN] = vrec

            holdPay = payValues[m, kN]
            holdRec = recValues[m, kN]

            # The floating value is clean and so must be the fixed value
            fixedLegValue = fixedLegValues[m, kN] - accrued[m]          
            floatLegValue = floatLegValues[m]               

            payExercise = max(floatLegValue - fixedLegValue, 0.0)
            recExercise = max(fixedLegValue - floatLegValue, 0.0)

            if m == expiryStep:

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

            elif exerciseTypeInt == 2 and flow > gSmall and m > expiryStep:

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

            elif exerciseTypeInt == 3 and m > expiryStep:

                raise FinError("American optionality not allowed.")

                ## Need to define floating value on all grid dates

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

    return payValues[0, jmax], recValues[0, jmax]

###############################################################################


@njit(fastmath=True, cache=True)
def americanBondOption_Tree_Fast(texp, tmat,
                                 strikePrice, faceAmount,
                                 couponTimes, couponFlows,
                                 exerciseTypeInt,
                                 _dfTimes, _dfValues,
                                 _treeTimes, _Q,
                                 _pu, _pm, _pd,
                                 _rt,
                                 _dt, _a):
    ''' Option that can be exercised at any time over the exercise period.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. '''

    DEBUG = False

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)
    numCoupons = len(couponTimes)

    # Tree flows go all the way out to the bond maturity date
    # Do not include first coupon as it is the previous coupon and is negative
    for i in range(1, numCoupons):
        tcpn = couponTimes[i]
        
        if tcpn < 0.0:
            raise FinError("Coupon times must be positive.")

        n = int(tcpn/_dt + 0.50)
        ttree = _treeTimes[n]
        df_flow = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    #
    # mappedTimes = np.zeros(0)
    # mappedAmounts = np.zeros(0)
    # for n in range(0, len(_treeTimes)):
    #     if treeFlows[n] > 0.0:
    #         mappedTimes = np.append(mappedTimes, _treeTimes[n])
    #         mappedAmounts = np.append(mappedAmounts, treeFlows[n])
    ###########################################################################

    accrued = np.zeros(numTimeSteps)
    for m in range(0, maturityStep+1):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, couponTimes, couponFlows)
        accrued[m] *= faceAmount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > gSmall:
            accrued[m] = treeFlows[m] * faceAmount

    if DEBUG:
        for i in range(0, expiryStep+1):
            print(i, treeFlows[i], accrued[i])

    #######################################################################

    callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    bondValues = np.zeros(shape=(numTimeSteps, numNodes))

    # Start with the value of the bond at maturity
    for k in range(0, numNodes):
        bondValues[maturityStep, k] = (1.0 + treeFlows[maturityStep]) \
            * faceAmount

    if DEBUG:
        fullPrice = bondValues[maturityStep, 0]
        cleanPrice = fullPrice - accrued[maturityStep]
        print(m, _treeTimes[m], accrued[m], fullPrice, cleanPrice, 0, 0)

    # Step back from maturity to expiry date but with no exercise allowed.
    for m in range(maturityStep-1, expiryStep, -1):

        nm = min(m, jmax)
        flow = treeFlows[m] * faceAmount

        for k in range(-nm, nm+1):
            kN = k + jmax
            r = _rt[m, kN]
            df = np.exp(-r * _dt)

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

        if DEBUG:
            print(m, _treeTimes[m], accrued[m], fullPrice, cleanPrice, 0, 0)

    # Now consider exercise of the option on and before the expiry date
    for m in range(expiryStep, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * faceAmount

        for k in range(-nm, nm+1):

            kN = k + jmax
            r = _rt[m, kN]
            df = np.exp(-r * _dt)
 
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

            vcall = 0.0
            vput = 0.0

            if k == jmax:
                vu = callOptionValues[m+1, kN]
                vm = callOptionValues[m+1, kN-1]
                vd = callOptionValues[m+1, kN-2]
                vcall = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
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
            elif k == -jmax:
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

            fullPrice = bondValues[m, kN]
            cleanPrice = fullPrice - accrued[m]
            callExercise = max(cleanPrice - strikePrice, 0.0)
            putExercise = max(strikePrice - cleanPrice, 0.0)

            holdCall = callOptionValues[m, kN]
            holdPut = putOptionValues[m, kN]

            if m == expiryStep:

                callOptionValues[m, kN] = max(callExercise, holdCall)
                putOptionValues[m, kN] = max(putExercise, holdPut)

            elif exerciseTypeInt == 3 and m < expiryStep:  # AMERICAN

                callOptionValues[m, kN] = max(callExercise, holdCall)
                putOptionValues[m, kN] = max(putExercise, holdPut)

        if DEBUG:
            print(m, _treeTimes[m], accrued[m], fullPrice, cleanPrice,
                  callExercise, putExercise)
            
    return callOptionValues[0, jmax], putOptionValues[0, jmax]

###############################################################################


@njit(fastmath=True, cache=True)
def callablePuttableBond_Tree_Fast(couponTimes, couponFlows,
                                   callTimes, callPrices,
                                   putTimes, putPrices, faceAmount,
                                   _sigma, _a, _Q,  # IS SIGMA USED ?
                                   _pu, _pm, _pd, _rt, _dt, _treeTimes,
                                   _dfTimes, _dfValues):
    ''' Value a bond with embedded put and call options that can be exercised
    at any time over the specified list of put and call dates.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. '''

    #######################################################################
    numTimeSteps, numNodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * _dt))
    tmat = couponTimes[-1]
    maturityStep = int(tmat/dt + 0.50)

    ###########################################################################
    # Map coupons onto tree while preserving their present value
    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)

    numCoupons = len(couponTimes)
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _treeTimes[n]
        df_flow = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    #######################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    #######################################################################

    mappedTimes = np.array([0.0])
    mappedAmounts = np.array([0.0])
    for n in range(1, len(_treeTimes)):
        if treeFlows[n] > 0.0:
            mappedTimes = np.append(mappedTimes,_treeTimes[n])
            mappedAmounts = np.append(mappedAmounts, treeFlows[n])

    #######################################################################

    accrued = np.zeros(numTimeSteps)
    for m in range(0, numTimeSteps):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)
        accrued[m] *= faceAmount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * faceAmount

    ###########################################################################
    # map call onto tree - must have no calls at high value
    ###########################################################################

    treeCallValue = np.ones(numTimeSteps) * faceAmount * 1000.0
    numCalls = len(callTimes)
    for i in range(0, numCalls):
        callTime = callTimes[i]
        n = int(callTime/dt + 0.50)
        treeCallValue[n] = callPrices[i]

    # map puts onto tree
    treePutValue = np.zeros(numTimeSteps)
    numPuts = len(putTimes)
    for i in range(0, numPuts):
        putTime = putTimes[i]
        n = int(putTime/dt + 0.50)
        treePutValue[n] = putPrices[i]

    ###########################################################################
    # Value the bond by backward induction starting at bond maturity
    ###########################################################################

    callPutBondValues = np.zeros(shape=(numTimeSteps, numNodes))
    bondValues = np.zeros(shape=(numTimeSteps, numNodes))

    DEBUG = False
    if DEBUG:
        df = 1.0
        px = 0.0
        for i in range(0, maturityStep+1):
            flow = treeFlows[i]
            t = _treeTimes[i]
            df = _uinterpolate(t, _dfTimes, _dfValues, interp)
            px += flow*df
        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = min(maturityStep, jmax)
    vcall = treeCallValue[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * faceAmount
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm+1):
        kN = k + jmax
        bondValues[m, kN] = (1.0 + treeFlows[m]) * faceAmount
        callPutBondValues[m, kN] = value

    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * faceAmount
        vcall = treeCallValue[m]
        vput = treePutValue[m]

        for k in range(-nm, nm+1):
            kN = k + jmax
            rt = _rt[m, kN]
            df = np.exp(-rt*dt)
            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = bondValues[m+1, kN]
                vm = bondValues[m+1, kN-1]
                vd = bondValues[m+1, kN-2]
            elif k == -jmax:
                vu = bondValues[m+1, kN+2]
                vm = bondValues[m+1, kN+1]
                vd = bondValues[m+1, kN]
            else:
                vu = bondValues[m+1, kN+1]
                vm = bondValues[m+1, kN]
                vd = bondValues[m+1, kN-1]

            v = (pu*vu + pm*vm + pd*vd) * df
            bondValues[m, kN] = v
            bondValues[m, kN] += flow

            if k == jmax:
                vu = callPutBondValues[m+1, kN]
                vm = callPutBondValues[m+1, kN-1]
                vd = callPutBondValues[m+1, kN-2]
            elif k == -jmax:
                vu = callPutBondValues[m+1, kN+2]
                vm = callPutBondValues[m+1, kN+1]
                vd = callPutBondValues[m+1, kN]
            else:
                vu = callPutBondValues[m+1, kN+1]
                vm = callPutBondValues[m+1, kN]
                vd = callPutBondValues[m+1, kN-1]

            vhold = (pu*vu + pm*vm + pd*vd) * df
            # Need to make add on coupons paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            callPutBondValues[m, kN] = value

    return {'bondwithoption': callPutBondValues[0, jmax],
            'bondpure': bondValues[0, jmax]}

###############################################################################


@njit(fastmath=True)
def buildTreeFast(a, sigma, treeTimes, numTimeSteps, discountFactors):
    ''' Calibrate the tree to a term structure of interest rates. '''

    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)
    dX = sigma * np.sqrt(3.0 * dt)
    jmax = ceil(0.1835/(a * dt))

    if jmax > 1000:
        raise FinError("Jmax > 1000. Increase a or dt.")

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
        jN = j + jmax

        if j == jmax:
            pu[jN] = 7.0/6.0 + 0.50 * (ajdt*ajdt - 3.0*ajdt)
            pm[jN] = -1.0/3.0 - ajdt * ajdt + 2.0 * ajdt
            pd[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt - ajdt)
        elif j == -jmax:
            pu[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt + ajdt)
            pm[jN] = -1.0/3.0 - ajdt * ajdt - 2.0 * ajdt
            pd[jN] = 7.0/6.0 + 0.50 * (ajdt * ajdt + 3.0 * ajdt)
        else:
            pu[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt - ajdt)
            pm[jN] = 2.0/3.0 - ajdt * ajdt
            pd[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt + ajdt)

    # Arrow-Debreu array
    Q = np.zeros(shape=(numTimeSteps+2, 2*jmax+1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(numTimeSteps+1)

    # Time zero is trivial for the Arrow-Debreu price
    Q[0, jmax] = 1.0

    # Estimate short rate over first year
    r0 = -np.log(discountFactors[1])/treeTimes[1]

    # We initialise x0 with value of log of r0
    x0 = np.log(r0)

    # Big loop over time steps
    for m in range(0, numTimeSteps + 1):

        nm = min(m, jmax)

        # Need to do drift adjustment which is non-linear and so requires
        # a root search algorithm to find value of x0.

        alpha[m] = searchRootDeriv(x0, nm, Q[m], discountFactors[m+1],
                                   dX, dt, jmax)

        x0 = alpha[m]

        for j in range(-nm, nm+1):
            jN = j + jmax
            X[m, jN] = alpha[m] + j*dX
            rt[m, jN] = np.exp(X[m, jN])

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):

            jN = j + jmax
            rdt = np.exp(X[m, jN]) * dt
            z = np.exp(-rdt)

            if j == jmax:
                Q[m+1, jN] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-2] += Q[m, jN] * pd[jN] * z
            elif j == -jmax:
                Q[m+1, jN] += Q[m, jN] * pd[jN] * z
                Q[m+1, jN+1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN+2] += Q[m, jN] * pu[jN] * z
            else:
                Q[m+1, jN+1] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pd[jN] * z

    return (Q, pu, pm, pd, rt, dt)

##########################################################################


class FinModelRatesBK():

    def __init__(self, 
                 sigma: float, 
                 a: float, 
                 numTimeSteps:int=100):
        ''' Constructs the Black Karasinski rate model. The speed of mean
        reversion a and volatility are passed in. The short rate process
        is given by d(log(r)) = (theta(t) - a*log(r)) * dt  + sigma * dW '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        if a < 1e-10:
            a = 1e-10

        self._a = a
        self._sigma = sigma

        if numTimeSteps < 3:
            raise FinError("Drift fitting requires at least 3 time steps")

        self._numTimeSteps = numTimeSteps

        self._Q = None
        self._rt = None
        self._treeTimes = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discountCurve = None

###############################################################################

    def buildTree(self, tmat, dfTimes, dfValues):

        if isinstance(dfTimes, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(dfValues, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = FinInterpTypes.FLAT_FWD_RATES.value

        treeMaturity = tmat * (self._numTimeSteps+1)/self._numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, self._numTimeSteps + 2)
        self._treeTimes = treeTimes

        dfTree = np.zeros(shape=(self._numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, self._numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = _uinterpolate(t, dfTimes, dfValues, interp)

        self._dfTimes = dfTimes
        self._dfs = dfValues
        
        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = buildTreeFast(self._a, self._sigma,
                            treeTimes, self._numTimeSteps, dfTree)

        return

###############################################################################

    def bondOption(self, texp, strikePrice, faceAmount,
                   couponTimes, couponFlows, exerciseType):
        ''' Value a bond option that has European or American exercise using
        the Black-Karasinski model. The model uses a trinomial tree. '''

        exerciseTypeInt = optionExerciseTypesToInt(exerciseType)

        tmat = couponTimes[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        callValue, putValue \
            = americanBondOption_Tree_Fast(texp, tmat,
                                           strikePrice, faceAmount,
                                           couponTimes, couponFlows,
                                           exerciseTypeInt,
                                           self._dfTimes, self._dfs,
                                           self._treeTimes, self._Q,
                                           self._pu, self._pm, self._pd,
                                           self._rt,
                                           self._dt, self._a)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudanSwaption(self, texp, tmat, strikePrice, faceAmount,
                         couponTimes, couponFlows, exerciseType):
        ''' Swaption that can be exercised on specific dates over the exercise
        period. Due to non-analytical bond price we need to extend tree out to
        bond maturity and take into account cash flows through time. '''
        
        exerciseTypeInt = optionExerciseTypesToInt(exerciseType)

        tmat = couponTimes[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        payValue, recValue \
            = bermudanSwaption_Tree_Fast(texp, tmat,
                                         strikePrice, faceAmount,
                                         couponTimes, couponFlows,
                                         exerciseTypeInt,
                                         self._dfTimes, self._dfs,
                                         self._treeTimes, self._Q,
                                         self._pu, self._pm, self._pd,
                                         self._rt,
                                         self._dt, self._a)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def callablePuttableBond_Tree(self,
                                  couponTimes, couponFlows,
                                  callTimes, callPrices,
                                  putTimes, putPrices,
                                  face):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

        callTimes = np.array(callTimes)
        putTimes = np.array(putTimes)

        callPrices = np.array(callPrices)
        putPrices = np.array(putPrices)

        v = callablePuttableBond_Tree_Fast(couponTimes, couponFlows,
                                           callTimes, callPrices,
                                           putTimes, putPrices, face,
                                           self._sigma, self._a,
                                           self._Q,
                                           self._pu, self._pm, self._pd,
                                           self._rt, self._dt,
                                           self._treeTimes,
                                           self._dfTimes, self._dfs)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = "Black-Karasinski Model\n"
        s += labelToString("Sigma", self._sigma)
        s += labelToString("a", self._a)
        s += labelToString("numTimeSteps", self._numTimeSteps)
        return s

###############################################################################
