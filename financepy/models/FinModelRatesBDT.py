##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64, int64

from ..finutils.FinError import FinError
from ..finutils.FinMath import accruedInterpolator
from ..market.curves.FinInterpolator import FinInterpTypes, _uinterpolate
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinGlobalTypes import FinExerciseTypes
from ..finutils.FinGlobalVariables import gSmall

interp = FinInterpTypes.FLAT_FWD_RATES.value

###############################################################################
# TODO: PUT CALL PARITY IS NOT EXACTLY OBSERVED FOR BERMUDAN SWAPTIONS WHEN 
#       VOL IS TURNED UP. SMALL EFFECT. $3 OUT OF $1m. SEE TESTFINBERMUDANSWAPTIONS
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


@njit(float64(float64, int64, float64[:, :], float64[:, :],
              float64, float64, float64),
      fastmath=True, cache=True)
def f(x0, m, Q, rt, dfEnd, dt, sigma):

    # x is the middle value on the short-rate on the tree
    midm = int(m/2)
    rt[m, midm] = x0

    for i in range(midm, 0, -1):
        rt[m, i-1] = rt[m, i] * np.exp(-2.0 * sigma * np.sqrt(dt))

    for i in range(midm + 1, m + 1, 1):
        rt[m, i] = rt[m, i-1] * np.exp(2.0 * sigma * np.sqrt(dt))

    sumInner = 0.0
    for i in range(0, m+1):
        r = rt[m, i]
        nextPeriodDf = (1.0 / ((1.0 + r)**dt))
        q = Q[m, i]
        sumInner += q * nextPeriodDf

    objFn = sumInner - dfEnd
    return objFn

###############################################################################


@njit(float64(float64, int64, float64[:, :], float64[:, :], float64, float64,
              float64), fastmath=True, cache=True)
def searchRoot(x0, m, Q, rt, dfEnd, dt, sigma):

    #    print("Searching for root", x0)
    max_iter = 10
    max_error = 1e-8

    x1 = x0 * 1.0001
    f0 = f(x0, m, Q, rt, dfEnd, dt, sigma)
    f1 = f(x1, m, Q, rt, dfEnd, dt, sigma)

    for _ in range(0, max_iter):

        df = f1 - f0

        if df == 0.0:
            raise FinError("Search for alpha fails due to zero derivative")

        x = x1 - f1 * (x1-x0)/df
        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, m, Q, rt, dfEnd, dt, sigma)

        if (abs(f1) <= max_error):
            return x1

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################


@njit(fastmath=True, cache=True)
def bermudanSwaption_Tree_Fast(texp, tmat,
                               strikePrice, faceAmount,
                               couponTimes, couponFlows,
                               exerciseTypeInt,
                               _dfTimes, _dfValues,
                               _treeTimes,
                               _Q, _rt, _dt):
    ''' Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to non-analytical bond
    price we need to extend tree out to bond maturity and take into account
    cash flows through time. '''

    pu = 0.50
    pd = 0.50

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    fixedLegFlows = np.zeros(numTimeSteps)
    floatLegValues = np.zeros(numTimeSteps)
    numCoupons = len(couponTimes)

    # Tree flows go all the way out to the swap maturity date
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _treeTimes[n]
        df_flow = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
        fixedLegFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree
        floatLegValues[n] = strikePrice # * df_flow / df_tree

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

    # Now step back to today considering early exercise on coupon dates
    for m in range(maturityStep-1, -1, -1):
        nm = m
        flow = fixedLegFlows[m] * faceAmount

        for k in range(0, nm+1):
            rt = _rt[m, k]
            df = np.exp(- rt * _dt)

            vu = fixedLegValues[m+1, k+1]
            vd = fixedLegValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            fixedLegValues[m, k] = v
            fixedLegValues[m, k] += flow

            vpay = 0.0
            vrec = 0.0

            vu = payValues[m+1, k+1]
            vd = payValues[m+1, k]
            vpay = (pu*vu + pd*vd) * df
            payValues[m, k] = vpay

            vu = recValues[m+1, k+1]
            vd = recValues[m+1, k]
            vrec = (pu*vu + pd*vd) * df
            recValues[m, k] = vrec

            holdPay = payValues[m, k]
            holdRec = recValues[m, k]

            # The floating value is clean and so must be the fixed value
            fixedLegValue = fixedLegValues[m, k] - accrued[m]          
            floatLegValue = floatLegValues[m]               

            payExercise = max(floatLegValue - fixedLegValue, 0.0)
            recExercise = max(fixedLegValue - floatLegValue, 0.0)

            if m == expiryStep:

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

            elif exerciseTypeInt == 2 and flow > gSmall and m > expiryStep:

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

            elif exerciseTypeInt == 3 and m > expiryStep:

                raise FinError("American optionality not completed.")

                ## Need to define floating value on all grid dates

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

    return payValues[0, 0], recValues[0, 0]

###############################################################################


@njit(fastmath=True, cache=True)
def americanBondOption_Tree_Fast(texp, tmat,
                                 strikePrice, faceAmount,
                                 couponTimes, couponFlows,
                                 exerciseTypeInt,
                                 _dfTimes, _dfValues,
                                 _treeTimes, _Q,
                                 _rt, _dt):
    ''' Option to buy or sell bond at a specified strike price that can be
    exercised over the exercise period depending on choice of exercise type.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. '''

    DEBUG = False

    ###########################################################################

    pu = 0.50
    pd = 0.50

    if DEBUG:
        print("COUPON TIMES", couponTimes)
        print("COUPON AMOUNTS", couponFlows)

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
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
            print(couponTimes)
            raise FinError("Coupon times must be positive.")

        n = int(tcpn/_dt + 0.50)
        ttree = _treeTimes[n]
        df_flow = _uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = _uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # mappedTimes = np.zeros(0)   # CHANGE
    # mappedAmounts = np.zeros(0)  # CHANGE
    # for n in range(0, len(_treeTimes)):
    #     if treeFlows[n] > 0.0:
    #         mappedTimes = np.append(mappedTimes, _treeTimes[n])
    #         mappedAmounts = np.append(mappedAmounts, treeFlows[n])
    # if DEBUG:
    #     print("MAPPED TIMES", mappedTimes)
    #     print("MAPPED AMOUNTS", mappedAmounts)
    #
    #    if mappedTimes[0] > gSmall:
    #        raise FinError("Mapped times [0] must be <= 0 for first coupon > 0")
    #
    ###########################################################################

    if DEBUG:
        print(_treeTimes)
        print(treeFlows)

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

    ###########################################################################

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

        nm = m
        flow = treeFlows[m] * faceAmount

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r * _dt)

            vu = bondValues[m+1, k+1]
            vd = bondValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            bondValues[m, k] = v
            bondValues[m, k] += flow
            fullPrice = bondValues[m, k]
            cleanPrice = fullPrice - accrued[m]

        if DEBUG:
            print(m, _treeTimes[m], accrued[m], fullPrice, cleanPrice, 0, 0)

    # Now step back to today from the expiry date considering early exercise
    for m in range(expiryStep, -1, -1):
        nm = m
        flow = treeFlows[m] * faceAmount

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r*_dt)

            vu = bondValues[m+1, k+1]
            vd = bondValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            bondValues[m, k] = v
            bondValues[m, k] += flow

            vcall = 0.0
            vput = 0.0

            vu = callOptionValues[m+1, k+1]
            vd = callOptionValues[m+1, k]
            vcall = (pu*vu + pd*vd) * df

            callOptionValues[m, k] = vcall

            vu = putOptionValues[m+1, k+1]
            vd = putOptionValues[m+1, k]
            vput = (pu*vu + pd*vd) * df

            putOptionValues[m, k] = vput

            fullPrice = bondValues[m, k]
            cleanPrice = fullPrice - accrued[m]
            callExercise = max(cleanPrice - strikePrice, 0.0)
            putExercise = max(strikePrice - cleanPrice, 0.0)

            holdCall = callOptionValues[m, k]
            holdPut = putOptionValues[m, k]

            if m == expiryStep: 

                callOptionValues[m, k] = max(callExercise, holdCall)
                putOptionValues[m, k] = max(putExercise, holdPut)

            elif exerciseTypeInt == 3 and m < expiryStep:  

                callOptionValues[m, k] = max(callExercise, holdCall)
                putOptionValues[m, k] = max(putExercise, holdPut)

        if DEBUG:
            print(m, _treeTimes[m], accrued[m], fullPrice, cleanPrice,
                  callExercise, putExercise)

    return callOptionValues[0, 0], putOptionValues[0, 0]

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

    pu = 0.50
    pd = 0.50

    #######################################################################
    numTimeSteps, numNodes = _Q.shape
    dt = _dt
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
            mappedTimes = np.append(mappedTimes, _treeTimes[n])
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
            px += flow * df
        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = maturityStep
    vcall = treeCallValue[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * faceAmount
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(0, nm+1):
        bondValues[m, k] = (1.0 + treeFlows[m]) * faceAmount
        callPutBondValues[m, k] = value

    for m in range(maturityStep-1, -1, -1):
        nm = m
        flow = treeFlows[m] * faceAmount
        vcall = treeCallValue[m]
        vput = treePutValue[m]

        for k in range(0, nm+1):

            rt = _rt[m, k]
            df = np.exp(-rt*dt)

            vu = bondValues[m+1, k+1]
            vd = bondValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            bondValues[m, k] = v
            bondValues[m, k] += flow

            vu = callPutBondValues[m+1, k+1]
            vd = callPutBondValues[m+1, k]

            vhold = (pu*vu + pd*vd) * df
            # Need to make add on coupons paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            callPutBondValues[m, k] = value

    return {'bondwithoption': callPutBondValues[0, 0],
            'bondpure': bondValues[0, 0]}

###############################################################################
###############################################################################


@njit(cache=True, fastmath=True)
def buildTreeFast(sigma, treeTimes, numTimeSteps, discountFactors):
    # Unlike the BK and HW Trinomial trees, this Tree is packed into the lower
    # diagonal of a square matrix because of its binomial nature. This means
    # that the indexing of the arrays is different.

    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)

    # The short rate goes out one step extra to have the final short rate
    # as it follows HW code but I am not sure this is needed. EXAMINE
    # NOTE HW code uses this to have short rate at expiry so it can use
    # analytical solutions for the zero coupon bond price
    # This is the BDT model so x = log(r)

    # Some implementations use continuous compounding. Does not affect results.
    # Useful if you want to make direct comparisons of the tree.
    CONT_COMPOUNDED = True

    Q = np.zeros(shape=(numTimeSteps+2, numTimeSteps+2))
    rt = np.zeros(shape=(numTimeSteps+2, numTimeSteps+2))

    if CONT_COMPOUNDED:
        r0 = -np.log(discountFactors[1])/dt
    else:
        r0 = (1.0 / discountFactors[1] - 1.0)/dt

    rt[0, 0] = r0
    Q[0, 0] = 1.0

    if CONT_COMPOUNDED:
        Q[1, 0] = 0.50 * np.exp(-rt[0, 0] * dt)
        Q[1, 1] = 0.50 * np.exp(-rt[0, 0] * dt)
    else:
        Q[1, 0] = 0.50 / ((1.0 + rt[0, 0]) ** dt)
        Q[1, 1] = 0.50 / ((1.0 + rt[0, 0]) ** dt)

    # The short rate goes out one step extra to have the final short rate
    # as it follows HW code but I am not sure this is needed. EXAMINE
    for m in range(1, numTimeSteps+1):

        dfEnd = discountFactors[m+1]
        searchRoot(r0, m, Q, rt, dfEnd, dt, sigma)

        if CONT_COMPOUNDED:
            Q[m+1, 0] = 0.50 * Q[m, 0] * np.exp(-rt[m, 0] * dt)
        else:
            Q[m+1, 0] = 0.50 * Q[m, 0] / ((1.0 + rt[m, 0]) ** dt)

        for n in range(1, m+1):

            if CONT_COMPOUNDED:
                Q[m+1, n] = 0.50 * Q[m, n-1] * np.exp(-rt[m, n-1] * dt) \
                    + 0.50 * Q[m, n] * np.exp(-rt[m, n] * dt)
            else:
                Q[m+1, n] = 0.50 * Q[m, n-1] / ((1.0 + rt[m, n-1])**dt) \
                    + 0.50 * Q[m, n] / ((1.0 + rt[m, n])**dt)

        if CONT_COMPOUNDED:
            Q[m+1, m+1] = 0.50 * Q[m, m] * np.exp(-rt[m, m] * dt)
        else:
            Q[m+1, m+1] = 0.50 * Q[m, m] / ((1.0 + rt[m, m])**dt)

    return (Q, rt, dt)

###############################################################################


class FinModelRatesBDT():

    def __init__(self, 
                 sigma: float, 
                 numTimeSteps:int=100):
        ''' Constructs the Black-Derman-Toy rate model in the case when the
        volatility is assumed to be constant. The short rate process simplifies
        and is given by d(log(r)) = theta(t) * dt + sigma * dW. Althopugh '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        self._sigma = sigma

        if numTimeSteps < 3:
            raise FinError("Drift fitting requires at least 3 time steps.")

        self._numTimeSteps = numTimeSteps

        self._Q = None
        self._rt = None
        self._treeTimes = None
        self._pu = 0.50
        self._pd = 0.50
        self._discountCurve = None

###############################################################################

    def buildTree(self, treeMat, dfTimes, dfValues):

        if isinstance(dfTimes, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(dfValues, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = FinInterpTypes.FLAT_FWD_RATES.value

        treeMaturity = treeMat * (self._numTimeSteps+1)/self._numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, self._numTimeSteps + 2)
        self._treeTimes = treeTimes

        dfTree = np.zeros(shape=(self._numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, self._numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = _uinterpolate(t, dfTimes, dfValues, interp)

        self._dfTimes = dfTimes
        self._dfs = dfValues

        self._Q, self._rt, self._dt \
            = buildTreeFast(self._sigma,
                            treeTimes, self._numTimeSteps, dfTree)

        return

###############################################################################

    def bondOption(self, texp, strikePrice, faceAmount,
                   couponTimes, couponFlows, exerciseType):
        ''' Value a bond option that can have European or American exercise
        using the Black-Derman-Toy model. The model uses a binomial tree. '''

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
                                           self._rt,
                                           self._dt)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudanSwaption(self, texp, tmat, strike, faceAmount,
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
                                         strike, faceAmount,
                                         couponTimes, couponFlows,
                                         exerciseTypeInt,
                                         self._dfTimes, self._dfs,
                                         self._treeTimes, self._Q,
                                         self._rt,
                                         self._dt)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def callablePuttableBond_Tree(self,
                                  couponTimes, couponFlows,
                                  callTimes, callPrices,
                                  putTimes, putPrices,
                                  faceAmount):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

        callTimes = np.array(callTimes)
        putTimes = np.array(putTimes)

        callPrices = np.array(callPrices)
        putPrices = np.array(putPrices)

        v = callablePuttableBond_Tree_Fast(couponTimes, couponFlows,
                                           callTimes, callPrices,
                                           putTimes, putPrices, faceAmount,
                                           self._sigma,
                                           self._Q,
                                           self._rt, self._dt,
                                           self._treeTimes,
                                           self._dfTimes, self._dfs)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = "Black-Derman-Toy Model\n"
        s += labelToString("Sigma", self._sigma)
        s += labelToString("numTimeSteps", self._numTimeSteps)
        return s

###############################################################################
