##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import jit, njit, float64, int64

from ..finutils.FinError import FinError
from ..finutils.FinMath import accruedInterpolator
from ..market.curves.FinInterpolate import FinInterpTypes, uinterpolate
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinOptionTypes import FinOptionExerciseTypes

interp = FinInterpTypes.FLAT_FORWARDS.value

###############################################################################


def optionExerciseTypesToInt(optionExerciseType):

    if optionExerciseType == FinOptionExerciseTypes.EUROPEAN:
        return 1
    if optionExerciseType == FinOptionExerciseTypes.BERMUDAN:
        return 2
    if optionExerciseType == FinOptionExerciseTypes.AMERICAN:
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

    for i in range(0, max_iter):

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


@jit(fastmath=True, cache=True)
def bermudanSwaption_Tree_Fast(texp, tmat, strikePrice,  face,
                               couponTimes, couponFlows,
                               exerciseType,
                               _dfTimes, _dfValues,
                               _treeTimes,
                               _Q, _rt, _dt):
    ''' Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to non-analytical bond
    price we need to extend tree out to bond maturity and take into account
    cash flows through time. '''

#    DEBUG = False
    pu = 0.50
    pd = 0.50

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)
    numCoupons = len(couponTimes)

    # if DEBUG:
    #     filename = "swaptionBDT_" + str(numTimeSteps) + ".txt"
    #     outfile = open(filename, "w")

    # Tree flows go all the way out to the bond maturity date
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(round(tcpn/_dt, 0))
        ttree = _treeTimes[n]
        df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    # if DEBUG:
    #     for i in range(0, len(_treeTimes)):
    #         t = _treeTimes[i]
    #         f = treeFlows[i]
    #         outfile.write("t: %9.5f  flow: %9.5f\n" % (t, f))

    ###########################################################################
    accrued = np.zeros(numTimeSteps)

    mappedTimes = [0.0]
    mappedAmounts = [0.0]
    for n in range(1, len(_treeTimes)):

        # I AM ENFORCING ZERO ACCRUED ON THE EXPIRY DATE OF THE OPTION
        # THIS NEEDS FURTHER WORK. I SHOULD CALCULATE THE EXPIRY ACCD EXACTLY
        # AND SET IT EQUAL TO THAT AMOUNT - !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        accdAtExpiry = 0.0
        if _treeTimes[n-1] < texp and _treeTimes[n] >= texp:
            mappedTimes.append(texp)
            mappedAmounts.append(accdAtExpiry)

        if treeFlows[n] > 0.0:
            mappedTimes.append(_treeTimes[n])
            mappedAmounts.append(treeFlows[n])

    for m in range(0, maturityStep+1):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)

        # if DEBUG:
        #     print("==>", m, ttree, accrued[m])

        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face

    # if DEBUG:
    #     print("TEXP", texp)
    #     print("TREE TIMES", _treeTimes)
    #     print("TREE FLOWS", treeFlows)
    #     print("MAPPED TIMES:", mappedTimes)
    #     print("MAPPED AMOUNTS:", mappedAmounts)
    #     print("ACCD TIMES", _treeTimes)
    #     print("ACCD AMOUNT", accrued)

    #######################################################################

    payValues = np.zeros(shape=(numTimeSteps, numNodes))
    recValues = np.zeros(shape=(numTimeSteps, numNodes))
    swapFixedLegValues = np.zeros(shape=(numTimeSteps, numNodes))
    cleanPrice = 0.0

    # Start with the value of the bond at maturity
    for k in range(0, numNodes-1):
        swapFixedLegValues[maturityStep, k] = \
            (1.0 + treeFlows[maturityStep]) * face

    # if DEBUG:
    #     outfile.write("m: %d k: %d EXPIRY: %9.5f\n" %
    #                   (maturityStep, k, swapFixedLegValues[maturityStep, k]))

    # Now step back to today considering early exercise on coupon dates
    for m in range(maturityStep-1, -1, -1):

        nm = m
        flow = treeFlows[m] * face

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(- r * _dt)

            vu = swapFixedLegValues[m+1, k+1]
            vd = swapFixedLegValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            swapFixedLegValues[m, k] = v
            swapFixedLegValues[m, k] += flow
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

            cleanPrice = swapFixedLegValues[m, k] - accrued[m]
            payExercise = max(strikePrice - cleanPrice, 0.0)
            recExercise = max(cleanPrice - strikePrice, 0.0)

            holdPay = payValues[m, k]
            holdRec = recValues[m, k]

            if exerciseType == 1 and m == expiryStep:
                # European style on expiry date

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

            elif exerciseType == 2 and m == expiryStep:
                # Bermudan style on expiry date (which has no coupon)

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

            elif exerciseType == 2 and flow > 0.0 and m >= expiryStep:
                # Bermudan style on coupon dates after expiry date

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

            elif exerciseType == 3 and m >= expiryStep:
                # American style on all dates after expiry date

                payValues[m, k] = max(payExercise, holdPay)
                recValues[m, k] = max(recExercise, holdRec)

        # if DEBUG:
        #     outfile.write("m: %5d %5d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n" % \
        #                   (m,
        #                    k,
        #                    _treeTimes[m],
        #                    treeFlows[m],
        #                    swapFixedLegValues[m, k],
        #                    cleanPrice,
        #                    accrued[m],
        #                    payValues[m, k],
        #                    recValues[m, k]))
                
    # if DEBUG:
    #     outfile.close()

    return payValues[0, 0], recValues[0, 0]

###############################################################################


@njit(fastmath=True, cache=True)
def americanBondOption_Tree_Fast(texp, tmat, strikePrice,  face,
                                 couponTimes, couponFlows,
                                 exerciseType,
                                 _dfTimes, _dfValues,
                                 _treeTimes, _Q, _rt, _dt):
    ''' Option that can be exercised at any time over the exercise period.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. '''

#    DEBUG = False
    pu = 0.50
    pd = 0.50

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)
    numCoupons = len(couponTimes)

    # if DEBUG:
    #     filename = "bondOptBDT_" + str(numTimeSteps) + ".txt"
    #     outfile = open(filename, "w")

    # Tree flows go all the way out to the bond maturity date
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(round(tcpn/_dt, 0))
        ttree = _treeTimes[n]
        df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    # if DEBUG:
    #     for i in range(0, len(_treeTimes)):
    #         t = _treeTimes[i]
    #         f = treeFlows[i]
    #         outfile.write("t: %9.5f  flow: %9.5f\n" % (t, f))

    ###########################################################################
    accrued = np.zeros(numTimeSteps)

    mappedTimes = [0.0]
    mappedAmounts = [0.0]
    for n in range(1, len(_treeTimes)):

        # I AM ENFORCING ZERO ACCRUED ON THE EXPIRY DATE OF THE OPTION
        # THIS NEEDS FURTHER WORK. I SHOULD CALCULATE THE EXPIRY ACCD EXACTLY
        # AND SET IT EQUAL TO THAT AMOUNT - !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        accdAtExpiry = 0.0
        if _treeTimes[n-1] < texp and _treeTimes[n] >= texp:
            mappedTimes.append(texp)
            mappedAmounts.append(accdAtExpiry)

        if treeFlows[n] > 0.0:
            mappedTimes.append(_treeTimes[n])
            mappedAmounts.append(treeFlows[n])

    for m in range(0, maturityStep+1):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)
#        print("==>", m, ttree, accrued[m])
        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face

    # if DEBUG:
    #     print("TEXP", texp)
    #     print("TREE TIMES", _treeTimes)
    #     print("TREE FLOWS", treeFlows)
    #     print("MAPPED TIMES:", mappedTimes)
    #     print("MAPPED AMOUNTS:", mappedAmounts)
    #     print("ACCD TIMES", _treeTimes)
    #     print("ACCD AMOUNT", accrued)

    #######################################################################

    callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    bondValues = np.zeros(shape=(numTimeSteps, numNodes))

    # Start with the value of the bond at maturity
    for k in range(0, numNodes-1):
        bondValues[maturityStep, k] = (1.0 + treeFlows[maturityStep]) * face

    # if DEBUG:
    #     outfile.write("m: %d k: %d EXPIRY: %9.5f\n" %
    #                   (maturityStep, k, bondValues[maturityStep, k]))

    for m in range(maturityStep-1, expiryStep-1, -1):

        nm = m
        flow = treeFlows[m] * face

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r * _dt)

            vu = bondValues[m+1, k+1]
            vd = bondValues[m+1, k]
            v = (pu*vu + pd*vd) * df

            bondValues[m, k] = v
            bondValues[m, k] += flow

        # if DEBUG:
        #     outfile.write("m: %d k: %d flow: %9.5f BOND VALUE: %9.5f\n" %
        #                   (m, k, flow, bondValues[m, k]))

    # Now consider exercise of the option on the expiry date
    # Start with the value of the bond at maturity and overwrite values
    nm = expiryStep
    for k in range(0, nm+1):
        accd = accrued[expiryStep]
        cleanPrice = bondValues[expiryStep, k] - accd
        callOptionValues[expiryStep, k] = max(cleanPrice - strikePrice, 0.0)
        putOptionValues[expiryStep, k] = max(strikePrice - cleanPrice, 0.0)

    # if DEBUG:
    #     outfile.write("k: %d cleanPrice: %9.5f accd: %9.5f \n" % (k, cleanPrice, accd))

    # if DEBUG:
    #     outfile.write("OPTION    m    k     treeT      treeF  bondVal   cleanPx     Accd    CallVal    PutVal\n")

    # Now step back to today considering early exercise
    for m in range(expiryStep-1, -1, -1):
        nm = m
        flow = treeFlows[m] * face

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r * _dt)

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

            cleanPrice = bondValues[m, k] - accrued[m]
            callExercise = max(cleanPrice - strikePrice, 0.0)
            putExercise = max(strikePrice - cleanPrice, 0.0)

            holdCall = callOptionValues[m, k]
            holdPut = putOptionValues[m, k]

            if exerciseType == 1:  # EUROPEAN SO DO NOTHING

                callOptionValues[m, k] = holdCall
                putOptionValues[m, k] = holdPut

            elif exerciseType == 2 and flow > 0.0:  # BERMUDAN - ON CPN DATES

                callOptionValues[m, k] = max(callExercise, holdCall)
                putOptionValues[m, k] = max(putExercise, holdPut)

            elif exerciseType == 3:  # AMERICAN

                callOptionValues[m, k] = max(callExercise, holdCall)
                putOptionValues[m, k] = max(putExercise, holdPut)

        # if DEBUG:
        #     outfile.write("CALL: %5d %5d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n" % \
        #             (m, k, _treeTimes[m],
        #              treeFlows[m],
        #              bondValues[m, k],\
        #              cleanPrice, accrued[m],\
        #              callOptionValues[m, k], putOptionValues[m, k]))

    # if DEBUG:
    #     outfile.close()

    return callOptionValues[0, 0], putOptionValues[0, 0]

###############################################################################


@njit(fastmath=True, cache=True)
def callablePuttableBond_Tree_Fast(couponTimes, couponFlows,
                                   callTimes, callPrices,
                                   putTimes, putPrices, face,
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
        n = int(round(tcpn/_dt, 0))
        ttree = _treeTimes[n]
        df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    #######################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    #######################################################################

    mappedTimes = [0.0]
    mappedAmounts = [0.0]
    for n in range(1, len(_treeTimes)):
        if treeFlows[n] > 0.0:
            mappedTimes.append(_treeTimes[n])
            mappedAmounts.append(treeFlows[n])

    #######################################################################

    accrued = np.zeros(numTimeSteps)
    for m in range(0, numTimeSteps):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)
        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face

    ###########################################################################
    # map call onto tree - must have no calls at high value
    ###########################################################################

    treeCallValue = np.ones(numTimeSteps) * face * 1000.0
    numCalls = len(callTimes)
    for i in range(0, numCalls):
        callTime = callTimes[i]
        n = int(round(callTime/dt, 0))
        treeCallValue[n] = callPrices[i]

    # map puts onto tree
    treePutValue = np.zeros(numTimeSteps)
    numPuts = len(putTimes)
    for i in range(0, numPuts):
        putTime = putTimes[i]
        n = int(round(putTime/dt, 0))
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
            df = uinterpolate(t, _dfTimes, _dfValues, interp)
            px += flow*df
        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = maturityStep
    vcall = treeCallValue[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * face
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(0, nm+1):
        bondValues[m, k] = (1.0 + treeFlows[m]) * face
        callPutBondValues[m, k] = value

    for m in range(maturityStep-1, -1, -1):
        nm = m
        flow = treeFlows[m] * face
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


@njit(fastmath=True)
def buildTreeFast(sigma, treeTimes, numTimeSteps, discountFactors):
    # Unlike the BK and HK Trinomial trees, this Tree is packed into the lower
    # diagonal of a square matrix because of its binomial nature. This means
    # that the indexing of the arrays is different.

    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)

    # The short rate goes out one step extra to have the final short rate
    # as it follows HW code but I am not sure this is needed. EXAMINE
    # NOTE HW code uses this to have short rate at expiry so it can use
    # analytical solutions for the zero coupon bond price
    # This is the BDT model so x = log(r)

    # I AM USING A FLAG BUT THIS SHOULD BECOME PERMANENT AS ALL OTHER FUNCTIONS
    # EXPECT THE SHORT RATE TO BE CONTINUOUSLY COMPOUNDED
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
        Q[1, 0] = 0.50/((1.0 + rt[0, 0]) ** dt)
        Q[1, 1] = 0.50/((1.0 + rt[0, 0]) ** dt)

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
            Q[m+1, m+1] = 0.50 * Q[m, m] * np.exp(- rt[m, m] * dt)
        else:
            Q[m+1, m+1] = 0.50 * Q[m, m] / ((1.0 + rt[m, m])**dt)

    return (Q, rt, dt)

###############################################################################


class FinModelRatesBDT():

    def __init__(self, sigma, numTimeSteps=100):
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

        return

###############################################################################

    def buildTree(self, treeMat, dfTimes, dfValues):

        if isinstance(dfTimes, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(dfValues, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = FinInterpTypes.FLAT_FORWARDS.value

        treeMaturity = treeMat * (self._numTimeSteps+1)/self._numTimeSteps
        treeTimes = np.linspace(0.0, treeMaturity, self._numTimeSteps + 2)
        self._treeTimes = treeTimes

        dfTree = np.zeros(shape=(self._numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, self._numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = uinterpolate(t, dfTimes, dfValues, interp)

        self._dfTimes = dfTimes
        self._dfValues = dfValues

        self._Q, self._rt, self._dt \
            = buildTreeFast(self._sigma,
                            treeTimes, self._numTimeSteps, dfTree)

        return

###############################################################################

    def bondOption(self, texp, strike,
                   face, couponTimes, couponFlows, exerciseType):
        ''' Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. '''

        exerciseTypeInt = optionExerciseTypesToInt(exerciseType)

        tmat = couponTimes[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        callValue, putValue \
            = americanBondOption_Tree_Fast(texp, tmat, strike, face,
                                           couponTimes, couponFlows,
                                           exerciseTypeInt,
                                           self._dfTimes, self._dfValues,
                                           self._treeTimes, self._Q,
                                           self._rt,
                                           self._dt)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudanSwaption(self,
                         texp,
                         tmat,
                         strike,
                         face,
                         couponTimes,
                         couponFlows,
                         exerciseType):
        ''' Swaption that can be exercised on specific dates over the exercise
        period. Due to non-analytical bond price we need to extend tree out to
        bond maturity and take into account cash flows through time. '''

        exerciseTypeInt = optionExerciseTypesToInt(exerciseType)

        if 1 == 0:
            print(texp, tmat, strike, face)
            print(couponTimes)
            print(couponFlows)
            print(exerciseType)

        if tmat != couponTimes[-1]:
            raise FinError("Tmat does not equal couponTimes[-1]")

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        payValue, recValue \
            = bermudanSwaption_Tree_Fast(texp, tmat, strike, face,
                                         couponTimes, couponFlows,
                                         exerciseTypeInt,
                                         self._dfTimes, self._dfValues,
                                         self._treeTimes, self._Q,
                                         self._rt,
                                         self._dt)

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
                                           self._sigma,
                                           self._Q,
                                           self._rt, self._dt,
                                           self._treeTimes,
                                           self._dfTimes, self._dfValues)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = "Black-Derman-Toy Model\n"
        s += labelToString("Sigma", self._sigma)
        return s

###############################################################################
