##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
from numba import njit, jit
from math import ceil

from ..finutils.FinError import FinError
from ..finutils.FinMath import N, accruedInterpolator
from ..market.curves.FinInterpolate import FinInterpTypes, uinterpolate
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinOptionTypes import FinOptionExerciseTypes

interp = FinInterpTypes.FLAT_FORWARDS.value

small = 1e-10

###############################################################################
###############################################################################
# dr = (theta(t) - r) dt + sigma * dW
###############################################################################
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


@njit(fastmath=True, cache=True)
def P_Fast(t, T, Rt, delta, pt, ptd, pT, _sigma, _a):
    ''' Forward discount factor as seen at some time t which may be in the
    future for payment at time T where Rt is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor. '''

    if abs(_a) < small:
        _a = small

    BtT = (1.0 - np.exp(-_a*(T-t)))/_a
    BtDelta = (1.0 - np.exp(-_a * delta))/_a

    term1 = np.log(pT/pt) - (BtT/BtDelta) * np.log(ptd/pt)

    term2 = (_sigma**2)*(1.0-np.exp(-2.0*_a*t)) \
        * BtT * (BtT - BtDelta)/(4.0*_a)

    logAhat = term1 - term2
    BhattT = (BtT/BtDelta) * delta
    p = np.exp(logAhat - BhattT * Rt)
    return p

###############################################################################


@njit(fastmath=True, cache=True)
def buildTree_Fast(a, sigma, treeTimes, numTimeSteps, discountFactors):
    ''' Fast tree construction using Numba. '''
    treeMaturity = treeTimes[-1]
    dt = treeMaturity / (numTimeSteps+1)
    dR = sigma * np.sqrt(3.0 * dt)
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
            sumQZ += Q[m, j+N] * np.exp(-rdt)
        alpha[m] = np.log(sumQZ/discountFactors[m+1]) / dt

        for j in range(-nm, nm+1):
            jN = j + N
            rt[m, jN] = alpha[m] + j*dR

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):
            jN = j + N
            rdt = rt[m, jN] * dt
            z = np.exp(-rdt)

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

###############################################################################


@njit(fastmath=True, cache=True)
def americanBondOption_Tree_Fast(texp,
                                 strikePrice,
                                 face,
                                 couponTimes, couponAmounts,
                                 americanExercise,
                                 _sigma,
                                 _a,
                                 _Q,
                                 _pu, _pm, _pd,
                                 _rt,
                                 _dt,
                                 _treeTimes,
                                 _dfTimes,
                                 _dfValues):
    ''' Value an option on a bond with coupons that can have European or
    American exercise. Some minor issues to do with handling coupons on
    the option expiry date need to be solved. Also this function should be
    moved out of the class so it can be sped up using NUMBA. '''

#    DEBUG = False

    numTimeSteps, numNodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * dt))
    expiryStep = int(texp/dt + 0.50)
    treeFlows = np.zeros(numTimeSteps)

#    if DEBUG:
#        filename = "bondOptHW_" + str(numTimeSteps) + ".txt"
#        outfile = open(filename, "w")

    # Want to add coupons before expiry to the grid so that we can value
    # their impact on the decision to exercise the option early
    numCoupons = len(couponTimes)

#    if DEBUG:
#        outfile.write("TEXP: %9.5f\n" % texp)

    # Flows that fall on the expiry date are not included
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        if tcpn <= texp:
            n = int(tcpn/dt)
            ttree = _treeTimes[n]
            df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
            df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
            treeFlows[n] += couponAmounts[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mappedTimes = np.array([0.0])
    mappedAmounts = np.array([0.0])
    for n in range(0, len(_treeTimes)):
        if treeFlows[n] > 0.0:
            np.append(mappedTimes, _treeTimes[n])
            np.append(mappedAmounts, treeFlows[n])

    # Need future cashflows which are exact time and size for accrued at texp
    for n in range(0, numCoupons):
        if couponTimes[n] > texp:
            np.append(mappedTimes, couponTimes[n])
            np.append(mappedAmounts, couponAmounts[n])

    ###########################################################################

    accrued = np.zeros(numTimeSteps)
    for m in range(0, expiryStep+1):
        ttree = _treeTimes[m]
        accrued[m] = accruedInterpolator(ttree, mappedTimes, mappedAmounts)
        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face

    ###########################################################################

#    if DEBUG:
#        outfile.write("COUPON TIMES FLOWS \n")
#        for i in range(0, len(couponTimes)):
#            t = couponTimes[i]
#            f = couponAmounts[i]
#            outfile.write("%9.5f  %9.5f\n" % (t, f*face))

        # outfile.write("MAPPED TIMES FLOWS \n")
        # for i in range(0, len(mappedTimes)):
        #     t = mappedTimes[i]
        #     f = mappedAmounts[i]
        #     outfile.write("%9.5f  %9.5f\n" % (t, f*face))

        # outfile.write("TREE TIMES, FLOWS AND ACCRUED \n")
        # for i in range(0, len(accrued)):
        #     t = _treeTimes[i]
        #     f = treeFlows[i]
        #     a = accrued[i]
        #     outfile.write("%9.5f  %9.5f  %9.5f\n" % (t, f*face, a))

    ###########################################################################

    callOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    putOptionValues = np.zeros(shape=(numTimeSteps, numNodes))
    bondValues = np.zeros(shape=(numTimeSteps, numNodes))

    ptexp = uinterpolate(texp, _dfTimes, _dfValues, interp)
    ptdelta = uinterpolate(texp+dt, _dfTimes, _dfValues, interp)

    flow = treeFlows[expiryStep]

    # if DEBUG:
    #     outfile.write("FUTURE FLOWS AT EXPIRY\n")

    cpn = 0.0
    zcb = 0.0

    ###########################################################################
    # As the HW model has a closed form solution for the bond price, I use
    # this fact to calculate the bond price at expiry on the tree nodes
    ###########################################################################

    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + jmax
        rt = _rt[expiryStep, kN]
        bondPrice = 0.0
        for i in range(0, numCoupons):
            tflow = couponTimes[i]
            if tflow > texp:
                ptflow = uinterpolate(tflow, _dfTimes, _dfValues, interp)
                zcb = P_Fast(texp, tflow, rt, dt, ptexp, ptdelta, ptflow,
                             _sigma, _a)
                cpn = couponAmounts[i]
                bondPrice += cpn * face * zcb

        bondPrice += zcb * face
        bondValues[expiryStep, kN] = bondPrice + flow * face

        # if DEBUG:
        #     outfile.write("NODE: %d BOND VALUE EXPIRY: %9.5f\n" %
        #                   (kN, bondValues[expiryStep, kN]))

    # Now consider exercise of the option on the expiry date
    # Start with the value of the bond at maturity and overwrite values
    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + jmax
        cleanPrice = bondValues[expiryStep, kN] - accrued[expiryStep]
        callOptionValues[expiryStep, kN] = max(cleanPrice - strikePrice, 0.0)
        putOptionValues[expiryStep, kN] = max(strikePrice - cleanPrice, 0.0)

    # if DEBUG:
    #     outfile.write("OPTION    m    k     treeT      treeF  bondVal   cleanPx     Accd    CallVal    PutVal\n")

    # Now step back to today considering early exercise
    for m in range(expiryStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face

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

                callExercise = max(cleanPrice - strikePrice, 0.0)
                putExercise = max(strikePrice - cleanPrice, 0.0)

                hold = callOptionValues[m, kN]
                callOptionValues[m, kN] = max(callExercise, hold)

                hold = putOptionValues[m, kN]
                putOptionValues[m, kN] = max(putExercise, hold)

    #     if DEBUG:
    #         outfile.write("STEP: %5d %5d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n" % \
    #                       (m,
    #                        k,
    #                        _treeTimes[m],
    #                        treeFlows[m],
    #                        bondValues[m, kN],
    #                        cleanPrice,
    #                        accrued[m],
    #                        callOptionValues[m, kN],
    #                        putOptionValues[m, kN]))

    # if DEBUG:
    #     outfile.close()

    return callOptionValues[0, jmax], putOptionValues[0, jmax]

###############################################################################


@jit(fastmath=True, cache=True)
def bermudanSwaption_Tree_Fast(texp, tmat, strikePrice, face,
                               couponTimes, couponFlows,
                               exerciseTypeInt,
                               _dfTimes, _dfValues,
                               _treeTimes, _Q, _pu, _pm, _pd, _rt, _dt, _a):
    ''' Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time. '''

    # DEBUG = False # True

    ###########################################################################

    numTimeSteps, numNodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)
    numCoupons = len(couponTimes)

    # if DEBUG:
    #     filename = "swaptionHW_" + str(numTimeSteps) + ".txt"
    #     outfile = open(filename, "w")

    # Tree flows go all the way out to the bond maturity date
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(round(tcpn/_dt, 0))
        ttree = _treeTimes[n]
        df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mappedTimes = [0.0]
    mappedAmounts = [0.0]
    for n in range(1, len(_treeTimes)):

        accdAtExpiry = 0.0
        if _treeTimes[n-1] < texp and _treeTimes[n] >= texp:
            mappedTimes.append(texp)
            mappedAmounts.append(accdAtExpiry)

        if treeFlows[n] > 0.0:
            mappedTimes.append(_treeTimes[n])
            mappedAmounts.append(treeFlows[n])

    ###########################################################################

    accrued = np.zeros(numTimeSteps)
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

    ###########################################################################

    payValues = np.zeros(shape=(numTimeSteps, numNodes))
    recValues = np.zeros(shape=(numTimeSteps, numNodes))
    swapValues = np.zeros(shape=(numTimeSteps, numNodes))

    # Start with the value of the bond at maturity
    for k in range(0, numNodes):
        swapValues[maturityStep, k] = (1.0 + treeFlows[maturityStep]) * face

    # if DEBUG:
    #     outfile.write("m: %d k: %d EXPIRY: %9.5f\n" %
    #                   (maturityStep, k, swapValues[maturityStep, k]))

    N = jmax

    # Now step back to today considering early exercise
    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face

        for k in range(-nm, nm+1):
            kN = k + N
            rt = _rt[m, kN]
            df = np.exp(-rt * _dt)
            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = swapValues[m+1, kN]
                vm = swapValues[m+1, kN-1]
                vd = swapValues[m+1, kN-2]
                v = (pu*vu + pm*vm + pd*vd) * df
                swapValues[m, kN] = v
            elif k == jmax:
                vu = swapValues[m+1, kN+2]
                vm = swapValues[m+1, kN+1]
                vd = swapValues[m+1, kN]
                v = (pu*vu + pm*vm + pd*vd) * df
                swapValues[m, kN] = v
            else:
                vu = swapValues[m+1, kN+1]
                vm = swapValues[m+1, kN]
                vd = swapValues[m+1, kN-1]
                v = (pu*vu + pm*vm + pd*vd) * df
                swapValues[m, kN] = v

            swapValues[m, kN] += flow
            vpay = 0.0
            vrec = 0.0

            if k == jmax:
                vu = payValues[m+1, kN]
                vm = payValues[m+1, kN-1]
                vd = payValues[m+1, kN-2]
                vpay = (pu*vu + pm*vm + pd*vd) * df
            elif k == jmax:
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
            elif k == jmax:
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

            cleanPrice = swapValues[m, kN] - accrued[m]
            payExercise = max(strikePrice - cleanPrice, 0.0)
            recExercise = max(cleanPrice - strikePrice, 0.0)

            holdPay = payValues[m, kN]
            holdRec = recValues[m, kN]

            if exerciseTypeInt == 1 and m == expiryStep:
                # European style on expiry date

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

            elif exerciseTypeInt == 2 and m == expiryStep:

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

            elif exerciseTypeInt == 2 and flow > 0.0 and m >= expiryStep:
                # Bermudan style on coupon dates after expiry date

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

            elif exerciseTypeInt == 3 and m >= expiryStep:
                # American style on all dates after expiry date

                payValues[m, kN] = max(payExercise, holdPay)
                recValues[m, kN] = max(recExercise, holdRec)

    #     if DEBUG:
    #         outfile.write("STEP: %5d %5d %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n" % \
    #                       (m,
    #                        k,
    #                        _treeTimes[m],
    #                        treeFlows[m],
    #                        swapValues[m, kN],
    #                        cleanPrice,
    #                        accrued[m],
    #                        payValues[m, kN],
    #                        recValues[m, kN]))

    # if DEBUG:
    #     outfile.close()

    return payValues[0, jmax], recValues[0, jmax]

###############################################################################
# TODO: CHECK ACCRUED AND COUPONS TO SEE IF IT WORKS FOR LOW TREE STEPS
###############################################################################


@njit(fastmath=True, cache=True)
def callablePuttableBond_Tree_Fast(couponTimes, couponAmounts,
                                   callTimes, callPrices,
                                   putTimes, putPrices, face,
                                   _sigma, _a, _Q,  # IS SIGMA USED ?
                                   _pu, _pm, _pd, _rt, _dt, _treeTimes,
                                   _dfTimes, _dfValues):
    ''' Value an option on a bond with coupons that can have European or
    American exercise. Some minor issues to do with handling coupons on
    the option expiry date need to be solved. Also this function should be
    moved out of the class so it can be sped up using NUMBA. '''

    if np.any(couponTimes < 0.0):
        raise FinError("No coupon times can be before the value date.")

    numTimeSteps, numNodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * dt))
    tmat = couponTimes[-1]
    maturityStep = int(tmat/dt + 0.50)

    ###########################################################################
    # Map coupons onto tree while preserving their present value
    ###########################################################################

    treeFlows = np.zeros(numTimeSteps)

    numCoupons = len(couponTimes)
    for i in range(0, numCoupons):
        tcpn = couponTimes[i]
        n = int(round(tcpn/dt, 0))
        ttree = _treeTimes[n]
        df_flow = uinterpolate(tcpn, _dfTimes, _dfValues, interp)
        df_tree = uinterpolate(ttree, _dfTimes, _dfValues, interp)
        treeFlows[n] += couponAmounts[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mappedTimes = [0.0]
    mappedAmounts = [0.0]
    for n in range(1, len(_treeTimes)):
        if treeFlows[n] > 0.0:
            mappedTimes.append(_treeTimes[n])
            mappedAmounts.append(treeFlows[n])

    ###########################################################################

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
    nm = min(maturityStep, jmax)
    vcall = treeCallValue[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * face
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm+1):
        kN = k + jmax
        bondValues[m, kN] = (1.0 + treeFlows[m]) * face
        callPutBondValues[m, kN] = value

    # Now step back to today considering early put and call
    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face
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
            elif k == jmax:
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
            elif k == jmax:
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

        if tcpn >= texp:  # CHECK IF IT SHOULD BE >=
            ptcpn = uinterpolate(tcpn, dfTimes, dfValues, interp)
            zcb = P_Fast(texp, tcpn, rt, dt, ptexp, ptdelta, ptcpn,
                         self._sigma, self._a)
            pv = pv + zcb * cpn
#            print("TCPN", tcpn, "ZCB", zcb, "CPN", cpn, "PV", pv)

    if tcpn >= texp:
        pv = pv + zcb

#    print("TCPN", tcpn, "ZCB", zcb, "PRI", 1.0, "PV", pv)

    accd = accruedInterpolator(texp, cpnTimes, cpnAmounts)

#    print("texp:", texp)
#    print("cpnTimes:", cpnTimes)
#    print("cpnAmounts:", cpnAmounts)

    pv_clean = pv - accd
    obj = face * pv_clean - strikePrice

#    print("FWD PRICE", rt, pv, accd, strikePrice, obj)
    return obj

###############################################################################


class FinModelRatesHW():

    def __init__(self,
                 sigma,
                 a,
                 numTimeSteps=100,
                 useJamshidian=True):
        ''' Constructs the Hull-White rate model. The speed of mean reversion
        a and volatility are passed in. The short rate process is given by
        dr = (theta(t) - ar) * dt  + sigma * dW. The model will switch to use
        Jamshidian's approach where possible unless the useJamshidian flag is
        set to false in which case it uses the trinomial Tree. '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        self._sigma = sigma
        self._a = a
        self._numTimeSteps = numTimeSteps
        self._useJamshidian = useJamshidian

        self._Q = None
        self._r = None
        self._treeTimes = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discountCurve = None
        self._treeBuilt = False

###############################################################################

    def optionOnZCB(self,
                    texp,
                    tmat,
                    strike,
                    face,
                    dfTimes,
                    dfValues):
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

        if abs(a) < small:
            a = small

        sigmap = (sigma/a) * (1.0 - np.exp(-a*(tmat-texp)))
        sigmap *= np.sqrt((1.0-np.exp(-2.0*a*texp))/2.0/a)

        if abs(sigmap) < small:
            sigmap = small

        h = np.log((face*ptmat)/(strike * ptexp)) / sigmap + sigmap / 2.0
        callValue = face * ptmat * N(h) - strike * ptexp * N(h-sigmap)
        putValue = strike * ptexp * N(-h+sigmap) - face * ptmat * N(-h)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def europeanBondOption_Jamshidian(self,
                                      texp,
                                      strike,
                                      face,
                                      cpnTimes,
                                      cpnAmounts,
                                      dfTimes,
                                      dfValues):
        ''' Valuation of a European bond option using the Jamshidian
        deconstruction of the bond into a strip of zero coupon bonds with the
        short rate that would make the bond option be at the money forward. '''

        argtuple = (self, texp, cpnTimes, cpnAmounts,
                    dfTimes, dfValues, strike, face)

        x0 = 0.05
        rstar = optimize.newton(fwdFullBondPrice, x0=x0, fprime=None,
                                args=argtuple, tol=1e-10, maxiter=50,
                                fprime2=None)

        # Now we price a series of zero coupon bonds using this short rate
        dt = 0.0001
        numCoupons = len(cpnTimes)

        ptexp = uinterpolate(texp, dfTimes, dfValues, interp)
        ptdelta = uinterpolate(texp+dt, dfTimes, dfValues, interp)

        callValue = 0.0
        putValue = 0.0

        for i in range(0, numCoupons):
            tcpn = cpnTimes[i]
            cpn = cpnAmounts[i]

            if tcpn > texp:  # coupons on the expiry date are ignored

                ptcpn = uinterpolate(tcpn, dfTimes, dfValues, interp)

                strike = P_Fast(texp, tcpn, rstar, dt, ptexp, ptdelta,
                                ptcpn, self._sigma, self._a)

                v = self.optionOnZCB(texp, tcpn, strike, 1.0,
                                     dfTimes, dfValues)

                call = v['call']
                put = v['put']

#                print(tcpn, ptcpn, strike, call, put)

                callValue += call * cpn * face
                putValue += put * cpn * face

        callValue += call * face
        putValue += put * face

#        print("CALL VAL", callValue, "PUT VAL", putValue)
        return {'call': callValue, 'put': putValue}

###############################################################################

    def europeanBondOption_Tree(self,
                                texp,
                                strikePrice,
                                face,
                                cpnTimes,
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

        callValue = 0.0
        putValue = 0.0
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
            putValue += q * putPayoff
            callValue += q * callPayoff

        return {'call': callValue, 'put': putValue}

###############################################################################

    def optionOnZeroCouponBond_Tree(self,
                                    texp,
                                    tmat,
                                    strikePrice,
                                    face):
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

        callValue = 0.0
        putValue = 0.0

        for k in range(0, numNodes):

            q = self._Q[expiryStep, k]
            rt = self._rt[expiryStep, k]

            zcb = P_Fast(texp, tmat, rt, dt, ptexp, ptdelta, ptmat,
                         self._sigma, self._a)

            putPayoff = max(strikePrice - zcb * face, 0.0)
            callPayoff = max(zcb * face - strikePrice, 0.0)
            putValue += q * putPayoff
            callValue += q * callPayoff

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudanSwaption(self, texp, tmat, strike, face,
                         couponTimes, couponFlows, exerciseType):
        ''' Swaption that can be exercised on specific dates over the exercise
        period. Due to non-analytical bond price we need to extend tree out to
        bond maturity and take into account cash flows through time. '''

        exerciseTypeInt = optionExerciseTypesToInt(exerciseType)

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
                                         self._pu, self._pm, self._pd,
                                         self._rt,
                                         self._dt, self._a)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def americanBondOption_Tree(self,
                                texp,
                                strikePrice,
                                face,
                                couponTimes,
                                couponAmounts,
                                americanExercise):
        ''' Value an option on a bond with coupons that can have European or
        American exercise. Some minor issues to do with handling coupons on
        the option expiry date need to be solved. '''

        callValue, putValue \
            = americanBondOption_Tree_Fast(texp, strikePrice, face,
                                           couponTimes, couponAmounts,
                                           americanExercise,
                                           self._sigma, self._a,
                                           self._Q,
                                           self._pu, self._pm, self._pd,
                                           self._rt, self._dt,
                                           self._treeTimes,
                                           self._dfTimes, self._dfValues)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def callablePuttableBond_Tree(self,
                                  couponTimes,
                                  couponAmounts,
                                  callTimes,
                                  callPrices,
                                  putTimes,
                                  putPrices,
                                  face):
        ''' Value an option on a bond with coupons that can have European or
        American exercise. Some minor issues to do with handling coupons on
        the option expiry date need to be solved. Also this function should be
        moved out of the class so it can be sped up using NUMBA. '''

        callTimes = np.array(callTimes)
        putTimes = np.array(putTimes)

        callPrices = np.array(callPrices)
        putPrices = np.array(putPrices)

        v = callablePuttableBond_Tree_Fast(couponTimes, couponAmounts,
                                           callTimes, callPrices,
                                           putTimes, putPrices, face,
                                           self._sigma, self._a,
                                           self._Q,
                                           self._pu, self._pm, self._pd,
                                           self._rt, self._dt,
                                           self._treeTimes,
                                           self._dfTimes, self._dfValues)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

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
        zeroRate = -np.log(p)/tmat
        return p, zeroRate

###############################################################################

    def buildTree(self, treeMat, dfTimes, dfValues):
        ''' Build the trinomial tree. '''

        if isinstance(dfTimes, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(dfValues, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        # I wish to add on an additional time to the tree so that the second
        # last time corresponds to a maturity treeMat. For this reason I scale
        # up the maturity date of the tree as follows
        treeMaturity = treeMat * (self._numTimeSteps+1)/self._numTimeSteps

        # The vector of times goes out to this maturity
        treeTimes = np.linspace(0.0, treeMaturity, self._numTimeSteps + 2)
        self._treeTimes = treeTimes

        dfTree = np.zeros(shape=(self._numTimeSteps+2))
        dfTree[0] = 1.0

        for i in range(1, self._numTimeSteps+2):
            t = treeTimes[i]
            dfTree[i] = uinterpolate(t, dfTimes, dfValues, interp)

        self._dfTimes = dfTimes
        self._dfValues = dfValues

        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = buildTree_Fast(self._a, self._sigma,
                             treeTimes, self._numTimeSteps, dfTree)

        return

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''
        s = "Hull-White Model\n"
        s += labelToString("Sigma", self._sigma)
        s += labelToString("a", self._a)
        return s

###############################################################################
