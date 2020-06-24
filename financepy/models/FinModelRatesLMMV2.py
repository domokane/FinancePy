##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from financepy.finutils.FinError import FinError
from financepy.finutils.FinMath import N
from numba import jit, njit, float64, int64, prange
from financepy.finutils.FinSobol import generateSobol

# TO DO: SHITED LOGNORMAL
# TO DO: TERMINAL MEASURE

###############################################################################

def printForwards(fwds):
    
    numPaths = len(fwds)
    numTimes = len(fwds[0])
    numFwds = len(fwds[0][0])

    if numPaths > 10:
        return

    for ip in range(0, numPaths):
        for it in range(0, numTimes):
            for ifwd in range(it, numFwds):
                print(ip, it, ifwd, fwds[ip][it][ifwd])

###############################################################################
                
@njit(float64(int64, int64, float64[:], float64[:], float64[:], float64[:, :]),
      cache=True, fastmath=True, parallel=True)
def LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, rho):
    ''' Implements Rebonato's approximation for the swap rate volatility to be
    used when pricing a swaption that expires in period a for a swap maturing
    at the end of period b taking into account the forward volatility term
    structure (zetas) and the forward-forward correlation matrix rho.. '''

    numPeriods = len(fwd0)
    p = np.zeros(numPeriods)
    p[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, numPeriods):
        p[ix] = p[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    wts = np.zeros(numPeriods)
    pv01ab = 0.0
    for k in range(a+1, b):
        pv01ab += taus[k] * p[k]

    sab = (p[a] - p[b-1])/pv01ab

    for i in range(a, b):
        wts[i] = taus[i] * p[i] / pv01ab

    swaptionVar = 0.0
    for i in range(a, b):
        for j in range(a, b):
            wti = wts[i]
            wtj = wts[j]
            fi = fwd0[i]
            fj = fwd0[j]
            intsigmaij = 0.0

            for k in range(0, a):
                intsigmaij += zetas[i] * zetas[j] * taus[k]

            term = wti * wtj * fi * fj * rho[i][j] * intsigmaij / (sab**2)
            swaptionVar += term

    taua = 0.0
    for i in range(0, a):
        taua += taus[i]

    taub = 0.0
    for i in range(0, b):
        taub += taus[i]

    swaptionVol = np.sqrt(swaptionVar/taua)
    return swaptionVol


###############################################################################


@njit(float64(int64, int64, float64[:], float64[:, :, :], float64[:]),
      cache=True, fastmath=True, parallel=True)
def LMMSimSwaptionVol(a, b, fwd0, fwds, taus):
    ''' Calculates the swap rate volatility using the forwards generated in the
    simulation to see how it compares to Rebonatto estimate. '''

    numPaths = len(fwds)
    numForwards = len(fwds[0])

    if a > numForwards:
        raise FinError("NumPeriods > numForwards")

    if a >= b:
        raise FinError("Swap maturity is before expiry date")

    fwdSwapRateMean = 0.0
    fwdSwapRateVar = 0.0

    for iPath in prange(0, numPaths):

        numeraire = 1.0

        for k in range(0, a):
            numeraire *= (1.0 + taus[k] * fwds[iPath, k, k])

        pv01 = 0.0
        df = 1.0

        for k in range(a, b):
            f = fwds[iPath, a, k]
            tau = taus[k]
            df = df / (1.0 + tau * f)
            pv01 = pv01 + tau * df

        fwdSwapRate = (1.0 - df) / pv01

        fwdSwapRateMean += fwdSwapRate
        fwdSwapRateVar += fwdSwapRate**2

    taua = 0.0
    for i in range(0, a):
        taua += taus[i]

    fwdSwapRateMean /= numPaths
    fwdSwapRateVar = fwdSwapRateVar/numPaths - fwdSwapRateMean**2
    fwdSwapRateVol = np.sqrt(fwdSwapRateVar/taua)
    fwdSwapRateVol /= fwdSwapRateMean
    return fwdSwapRateVol

###############################################################################

@njit(float64[:, :](int64, int64, int64, float64[:, :, :]),
      cache=True, fastmath=True, parallel=True)
def LMMFwdFwdCorrelation(numForwards, numPaths, iTime, fwds):
    ''' Extract forward forward correlation matrix at some future time index
    from the simulated forward rates and return the matrix. '''

    size = numForwards - iTime
    fwdCorr = np.zeros((size, size))

    for iFwd in range(iTime, numForwards):
        for jFwd in range(iFwd, numForwards):

            sumfwdi = 0.0
            sumfwdj = 0.0
            sumfwdifwdi = 0.0
            sumfwdifwdj = 0.0
            sumfwdjfwdj = 0.0

            for p in prange(0, numPaths):
                dfwdi = fwds[p, iTime, iFwd] - fwds[p, iTime-1, iFwd]
                dfwdj = fwds[p, iTime, jFwd] - fwds[p, iTime-1, jFwd]
                sumfwdi += dfwdi
                sumfwdj += dfwdj
                sumfwdifwdi += dfwdi * dfwdi
                sumfwdifwdj += dfwdi * dfwdj
                sumfwdjfwdj += dfwdj * dfwdj

            avgfwdi = sumfwdi / numPaths
            avgfwdj = sumfwdj / numPaths
            avgfwdifwdi = sumfwdifwdi / numPaths
            avgfwdifwdj = sumfwdifwdj / numPaths
            avgfwdjfwdj = sumfwdjfwdj / numPaths

            covii = avgfwdifwdi - avgfwdi * avgfwdi
            covjj = avgfwdjfwdj - avgfwdj * avgfwdj
            covij = avgfwdifwdj - avgfwdi * avgfwdj
            corr = covij / np.sqrt(covii*covjj)

            if abs(covii*covjj) > 1e-20:
                fwdCorr[iFwd-iTime][jFwd-iTime] = corr
                fwdCorr[jFwd-iTime][iFwd-iTime] = corr
            else:
                fwdCorr[iFwd-iTime][jFwd-iTime] = 0.0
                fwdCorr[jFwd-iTime][iFwd-iTime] = 0.0

    return fwdCorr

###############################################################################

@njit(float64[:](float64[:], float64[:], int64, float64, float64[:]),
      cache=True, fastmath=True)
def priceCapsBlack(fwd0, volCaplet, p, K, taus):
    ''' Price a strip of capfloorlets using Black's model using the time grid
    of the LMM model. The prices can be compared with the LMM model prices. '''

    caplet = np.zeros(p+1)
    discFwd = np.zeros(p+1)

    if K <= 0.0:
        raise FinError("Negative strike not allowed.")

    # Set up initial term structure
    discFwd[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for i in range(1, p):
        discFwd[i] = discFwd[i-1] / (1.0 + fwd0[i] * taus[i])

    # Price ATM caplets
    texp = 0.0

    for i in range(1, p):  # 1 to p-1

        K = fwd0[i]
        texp += taus[i]
        vol = volCaplet[i]
        F = fwd0[i]
        d1 = (np.log(F/K) + vol * vol * texp / 2.0) / vol / np.sqrt(texp)
        d2 = d1 - vol * np.sqrt(texp)
        caplet[i] = (F * N(d1) - K * N(d2)) * taus[i] * discFwd[i]

    return caplet

###############################################################################


@njit(float64[:, :](float64[:, :], int64), cache=True, fastmath=True)
def subMatrix(t, N):
    ''' Returns a submatrix of correlation matrix at later time step in the LMM
    simulation which is then used to generate correlated Gaussian RVs. '''

    lent = len(t)
    result = np.zeros((lent-N-1, lent-N-1))

    for i in range(N + 1, lent):
        for j in range(N + 1, lent):
            result[i - N - 1][j - N - 1] = t[i][j]

    return result

###############################################################################


@njit(float64[:, :](float64[:, :]), cache=True, fastmath=True)
def CholeskyNP(rho):
    ''' Numba-compliant wrapper around Numpy cholesky function. '''
    chol = np.linalg.cholesky(rho)
    return chol

###############################################################################


@jit(float64[:, :, :](int64, int64, float64[:], float64[:], float64[:, :],
                      float64[:], int64),
     cache=True, fastmath=True, parallel=True)
def LMMSimulateFwdsNF(numForwards, numPaths, fwd0, zetas, correl, taus, seed):
    ''' Arbitrage-free simulation of forward Libor curves in the spot measure
    given an initial forward curve, volatility term structure and correlation
    structure. The number of forwards at time 0 is given. The 3D matrix of
    forward rates by path, time and forward point is returned. '''

    np.random.seed(seed)

    # Even number of paths for antithetics
    numPaths = 2 * int(numPaths/2)
    halfNumPaths = int(numPaths/2)

    fwd = np.empty((numPaths, numForwards, numForwards))
    fwdB = np.zeros(numForwards)

    discFwd = np.zeros(numForwards)

    # Set up initial term structure
    discFwd[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, numForwards):
        discFwd[ix] = discFwd[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    corr = [None]  # from 0 to p-1
    factors = [None]  # from 0 to p-1

    for ix in range(1, numForwards):  # from 1 to p-1
        matrix = subMatrix(correl, ix - 1)
        corr.append(matrix)
        chol = CholeskyNP(matrix)
        factors.append(chol)

    ###########################################################################
    # I HAVE PROBLEMS AS THE PARALLELISATION CHANGES THE OUTPUT IF RANDS ARE
    # CALCULATED INSIDE THE MAIN LOOP SO I CALCULATE THEM NOW
    ###########################################################################

    gMatrix = np.empty((numPaths, numForwards, numForwards))
    for iPath in range(0, halfNumPaths):
        for j in range(1, numForwards):
            for k in range(0, numForwards-j):
                g = np.random.normal()
                # ANTITHETICS
                gMatrix[iPath, j, k] = g
                gMatrix[iPath + halfNumPaths, j, k] = -g

#    gMatrix = getSobol3(numPaths)

    avgg = 0.0
    stdg = 0.0

    for iPath in prange(0, numPaths):

        # Initial value of forward curve at time 0
        for iFwd in range(0, numForwards):
            fwd[iPath, 0, iFwd] = fwd0[iFwd]

        for j in range(1, numForwards):  # TIME LOOP

            dt = taus[j]
            sqrtdt = np.sqrt(dt)

            for i in range(j, numForwards):  # FORWARDS LOOP

                zi = zetas[i]

                muA = 0.0
                for k in range(j, i+1):
                    rho = corr[j][k-j, i-j]
                    fk = fwd[iPath, j-1, k]
                    zk = zetas[k]
                    tk = taus[k]
                    muA += zi * fk * tk * zk * rho / (1.0 + fk * tk)

                w = 0.0
                for k in range(0, numForwards-j):
                    f = factors[j][i-j, k]
                    w = w + f * gMatrix[iPath, j, k]

                avgg += w
                stdg += w*w

                fwdB[i] = fwd[iPath, j-1, i] \
                    * np.exp(muA * dt - 0.5 * (zi**2) * dt + zi * w * sqrtdt)

                muB = 0.0
                for k in range(j, i+1):
                    rho = corr[j][k-j, i-j]
                    fk = fwdB[k]
                    zk = zetas[k]
                    tk = taus[k]
                    muB += zi * fk * tk * zk * rho / (1.0 + fk * tk)

                muAvg = 0.5*(muA + muB)
                x = np.exp(muAvg * dt - 0.5 * (zi**2) * dt + zi * w * sqrtdt)
                fwd[iPath, j, i] = fwd[iPath, j-1, i] * x

#    avgg = 2.0 * avgg / (numPaths*numForwards*(numForwards-1))
#    stdg = 2.0 * stdg/(numPaths*numForwards*(numForwards-1)) - avgg*avgg          
#    print("AVGG:", avgg, "STDG", stdg)

    return fwd

###############################################################################


@njit(float64[:, :, :](int64, int64, float64[:], float64[:], float64[:],
                       int64), cache=True, fastmath=True, parallel=True)
def LMMSimulateFwds1F(numForwards, numPaths, fwd0, zetas, taus, seed):
    ''' One factor Arbitrage-free simulation of forward Libor curves in the
    spot measure following Hull Page 768. Given an initial forward curve,
    volatility term structure. The 3D matrix of forward rates by path, time
    and forward point is returned. '''

    np.random.seed(seed)
    fwdB = np.zeros(numForwards)

    # Even number of paths for antithetics
    numPaths = 2 * int(numPaths/2)
    halfNumPaths = int(numPaths/2)

    fwd = np.empty((numPaths, numForwards, numForwards))
    fwdB = np.zeros(numForwards)

    discFwd = np.zeros(numForwards)

    # Set up initial term structure
    discFwd[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, numForwards):
        discFwd[ix] = discFwd[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    numTimes = numForwards
    gMatrix = np.empty((numPaths, numTimes))
    for iPath in range(0, halfNumPaths):
        for j in range(0, numTimes):
            g = np.random.normal()
            gMatrix[iPath, j] = g
            gMatrix[iPath + halfNumPaths, j] = -g

    for iPath in prange(0, numPaths):
        # Initial value of forward curve at time 0
        for iFwd in range(0, numForwards):
            fwd[iPath, 0, iFwd] = fwd0[iFwd]

        for j in range(0, numForwards-1):  # TIME LOOP
            dtj = taus[j]
            sqrtdtj = np.sqrt(dtj)
            w = gMatrix[iPath, j]

            for k in range(j, numForwards):  # FORWARDS LOOP
                zkj = zetas[k-j-1]
                mu = 0.0
                for i in range(j+1, k+1):
                    fi = fwd[iPath, j, i]
                    zij = zetas[i-j-1]
                    ti = taus[i]
                    mu += zkj * fi * ti * zij / (1.0 + fi * ti)

#                    print(iPath, j, k, mu, fi, zij, ti, dt, zkj, w)
                fwdB[k] = fwd[iPath, j, k] \
                   * np.exp(mu * dtj - 0.5*(zkj**2) * dtj + zkj * w * sqrtdtj)

                muB = 0.0 
                for i in range(j+1, k+1):
                    fi = fwdB[k]
                    zij = zetas[i-j-1]
                    ti = taus[i]
                    mu += zkj * fi * ti * zij / (1.0 + fi * ti)

                muC = 0.5*(mu+muB)

                x = np.exp(muC * dtj - 0.5 * (zkj**2) * dtj + zkj * w * sqrtdtj)
                fwd[iPath, j+1, k] = fwd[iPath, j, k] * x

    return fwd

###############################################################################


@njit(float64[:](int64, int64, float64, float64[:], float64[:, :, :],
                 float64[:], int64),
      cache=True, fastmath=True, parallel=True)
def LMMCapFlrPricer(numPeriods, numPaths, K, fwd0, fwds, taus, isCap):
    ''' Function to price a strip of cap or floorlets in accordance with the
    simulated forward curve dynamics. '''

    maxPaths = len(fwds)
    maxForwards = len(fwds[0])

    if numPeriods > maxForwards:
        raise FinError("NumPeriods > numForwards")

    if numPaths > maxPaths:
        raise FinError("NumPaths > MaxPaths")

    discFactor = np.zeros(maxForwards)
    capFlrLets = np.zeros(maxForwards)
    capFlrLetValues = np.zeros(maxForwards)
    numeraire = np.zeros(maxForwards)

    # Set up initial term structure
    discFactor[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, maxForwards):
        discFactor[ix] = discFactor[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for iPath in range(0, numPaths):

        periodRoll = 1.0
        libor = fwds[iPath, 0, 0]
        capFlrLets[0] = max(K - libor, 0) * taus[0]

        for j in range(1, numPeriods):  # TIME LOOP

            libor = fwds[iPath, j, j]
            if j == 1:
                if isCap == 0:
                    capFlrLets[j] = max(K - libor, 0) * taus[j]
                else:
                    capFlrLets[j] = max(libor - K, 0) * taus[j]

                numeraire[0] = 1.0 / discFactor[0]
            else:
                if isCap == 1:
                    capFlrLets[j] = max(libor - K, 0) * taus[j]
                elif isCap == 0:
                    capFlrLets[j] = max(K - libor, 0) * taus[j]
                else:
                    raise FinError("isCap should be 0 or 1")

            periodRoll = (1.0 + libor * taus[j])
            numeraire[j] = numeraire[j - 1] * periodRoll

        for iFwd in range(0, numPeriods):
            capFlrLetValues[iFwd] += capFlrLets[iFwd] / numeraire[iFwd]

    for iFwd in range(0, numPeriods):
        capFlrLetValues[iFwd] *= 100.0 / numPaths

    return capFlrLetValues

###############################################################################

@njit(float64(float64, int64, int64, float64[:], float64[:, :, :],
              float64[:]), cache=True, fastmath=True, parallel=True)
def LMMSwapPricer(cpn, numPeriods, numPaths, fwd0, fwds, taus):
    ''' Function to reprice a basic swap using the simulated forward Libors.
    '''

    maxPaths = len(fwds)
    maxForwards = len(fwds[0])

    if numPeriods > maxForwards:
        raise FinError("NumPeriods > numForwards")

    if numPaths > maxPaths:
        raise FinError("NumPaths > MaxPaths")

    discFactor = np.zeros(maxForwards)
    numeraire = np.zeros(maxForwards)
    sumFixed = 0.0
    sumFloat = 0.0
    fixedFlows = np.zeros(maxForwards)
    floatFlows = np.zeros(maxForwards)

    # Set up initial term structure
    discFactor[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, maxForwards):
        discFactor[ix] = discFactor[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for iPath in range(0, numPaths):

        periodRoll = 1.0
        libor = fwds[iPath, 0, 0]
        floatFlows[0] = libor * taus[0]
        fixedFlows[0] = cpn * taus[0]
        numeraire[0] = 1.0 / discFactor[0]

        for j in range(1, numPeriods):  # TIME LOOP

            libor = fwds[iPath, j, j]

            if j == 1:
                fixedFlows[j] = cpn * taus[j]
                floatFlows[j] = libor * taus[j]
            else:
                fixedFlows[j] = fixedFlows[j-1] * periodRoll + cpn * taus[j]
                floatFlows[j] = floatFlows[j-1] * periodRoll + libor * taus[j]

            periodRoll = (1.0 + libor * taus[j])
            numeraire[j] = numeraire[j - 1] * periodRoll

        for iFwd in range(0, numPeriods):
            sumFloat += floatFlows[iFwd] / numeraire[iFwd]
            sumFixed += fixedFlows[iFwd] / numeraire[iFwd]

    sumFloat /= numPaths
    sumFixed /= numPaths
    v = sumFixed - sumFloat
    pv01 = sumFixed/cpn
    swapRate = sumFloat/pv01

    print("FLOAT LEG:", sumFloat)
    print("FIXED LEG:", sumFixed)
    print("SWAP RATE:", swapRate)
    print("NET VALUE:", v)
    return v

###############################################################################
@njit(float64(float64, int64, int64, int64, float64[:], float64[:, :, :],
              float64[:], int64), cache=True, fastmath=True, parallel=False)
def LMMSwaptionPricer(strike, a, b, numPaths, fwd0, fwds, taus, isPayer):
    ''' Function to price a European swaption using the simulated forward
    curves. '''

    maxPaths = len(fwds)
    maxForwards = len(fwds[0])

    if a > maxForwards:
        raise FinError("NumPeriods > numForwards")

    if a >= b:
        raise FinError("Swap maturity is before expiry date")

    if numPaths > maxPaths:
        raise FinError("NumPaths > MaxPaths")

    discFactor = np.zeros(maxForwards)
    pv01 = np.zeros(maxForwards)

    # Set up initial term structure
    discFactor[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, b):
        discFactor[ix] = discFactor[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    sumPayRecSwaption = 0.0

    for iPath in range(0, numPaths):

        numeraire = 1.0
        for k in range(0, a):
            numeraire *= (1.0 + taus[k] * fwds[iPath, k, k])

        pv01 = 0.0
        df = 1.0

        # Value the swap as if we were at time a with forward curve known
        for k in range(a, b):
            f = fwds[iPath, a, k]
            tau = taus[k]
            df = df / (1.0 + tau * f)
            pv01 = pv01 + tau * df

        fwdSwapRate = (1.0 - df) / pv01

        if isPayer == 1:
            payRecSwaption = max(fwdSwapRate - strike, 0.0) * pv01
        elif isPayer == 0:
            payRecSwaption = max(strike - fwdSwapRate, 0.0) * pv01
        else:
            raise FinError("Unknown payRecSwaption value - must be 0 or 1")

        sumPayRecSwaption += payRecSwaption / numeraire

    payRecPrice = sumPayRecSwaption / numPaths
    return payRecPrice

###############################################################################


@njit(float64[:](float64, int64, int64, float64[:], float64[:, :, :],
              float64[:]), cache=True, fastmath=True, parallel=True)
def LMMRatchetCapletPricer(spread, numPeriods, numPaths, fwd0, fwds, taus):
    ''' Price a ratchet using the simulated Libor rates.'''

    maxPaths = len(fwds)
    maxTimes = len(fwds[0])
    maxForwards = len(fwds[0][0])

    if numPeriods > maxForwards:
        raise FinError("NumPeriods > numForwards")

    if numPaths > maxPaths:
        raise FinError("NumPaths > MaxPaths")

    discFactor = np.zeros(maxForwards)
    numeraire = np.zeros(maxForwards)
    ratchetCaplets = np.zeros(maxForwards)
    ratchetCapletValues = np.zeros(maxForwards)

    # Set up initial term structure
    discFactor[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, maxForwards):
        discFactor[ix] = discFactor[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for iPath in range(0, numPaths):

        periodRoll = 1.0
        libor = fwds[iPath, 0, 0]
        ratchetCaplets[0] = 0.0

        for j in range(1, numPeriods):  # TIME LOOP

            prevLibor = libor
            K = prevLibor + spread
            libor = fwds[iPath, j, j]

            if j == 1:
                ratchetCaplets[j] = max(libor - K, 0.0) * taus[j]
                numeraire[0] = 1.0 / discFactor[0]
            else:
                ratchetCaplets[j] = max(libor - K, 0.0) * taus[j]

            periodRoll = (1.0 + libor * taus[j])
            numeraire[j] = numeraire[j - 1] * periodRoll

#            print("iPath:", iPath, "j:", j, "prevL:", prevLibor, "K:", K, "L:", libor, "RC[i]", ratchetCap[j])

        for iFwd in range(0, numPeriods):
            ratchetCapletValues[iFwd] += ratchetCaplets[iFwd] / numeraire[iFwd]

    notional = 100.0
    for iFwd in range(0, numPeriods):
        ratchetCapletValues[iFwd] *= notional / numPaths

    return ratchetCapletValues

###############################################################################

@njit(float64[:](float64, int64, int64, float64[:], float64[:, :, :],
              float64[:]), cache=True, fastmath=True, parallel=True)
def LMMStickyCapletPricer(spread, numPeriods, numPaths, fwd0, fwds, taus):
    ''' Price a sticky cap using the simulated Libor rates. '''

    maxPaths = len(fwds)
    maxTimes = len(fwds[0])
    maxForwards = len(fwds[0][0])

    if numPeriods > maxForwards:
        raise FinError("NumPeriods > numForwards")

    if numPaths > maxPaths:
        raise FinError("NumPaths > MaxPaths")

    discFactor = np.zeros(maxForwards)
    numeraire = np.zeros(maxForwards)
    stickyCaplets = np.zeros(maxForwards)
    stickyCapletValues = np.zeros(maxForwards)

    # Set up initial term structure
    discFactor[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, maxForwards):
        discFactor[ix] = discFactor[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for iPath in range(0, numPaths):

        periodRoll = 1.0
        libor = fwds[iPath, 0, 0]
        stickyCaplets[0] = 0.0
        K = libor

        for j in range(1, numPeriods):  # TIME LOOP

            prevLibor = libor
            K = min(prevLibor, K) + spread
            libor = fwds[iPath, j, j]

            if j == 1:
                stickyCaplets[j] = max(libor-K, 0.0) * taus[j]
                numeraire[0] = 1.0 / discFactor[0]
            else:
                stickyCaplets[j] = max(libor - K, 0.0) * taus[j]

            periodRoll = (1.0 + libor * taus[j])
            numeraire[j] = numeraire[j - 1] * periodRoll

        for iFwd in range(0, numPeriods):
            stickyCapletValues[iFwd] += stickyCaplets[iFwd] / numeraire[iFwd]

    notional = 100.0
    for iFwd in range(0, numPeriods):
        stickyCapletValues[iFwd] *= notional / numPaths

    return stickyCapletValues

###############################################################################

