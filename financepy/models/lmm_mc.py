##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
import numpy as np
from numba import njit, float64, int64  # , prange DOES NOT WORK ON GITHUB

from ..utils.error import FinError
from ..utils.math import N
from ..utils.math import norminvcdf
from ..models.sobol import get_uniform_sobol

# TO DO: SHIFTED LOGNORMAL
# TO DO: TERMINAL MEASURE
# TO DO:: CALIBRATION

USE_PARALLEL = False

###############################################################################

""" This module manages the Ibor Market Model and so stores a specific MC
    forward rate simulation of a 3D matrix of num_paths x num_fwds
    x (num_fwds-1)/2 elements. This is a lognormal model although a shifted
    Lognormal rate is also allowed. Implementations include 1 factor, M factor
    where the volatility curve per factor is provided and a full N-factor corr-
    elation matrix where a Cholesky is done to decompose the N factors. """

###############################################################################


class ModelLMMModelTypes(Enum):
    LMM_ONE_FACTOR = 1
    LMM_HW_M_FACTOR = 2
    LMM_FULL_N_FACTOR = 3

###############################################################################


def lmm_print_forwards(fwds):
    """ Helper function to display the simulated Ibor rates. """

    num_paths = len(fwds)
    num_times = len(fwds[0])
    num_fwds = len(fwds[0][0])

    if num_paths > 10:
        return

    for ip in range(0, num_paths):
        for it in range(0, num_times):

            print("Path: %3d Time: %3d" % (ip, it), end=""),

            for ifwd in range(0, it):
                print("%8s" % ("-"), end=""),

            for ifwd in range(it, num_fwds):
                print("%8.4f" % (fwds[ip][it][ifwd]*100.0), end=""),

            print("")


###############################################################################


@njit(float64(int64, int64, float64[:], float64[:], float64[:], float64[:, :]),
      cache=True, fastmath=True)
def lmm_swaption_vol_approx(a, b, fwd0, taus, zetas, rho):
    """ Implements Rebonato's approximation for the swap rate volatility to be
    used when pricing a swaption that expires in period a for a swap maturing
    at the end of period b taking into account the forward volatility term
    structure (zetas) and the forward-forward correlation matrix rho.. """

    num_periods = len(fwd0)

#    if len(taus) != num_periods:
#        raise FinError("Tau vector must have length" + str(num_periods))

#    if len(zetas) != num_periods:
#        raise FinError("Tau vector must have length" + str(num_periods))

#    if len(rho) != num_periods:
#        raise FinError("Rho matrix must have length" + str(num_periods))

#    if len(rho[0]) != num_periods:
#        raise FinError("Rho matrix must have height" + str(num_periods))

    if b > num_periods:
        raise FinError("Swap maturity beyond num_periods.")

    if a == b:
        raise FinError("Swap maturity on swap expiry date")

    p = np.zeros(num_periods)
    p[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, num_periods):
        p[ix] = p[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    wts = np.zeros(num_periods)
    pv01ab = 0.0
    for k in range(a+1, b):
        pv01ab += taus[k] * p[k]

    sab = (p[a] - p[b-1])/pv01ab

    for i in range(a, b):
        wts[i] = taus[i] * p[i] / pv01ab

    swaption_var = 0.0
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
            swaption_var += term

    tau_a = 0.0
    for i in range(0, a):
        tau_a += taus[i]

    tau_b = 0.0
    for i in range(0, b):
        tau_b += taus[i]

    swaption_vol = np.sqrt(swaption_var/tau_a)
    return swaption_vol


###############################################################################


@njit(float64(int64, int64, float64[:], float64[:, :, :], float64[:]),
      cache=True, fastmath=True)
def lmm_sim_swaption_vol(a, b, fwd0, fwds, taus):
    """ Calculates the swap rate volatility using the forwards generated in the
    simulation to see how it compares to Rebonatto estimate. """

    num_paths = len(fwds)
    num_fwds = len(fwds[0])

    if a > num_fwds:
        raise FinError("NumPeriods > num_fwds")

    if a >= b:
        raise FinError("Swap maturity is before expiry date")

    fwd_swap_rate_mean = 0.0
    fwd_swap_rate_var = 0.0

    for i_path in range(0, num_paths):  # changed from prange

        numeraire = 1.0

        for k in range(0, a):
            numeraire *= (1.0 + taus[k] * fwds[i_path, k, k])

        pv01 = 0.0
        df = 1.0

        for k in range(a, b):

            f = fwds[i_path, a, k]
            tau = taus[k]
            df = df / (1.0 + tau * f)
            pv01 = pv01 + tau * df

        fwd_swap_rate = (1.0 - df) / pv01

        fwd_swap_rate_mean += fwd_swap_rate
        fwd_swap_rate_var += fwd_swap_rate**2

    taua = 0.0
    for i in range(0, a):
        taua += taus[i]

    fwd_swap_rate_mean /= num_paths
    fwd_swap_rate_var = fwd_swap_rate_var/num_paths - fwd_swap_rate_mean**2
    fwd_swap_rate_vol = np.sqrt(fwd_swap_rate_var/taua)
    fwd_swap_rate_vol /= fwd_swap_rate_mean
    return fwd_swap_rate_vol

###############################################################################


@njit(float64[:, :](int64, int64, int64, float64[:, :, :]),
      cache=True, fastmath=True)
def lmm_fwd_fwd_correlation(num_fwds, num_paths, i_time, fwds):
    """ Extract forward forward correlation matrix at some future time index
    from the simulated forward rates and return the matrix. """

    size = num_fwds - i_time
    fwd_corr = np.zeros((size, size))

    for i_fwd in range(i_time, num_fwds):
        for j_fwd in range(i_fwd, num_fwds):

            sumfwdi = 0.0
            sumfwdj = 0.0
            sumfwdifwdi = 0.0
            sumfwdifwdj = 0.0
            sumfwdjfwdj = 0.0

            for p in range(0, num_paths):  # changed from prange
                dfwdi = fwds[p, i_time, i_fwd] - fwds[p, i_time-1, i_fwd]
                dfwdj = fwds[p, i_time, j_fwd] - fwds[p, i_time-1, j_fwd]
                sumfwdi += dfwdi
                sumfwdj += dfwdj
                sumfwdifwdi += dfwdi * dfwdi
                sumfwdifwdj += dfwdi * dfwdj
                sumfwdjfwdj += dfwdj * dfwdj

            avgfwdi = sumfwdi / num_paths
            avgfwdj = sumfwdj / num_paths
            avgfwdifwdi = sumfwdifwdi / num_paths
            avgfwdifwdj = sumfwdifwdj / num_paths
            avgfwdjfwdj = sumfwdjfwdj / num_paths

            covii = avgfwdifwdi - avgfwdi * avgfwdi
            covjj = avgfwdjfwdj - avgfwdj * avgfwdj
            covij = avgfwdifwdj - avgfwdi * avgfwdj
            corr = covij / np.sqrt(covii*covjj)

            if abs(covii*covjj) > 1e-20:
                fwd_corr[i_fwd-i_time][j_fwd-i_time] = corr
                fwd_corr[j_fwd-i_time][i_fwd-i_time] = corr
            else:
                fwd_corr[i_fwd-i_time][j_fwd-i_time] = 0.0
                fwd_corr[j_fwd-i_time][i_fwd-i_time] = 0.0

    return fwd_corr

###############################################################################


@njit(float64[:](float64[:], float64[:], int64, float64, float64[:]),
      cache=True, fastmath=True)
def lmm_price_caps_black(fwd0, vol_caplet, p, k, taus):
    """ Price a strip of capfloorlets using Black's model using the time grid
    of the LMM model. The prices can be compared with the LMM model prices. """

    caplet = np.zeros(p+1)
    disc_fwd = np.zeros(p+1)

    if k <= 0.0:
        raise FinError("Negative strike not allowed.")

    # Set up initial term structure
    disc_fwd[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for i in range(1, p):
        disc_fwd[i] = disc_fwd[i-1] / (1.0 + fwd0[i] * taus[i])

    # Price ATM caplets
    t_exp = 0.0

    for i in range(1, p):  # 1 to p-1

        k = fwd0[i]
        t_exp += taus[i]
        vol = vol_caplet[i]
        f = fwd0[i]
        d1 = (np.log(f/k) + vol * vol * t_exp / 2.0) / vol / np.sqrt(t_exp)
        d2 = d1 - vol * np.sqrt(t_exp)
        caplet[i] = (f * N(d1) - k * N(d2)) * taus[i] * disc_fwd[i]

    return caplet

###############################################################################


@njit(float64[:, :](float64[:, :], int64), cache=True, fastmath=True)
def sub_matrix(t, N):
    """ Returns a submatrix of correlation matrix at later time step in the LMM
    simulation which is then used to generate correlated Gaussian RVs. """

    lent = len(t)
    result = np.zeros((lent-N-1, lent-N-1))

    for i in range(N + 1, lent):
        for j in range(N + 1, lent):
            result[i - N - 1][j - N - 1] = t[i][j]

    return result

###############################################################################


@njit(float64[:, :](float64[:, :]), cache=True, fastmath=True)
def cholesky_np(rho):
    """ Numba-compliant wrapper around Numpy cholesky function. """
    chol = np.linalg.cholesky(rho)
    return chol

###############################################################################


@njit(float64[:, :, :](int64, int64, float64[:], float64[:], float64[:, :],
                       float64[:], int64), cache=True, fastmath=True)
def lmm_simulate_fwds_nf(num_fwds, num_paths, fwd0, zetas, correl, taus, seed):
    """ Full N-Factor Arbitrage-free simulation of forward Ibor discount in the
    spot measure given an initial forward curve, volatility term structure and
    full rank correlation structure. Cholesky decomposition is used to extract
    the factor weights. The number of forwards at time 0 is given. The 3D
    matrix of forward rates by path, time and forward point is returned.
    WARNING: NEED TO CHECK THAT CORRECT VOLATILITY IS BEING USED (OFF BY ONE
    BUG NEEDS TO BE RULED OUT) """

    np.random.seed(seed)

    # Even number of paths for antithetics
    num_paths = 2 * int(num_paths/2)
    half_num_paths = int(num_paths/2)

    fwd = np.empty((num_paths, num_fwds, num_fwds))
    fwd_b = np.zeros(num_fwds)

    disc_fwd = np.zeros(num_fwds)

    # Set up initial term structure
    disc_fwd[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, num_fwds):
        disc_fwd[ix] = disc_fwd[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    corr = [None]  # from 0 to p-1
    factors = [None]  # from 0 to p-1

    for ix in range(1, num_fwds):  # from 1 to p-1
        matrix = sub_matrix(correl, ix - 1)
        corr.append(matrix)
        chol = cholesky_np(matrix)
        factors.append(chol)

    ###########################################################################
    # I HAVE PROBLEMS AS THE PARALLELISATION CHANGES THE OUTPUT IF RANDS ARE
    # CALCULATED INSIDE THE MAIN LOOP SO I CALCULATE THEM NOW
    ###########################################################################

    if 1 == 1:
        g_matrix = np.empty((num_paths, num_fwds, num_fwds))
        for i_path in range(0, half_num_paths):
            for j in range(1, num_fwds):
                for k in range(0, num_fwds-j):
                    g = np.random.normal()
                    # ANTITHETICS
                    g_matrix[i_path, j, k] = g
                    g_matrix[i_path + half_num_paths, j, k] = -g

    avgg = 0.0
    stdg = 0.0

    for i_path in range(0, num_paths):

        # Initial value of forward curve at time 0
        for i_fwd in range(0, num_fwds):
            fwd[i_path, 0, i_fwd] = fwd0[i_fwd]

        for j in range(1, num_fwds):  # TIME LOOP

            dt = taus[j]
            sqrt_dt = np.sqrt(dt)

            for i in range(j, num_fwds):  # FORWARDS LOOP

                zi = zetas[i]

                mu_a = 0.0
                for k in range(j, i+1):
                    rho = corr[j][k-j, i-j]
                    fk = fwd[i_path, j-1, k]
                    zk = zetas[k]
                    tk = taus[k]
                    mu_a += zi * fk * tk * zk * rho / (1.0 + fk * tk)

                w = 0.0
                for k in range(0, num_fwds-j):
                    f = factors[j][i-j, k]
                    w = w + f * g_matrix[i_path, j, k]

                avgg += w
                stdg += w*w

                fwd_b[i] = fwd[i_path, j-1, i] \
                    * np.exp(mu_a * dt - 0.5 * (zi**2) * dt + zi * w * sqrt_dt)

                mu_b = 0.0
                for k in range(j, i+1):
                    rho = corr[j][k-j, i-j]
                    fk = fwd_b[k]
                    zk = zetas[k]
                    tk = taus[k]
                    mu_b += zi * fk * tk * zk * rho / (1.0 + fk * tk)

                mu_avg = 0.5*(mu_a + mu_b)
                x = np.exp(mu_avg * dt - 0.5 * (zi**2) * dt + zi * w * sqrt_dt)
                fwd[i_path, j, i] = fwd[i_path, j-1, i] * x

    return fwd

###############################################################################


@njit(float64[:, :, :](int64, int64, int64, float64[:], float64[:], float64[:],
                       int64, int64), cache=True, fastmath=True)
def lmm_simulate_fwds_1f(num_fwds, num_paths, numeraire_index, fwd0, gammas,
                         taus, use_sobol, seed):
    """ One factor Arbitrage-free simulation of forward Ibor discount in the
    spot measure following Hull Page 768. Given an initial forward curve,
    volatility term structure. The 3D matrix of forward rates by path, time
    and forward point is returned. This function is kept mainly for its
    simplicity and speed.

    NB: The Gamma volatility has an initial entry of zero. This differs from
    Hull's indexing by one and so is why I do not subtract 1 from the index as
    Hull does in his equation 32.14.

    The Number of Forwards is the number of points on the initial curve to the
    trade maturity date.

    But be careful: a cap that matures in 10 years with quarterly caplets has
    40 forwards BUT the last forward to reset occurs at 9.75 years. You should
    not simulate beyond this time. If you give the model 10 years as in the
    Hull examples, you need to simulate 41 (or in this case 11) forwards as the
    final cap or ratchet has its reset in 10 years. """

    if len(gammas) != num_fwds:
        raise FinError("Gamma vector does not have right number of forwards")

    if len(fwd0) != num_fwds:
        raise FinError("The length of fwd0 is not equal to num_fwds")

    if len(taus) != num_fwds:
        raise FinError("The length of Taus is not equal to num_fwds")

    np.random.seed(seed)
    # Even number of paths for antithetics
    num_paths = 2 * int(num_paths/2)
    half_num_paths = int(num_paths/2)
    fwd = np.empty((num_paths, num_fwds, num_fwds))
    fwdB = np.zeros(num_fwds)

    num_times = num_fwds

    if use_sobol == 1:
        num_dimensions = num_times
        rands = get_uniform_sobol(half_num_paths, num_dimensions)
        g_matrix = np.empty((num_paths, num_times))
        for i_path in range(0, half_num_paths):
            for j in range(0, num_times):
                u = rands[i_path, j]
                g = norminvcdf(u)
                g_matrix[i_path, j] = g
                g_matrix[i_path + half_num_paths, j] = -g
    elif use_sobol == 0:
        g_matrix = np.empty((num_paths, num_times))
        for i_path in range(0, half_num_paths):
            for j in range(0, num_times):
                g = np.random.normal()
                g_matrix[i_path, j] = g
                g_matrix[i_path + half_num_paths, j] = -g
    else:
        raise FinError("Use Sobol must be 0 or 1")

    for i_path in range(0, num_paths):  # changed from prange
        # Initial value of forward curve at time 0
        for i_fwd in range(0, num_fwds):
            fwd[i_path, 0, i_fwd] = fwd0[i_fwd]

        for j in range(0, num_fwds-1):  # TIME LOOP
            dtj = taus[j]
            sqrt_dtj = np.sqrt(dtj)
            w = g_matrix[i_path, j]

            for k in range(j, num_fwds):  # FORWARDS LOOP
                zkj = gammas[k-j]
                muA = 0.0

                for i in range(j+1, k+1):
                    fi = fwd[i_path, j, i]
                    zij = gammas[i-j]
                    ti = taus[i]
                    muA += zkj * fi * ti * zij / (1.0 + fi * ti)

                # predictor corrector
                x = np.exp(muA * dtj - 0.5*(zkj**2) * dtj + zkj * w * sqrt_dtj)
                fwdB[k] = fwd[i_path, j, k] * x

                muB = 0.0
                for i in range(j+1, k+1):
                    fi = fwdB[k]
                    zij = gammas[i-j]
                    ti = taus[i]
                    muB += zkj * fi * ti * zij / (1.0 + fi * ti)

                muC = 0.5*(muA+muB)

                x = np.exp(muC*dtj - 0.5 * (zkj**2) * dtj + zkj * w * sqrt_dtj)
                fwd[i_path, j+1, k] = fwd[i_path, j, k] * x

    return fwd

###############################################################################


@njit(float64[:, :, :](int64, int64, int64, int64, float64[:], float64[:, :],
                       float64[:], int64, int64), cache=True, fastmath=True)
def lmm_simulate_fwds_mf(num_fwds, num_factors, num_paths, numeraire_index,
                         fwd0, lambdas, taus, use_sobol, seed):
    """ Multi-Factor Arbitrage-free simulation of forward Ibor discount in the
    spot measure following Hull Page 768. Given an initial forward curve,
    volatility factor term structure. The 3D matrix of forward rates by path,
    time and forward point is returned. """

    np.random.seed(seed)

    if len(lambdas) != num_factors:
        raise FinError("Lambda does not have the right number of factors")

    if len(lambdas[0]) != num_fwds:
        raise FinError("Lambda does not have the right number of forwards")

    # Even number of paths for antithetics
    num_paths = 2 * int(num_paths/2)
    half_num_paths = int(num_paths/2)
    fwd = np.empty((num_paths, num_fwds, num_fwds))
    fwdB = np.zeros(num_fwds)

    num_times = num_fwds

    if use_sobol == 1:
        num_dimensions = num_times * num_factors
        rands = get_uniform_sobol(half_num_paths, num_dimensions)
        g_matrix = np.empty((num_paths, num_times, num_factors))
        for i_path in range(0, half_num_paths):
            for j in range(0, num_times):
                for q in range(0, num_factors):
                    col = j*num_factors + q
                    u = rands[i_path, col]
                    g = norminvcdf(u)
                    g_matrix[i_path, j, q] = g
                    g_matrix[i_path + half_num_paths, j, q] = -g
    elif use_sobol == 0:
        g_matrix = np.empty((num_paths, num_times, num_factors))
        for i_path in range(0, half_num_paths):
            for j in range(0, num_times):
                for q in range(0, num_factors):
                    g = np.random.normal()
                    g_matrix[i_path, j, q] = g
                    g_matrix[i_path + half_num_paths, j, q] = -g
    else:
        raise FinError("Use Sobol must be 0 or 1.")

    for i_path in range(0, num_paths):
        # Initial value of forward curve at time 0
        for i_fwd in range(0, num_fwds):
            fwd[i_path, 0, i_fwd] = fwd0[i_fwd]

        for j in range(0, num_fwds-1):  # TIME LOOP
            dtj = taus[j]
            sqrt_dtj = np.sqrt(dtj)

            for k in range(j, num_fwds):  # FORWARDS LOOP

                muA = 0.0
                for i in range(j+1, k+1):
                    fi = fwd[i_path, j, i]
                    ti = taus[i]
                    zz = 0.0
                    for q in range(0, num_factors):
                        zij = lambdas[q][i-j]
                        zkj = lambdas[q][k-j]
                        zz += zij * zkj
                    muA += fi * ti * zz / (1.0 + fi * ti)

                itoTerm = 0.0
                for q in range(0, num_factors):
                    itoTerm += lambdas[q][k-j] * lambdas[q][k-j]

                random_term = 0.0
                for q in range(0, num_factors):
                    wq = g_matrix[i_path, j, q]
                    random_term += lambdas[q][k-j] * wq
                random_term *= sqrt_dtj

                x = np.exp(muA * dtj - 0.5 * itoTerm * dtj + random_term)
                fwdB[k] = fwd[i_path, j, k] * x

                muB = 0.0
                for i in range(j+1, k+1):
                    fi = fwdB[k]
                    ti = taus[i]
                    zz = 0.0
                    for q in range(0, num_factors):
                        zij = lambdas[q][i-j]
                        zkj = lambdas[q][k-j]
                        zz += zij * zkj
                    muB += fi * ti * zz / (1.0 + fi * ti)

                muC = 0.5 * (muA + muB)

                x = np.exp(muC * dtj - 0.5 * itoTerm * dtj + random_term)
                fwd[i_path, j+1, k] = fwd[i_path, j, k] * x

    return fwd

###############################################################################


@njit(float64[:](int64, int64, float64, float64[:], float64[:, :, :],
                 float64[:], int64),
      cache=True, fastmath=True)
def lmm_cap_flr_pricer(num_fwds, num_paths, K, fwd0, fwds, taus, is_cap):
    """ Function to price a strip of cap or floorlets in accordance with the
    simulated forward curve dynamics. """

    max_paths = len(fwds)
    max_fwds = len(fwds[0])

    if num_fwds > max_fwds:
        raise FinError("num_fwds > max_fwds")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(num_fwds)
    capFlrLets = np.zeros(num_fwds-1)
    capFlrLetValues = np.zeros(num_fwds-1)
    numeraire = np.zeros(num_fwds)

    for i_path in range(0, num_paths):

        period_roll = 1.0
        libor = fwds[i_path, 0, 0]
        capFlrLets[0] = max(K - libor, 0.0) * taus[0]

        # Now loop over the caplets starting with one that fixes immediately
        # but which may have intrinsic value that cannot be ignored.
        for j in range(0, num_fwds):

            libor = fwds[i_path, j, j]
            if j == 1:
                if is_cap == 0:
                    capFlrLets[j] = max(K - libor, 0.0) * taus[j]
                else:
                    capFlrLets[j] = max(libor - K, 0.0) * taus[j]

                numeraire[0] = 1.0 / df[0]
            else:
                if is_cap == 1:
                    capFlrLets[j] = max(libor - K, 0.0) * taus[j]
                elif is_cap == 0:
                    capFlrLets[j] = max(K - libor, 0.0) * taus[j]
                else:
                    raise FinError("is_cap should be 0 or 1")

            period_roll = 1.0 + libor * taus[j]
            numeraire[j] = numeraire[j - 1] * period_roll

        for i_fwd in range(0, num_fwds):
            denom = abs(numeraire[i_fwd]) + 1e-12
            capFlrLetValues[i_fwd] += capFlrLets[i_fwd] / denom

    for i_fwd in range(0, num_fwds):
        capFlrLetValues[i_fwd] /= num_paths

    return capFlrLetValues

###############################################################################


@njit(float64(float64, int64, int64, float64[:], float64[:, :, :],
              float64[:]), cache=True, fastmath=True)
def lmm_swap_pricer(cpn, num_periods, num_paths, fwd0, fwds, taus):
    """ Function to reprice a basic swap using the simulated forward Ibors.
    """

    max_paths = len(fwds)
    max_fwds = len(fwds[0])

    if num_periods > max_fwds:
        raise FinError("NumPeriods > num_fwds")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(max_fwds)
    numeraire = np.zeros(max_fwds)
    sum_fixed = 0.0
    sun_float = 0.0
    fixed_flows = np.zeros(max_fwds)
    float_flows = np.zeros(max_fwds)

    # Set up initial term structure
    df[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, max_fwds):
        df[ix] = df[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for i_path in range(0, num_paths):

        period_roll = 1.0
        libor = fwds[i_path, 0, 0]
        float_flows[0] = libor * taus[0]
        fixed_flows[0] = cpn * taus[0]
        numeraire[0] = 1.0 / df[0]

        for j in range(1, num_periods):  # TIME LOOP

            libor = fwds[i_path, j, j]

            if j == 1:
                fixed_flows[j] = cpn * taus[j]
                float_flows[j] = libor * taus[j]
            else:
                fixed_flows[j] = fixed_flows[j-1] * period_roll + cpn * taus[j]
                float_flows[j] = float_flows[j-1] * \
                    period_roll + libor * taus[j]

            period_roll = 1.0 + libor * taus[j]
            numeraire[j] = numeraire[j - 1] * period_roll

        for i_fwd in range(0, num_periods):
            sun_float += float_flows[i_fwd] / numeraire[i_fwd]
            sum_fixed += fixed_flows[i_fwd] / numeraire[i_fwd]

    sun_float /= num_paths
    sum_fixed /= num_paths
    v = sum_fixed - sun_float
    pv01 = sum_fixed/cpn
    swap_rate = sun_float/pv01

    print("FLOAT LEG:", sun_float)
    print("FIXED LEG:", sum_fixed)
    print("SWAP RATE:", swap_rate)
    print("NET VALUE:", v)
    return v

###############################################################################


@njit(float64(float64, int64, int64, int64, float64[:], float64[:, :, :],
              float64[:], int64), cache=True, fastmath=True)
def lmm_swaption_pricer(strike, a, b, num_paths, fwd0, fwds, taus, is_payer):
    """ Function to price a European swaption using the simulated forward
    discount. """

    max_paths = len(fwds)
    max_fwds = len(fwds[0])

    if a > max_fwds:
        raise FinError("NumPeriods > num_fwds")

    if a >= b:
        raise FinError("Swap maturity is before expiry date")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(max_fwds)
#    pv01 = np.zeros(max_fwds)

    # Set up initial term structure
    df[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, b):
        df[ix] = df[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    sumPayRecSwaption = 0.0

    for i_path in range(0, num_paths):

        numeraire = 1.0
        for k in range(0, a):
            numeraire *= (1.0 + taus[k] * fwds[i_path, k, k])

        pv01 = 0.0
        df = 1.0

        # Value the swap as if we were at time a with forward curve known
        for k in range(a, b):
            f = fwds[i_path, a, k]
            tau = taus[k]
            df = df / (1.0 + tau * f)
            pv01 = pv01 + tau * df

        fwd_swap_rate = (1.0 - df) / pv01

        if is_payer == 1:
            payRecSwaption = max(fwd_swap_rate - strike, 0.0) * pv01
        elif is_payer == 0:
            payRecSwaption = max(strike - fwd_swap_rate, 0.0) * pv01
        else:
            raise FinError("Unknown payRecSwaption value - must be 0 or 1")

        sumPayRecSwaption += payRecSwaption / (abs(numeraire) + 1e-10)

    payRecPrice = sumPayRecSwaption / num_paths
    return payRecPrice

###############################################################################


@njit(float64[:](float64, int64, int64, float64[:], float64[:, :, :],
                 float64[:]), cache=True, fastmath=True)
def lmm_ratchet_caplet_pricer(spd, num_periods, num_paths, fwd0, fwds, taus):
    """ Price a ratchet using the simulated Ibor rates."""

    max_paths = len(fwds)
    max_fwds = len(fwds[0][0])

    if num_periods > max_fwds:
        raise FinError("NumPeriods > num_fwds")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(max_fwds)
    numeraire = np.zeros(max_fwds)
    rachet_caplets = np.zeros(max_fwds)
    rachet_caplet_values = np.zeros(max_fwds)

    # Set up initial term structure
    df[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, max_fwds):
        df[ix] = df[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for i_path in range(0, num_paths):

        period_roll = 1.0
        libor = fwds[i_path, 0, 0]
        rachet_caplets[0] = 0.0

        for j in range(1, num_periods):  # TIME LOOP

            prevIbor = libor
            K = prevIbor + spd
            libor = fwds[i_path, j, j]

            if j == 1:
                rachet_caplets[j] = max(libor - K, 0.0) * taus[j]
                numeraire[0] = 1.0 / df[0]
            else:
                rachet_caplets[j] = max(libor - K, 0.0) * taus[j]

            period_roll = 1.0 + libor * taus[j]
            numeraire[j] = numeraire[j - 1] * period_roll

        for i_fwd in range(0, num_periods):
            rachet_caplet_values[i_fwd] += rachet_caplets[i_fwd] / \
                numeraire[i_fwd]

    for i_fwd in range(0, num_periods):
        rachet_caplet_values[i_fwd] /= num_paths

    return rachet_caplet_values

###############################################################################


@njit(float64(int64, float64, int64, int64, float64[:], float64[:, :, :],
              float64[:]), cache=True, fastmath=True)
def lmm_flexi_cap_pricer(maxCaplets, K, num_periods, num_paths,
                         fwd0, fwds, taus):
    """ Price a flexicap using the simulated Ibor rates."""

    max_paths = len(fwds)
    max_fwds = len(fwds[0][0])

    if num_periods > max_fwds:
        raise FinError("NumPeriods > num_fwds")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(max_fwds)
    numeraire = np.zeros(max_fwds)
    flexi_caplets = np.zeros(max_fwds)
    flexi_caplet_values = np.zeros(max_fwds)

    # Set up initial term structure
    df[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, max_fwds):
        df[ix] = df[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for i_path in range(0, num_paths):

        period_roll = 1.0
        libor = fwds[i_path, 0, 0]
        flexi_caplets[0] = 0.0

        num_caplets_left = maxCaplets

        for j in range(1, num_periods):  # TIME LOOP

            libor = fwds[i_path, j, j]

            if j == 1:
                if libor > K and num_caplets_left > 0:
                    flexi_caplets[j] = max(libor - K, 0.0) * taus[j]
                    num_caplets_left -= 1
                numeraire[0] = 1.0 / df[0]
            else:
                if libor > K and num_caplets_left > 0:
                    flexi_caplets[j] = max(libor - K, 0.0) * taus[j]
                    num_caplets_left -= 1

            period_roll = 1.0 + libor * taus[j]
            numeraire[j] = numeraire[j - 1] * period_roll

        for i_fwd in range(0, num_periods):
            flexi_caplet_values[i_fwd] += flexi_caplets[i_fwd] / numeraire[i_fwd]

    for i_fwd in range(0, num_periods):
        flexi_caplet_values[i_fwd] /= num_paths

    flexi_cap_value = 0.0
    for i_fwd in range(0, num_periods):
        flexi_cap_value += flexi_caplet_values[i_fwd]

    return flexi_cap_value

###############################################################################


@njit(float64[:](float64, int64, int64, float64[:], float64[:, :, :],
                 float64[:]), cache=True, fastmath=True)
def lmm_sticky_caplet_pricer(spread, num_periods, num_paths, fwd0, fwds, taus):
    """ Price a sticky cap using the simulated Ibor rates. """

    max_paths = len(fwds)
    max_fwds = len(fwds[0][0])

    if num_periods > max_fwds:
        raise FinError("NumPeriods > num_fwds")

    if num_paths > max_paths:
        raise FinError("NumPaths > MaxPaths")

    df = np.zeros(max_fwds)
    numeraire = np.zeros(max_fwds)
    stickyCaplets = np.zeros(max_fwds)
    stickyCapletValues = np.zeros(max_fwds)

    # Set up initial term structure
    df[0] = 1.0 / (1.0 + fwd0[0] * taus[0])
    for ix in range(1, max_fwds):
        df[ix] = df[ix-1] / (1.0 + fwd0[ix] * taus[ix])

    for i_path in range(0, num_paths):

        period_roll = 1.0
        libor = fwds[i_path, 0, 0]
        stickyCaplets[0] = 0.0
        K = libor

        for j in range(1, num_periods):  # TIME LOOP

            prevIbor = libor
            K = min(prevIbor, K) + spread
            libor = fwds[i_path, j, j]

            if j == 1:
                stickyCaplets[j] = max(libor-K, 0.0) * taus[j]
                numeraire[0] = 1.0 / df[0]
            else:
                stickyCaplets[j] = max(libor - K, 0.0) * taus[j]

            period_roll = (1.0 + libor * taus[j])
            numeraire[j] = numeraire[j - 1] * period_roll

        for i_fwd in range(0, num_periods):
            stickyCapletValues[i_fwd] += stickyCaplets[i_fwd] / \
                numeraire[i_fwd]

    for i_fwd in range(0, num_periods):
        stickyCapletValues[i_fwd] /= num_paths

    return stickyCapletValues

###############################################################################
