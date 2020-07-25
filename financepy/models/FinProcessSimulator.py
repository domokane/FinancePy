##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, exp, log
from enum import Enum

from numba import njit, float64, int64
import numpy as np

from ..finutils.FinError import FinError
from ..finutils.FinMath import norminvcdf
from ..finutils.FinHelperFunctions import labelToString

###############################################################################


class FinProcessTypes(Enum):
    GBM = 1
    CIR = 2
    HESTON = 3
    VASICEK = 4
    CEV = 5
    JUMP_DIFFUSION = 6

###############################################################################


class FinProcessSimulator():

    def __init__(self):
        pass

    def getProcess(
            self,
            processType,
            t,
            modelParams,
            numAnnSteps,
            numPaths,
            seed):

        if processType == FinProcessTypes.GBM:
            (stockPrice, drift, volatility, scheme) = modelParams
            paths = getGBMPaths(numPaths, numAnnSteps, t, drift,
                                stockPrice, volatility, scheme.value, seed)
            return paths

        elif processType == FinProcessTypes.HESTON:
            (stockPrice, drift, v0, kappa, theta, sigma, rho, scheme) = modelParams
            paths = getHestonPaths(numPaths,
                                   numAnnSteps,
                                   t,
                                   drift,
                                   stockPrice,
                                   v0,
                                   kappa,
                                   theta,
                                   sigma,
                                   rho,
                                   scheme.value,
                                   seed)
            return paths

        elif processType == FinProcessTypes.VASICEK:
            (r0, kappa, theta, sigma, scheme) = modelParams
            paths = getVasicekPaths(
                numPaths,
                numAnnSteps,
                t,
                r0,
                kappa,
                theta,
                sigma,
                scheme.value,
                seed)
            return paths

        elif processType == FinProcessTypes.CIR:
            (r0, kappa, theta, sigma, scheme) = modelParams
            paths = getCIRPaths(numPaths, numAnnSteps, t,
                                r0, kappa, theta, sigma, scheme.value, seed)
            return paths

        else:
            raise FinError("Unknown process" + str(processType))

###############################################################################


class FinHestonNumericalScheme(Enum):
    EULER = 1
    EULERLOG = 2
    QUADEXP = 3

###############################################################################


@njit(float64[:, :](int64, int64, float64, float64, float64, float64, float64,
                    float64, float64, float64, int64, int64),
      cache=True, fastmath=True)
def getHestonPaths(numPaths,
                   numAnnSteps,
                   t,
                   drift,
                   s0,
                   v0,
                   kappa,
                   theta,
                   sigma,
                   rho,
                   scheme,
                   seed):

    np.random.seed(seed)
    dt = 1.0 / numAnnSteps
    numSteps = int(t / dt)
    sPaths = np.empty(shape=(numPaths, numSteps + 1))
    sPaths[:, 0] = s0
    sdt = sqrt(dt)
    rhohat = sqrt(1.0 - rho * rho)
    sigma2 = sigma * sigma

    if scheme == FinHestonNumericalScheme.EULER.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0, numPaths):
            s = s0
            v = v0
            for iStep in range(1, numSteps + 1):
                z1 = np.random.normal(0.0, 1.0) * sdt
                z2 = np.random.normal(0.0, 1.0) * sdt
                zV = z1
                zS = rho * z1 + rhohat * z2
                vplus = max(v, 0.0)
                rtvplus = sqrt(vplus)
                v += kappa * (theta - vplus) * dt + sigma * \
                    rtvplus * zV + 0.25 * sigma2 * (zV * zV - dt)
                s += drift * s * dt + rtvplus * s * \
                    zS + 0.5 * s * vplus * (zV * zV - dt)
                sPaths[iPath, iStep] = s

    elif scheme == FinHestonNumericalScheme.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0, numPaths):
            x = log(s0)
            v = v0
            for iStep in range(1, numSteps + 1):
                zV = np.random.normal(0.0, 1.0) * sdt
                zS = rho * zV + rhohat * np.random.normal(0.0, 1.0) * sdt
                vplus = max(v, 0.0)
                rtvplus = sqrt(vplus)
                x += (drift - 0.5 * vplus) * dt + rtvplus * zS
                v += kappa * (theta - vplus) * dt + sigma * \
                    rtvplus * zV + sigma2 * (zV * zV - dt) / 4.0
                sPaths[iPath, iStep] = exp(x)

    elif scheme == FinHestonNumericalScheme.QUADEXP.value:
        # Due to Leif Andersen(2006)
        Q = exp(-kappa * dt)
        psic = 1.50
        gamma1 = 0.50
        gamma2 = 0.50
        K0 = -rho * kappa * theta * dt / sigma
        K1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho / sigma
        K2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho / sigma
        K3 = gamma1 * dt * (1.0 - rho * rho)
        K4 = gamma2 * dt * (1.0 - rho * rho)
        A = K2 + 0.5 * K4
        mu = drift
        c1 = sigma2 * Q * (1.0 - Q) / kappa
        c2 = theta * sigma2 * ((1.0 - Q)**2) / 2.0 / kappa

        for iPath in range(0, numPaths):
            x = log(s0)
            vn = v0
            for iStep in range(1, numSteps + 1):
                zV = np.random.normal(0, 1)
                zS = rho * zV + rhohat * np.random.normal(0, 1)
                m = theta + (vn - theta) * Q
                m2 = m * m
                s2 = c1 * vn + c2
                psi = s2 / m2
                u = np.random.uniform(0.0, 1.0)

                if psi <= psic:
                    b2 = 2.0 / psi - 1.0 + \
                        sqrt((2.0 / psi) * (2.0 / psi - 1.0))
                    a = m / (1.0 + b2)
                    b = sqrt(b2)
                    zV = norminvcdf(u)
                    vnp = a * ((b + zV)**2)
                    d = (1.0 - 2.0 * A * a)
                    M = exp((A * b2 * a) / d) / sqrt(d)
                    K0 = -log(M) - (K1 + 0.5 * K3) * vn
                else:
                    p = (psi - 1.0) / (psi + 1.0)
                    beta = (1.0 - p) / m

                    if u <= p:
                        vnp = 0.0
                    else:
                        vnp = log((1.0 - p) / (1.0 - u)) / beta

                    M = p + beta * (1.0 - p) / (beta - A)
                    K0 = -log(M) - (K1 + 0.5 * K3) * vn

                x += mu * dt + K0 + (K1 * vn + K2 * vnp) + \
                    sqrt(K3 * vn + K4 * vnp) * zS
                sPaths[iPath, iStep] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return sPaths

###############################################################################


class FinGBMNumericalScheme(Enum):
    NORMAL = 1
    ANTITHETIC = 2

###############################################################################

@njit(float64[:, :](int64, int64, float64, float64, float64,
                    float64, int64, int64), cache=True, fastmath=True)
def getGBMPaths(numPaths, numAnnSteps, t, mu, stockPrice, sigma, scheme, seed):

    np.random.seed(seed)
    dt = 1.0 / numAnnSteps
    numTimeSteps = int(t / dt + 0.50)
    vsqrtdt = sigma * sqrt(dt)
    m = exp((mu - sigma * sigma / 2.0) * dt)

    if scheme == FinGBMNumericalScheme.NORMAL.value:

        Sall = np.empty((numPaths, numTimeSteps + 1))
        Sall[:, 0] = stockPrice
        for it in range(1, numTimeSteps + 1):
            g1D = np.random.standard_normal((numPaths))
            for ip in range(0, numPaths):
                w = np.exp(g1D[ip] * vsqrtdt)
                Sall[ip, it] = Sall[ip, it - 1] * m * w

    elif scheme == FinGBMNumericalScheme.ANTITHETIC.value:

        Sall = np.empty((2 * numPaths, numTimeSteps + 1))
        Sall[:, 0] = stockPrice
        for it in range(1, numTimeSteps + 1):
            g1D = np.random.standard_normal((numPaths))
            for ip in range(0, numPaths):
                w = np.exp(g1D[ip] * vsqrtdt)
                Sall[ip, it] = Sall[ip, it - 1] * m * w
                Sall[ip + numPaths, it] = Sall[ip + numPaths, it - 1] * m / w

    else:

        raise FinError("Unknown FinGBMNumericalScheme")

#    m = np.mean(Sall[:, -1])
#    v = np.var(Sall[:, -1]/Sall[:, 0])
#    print("GBM", numPaths, numAnnSteps, t, mu, stockPrice, sigma, scheme, m,v)

    return Sall

###############################################################################


class FinVasicekNumericalScheme(Enum):
    NORMAL = 1
    ANTITHETIC = 2

###############################################################################

@njit(float64[:, :](int64, int64, float64, float64, float64,
                    float64, float64, int64, int64), cache=True, fastmath=True)
def getVasicekPaths(numPaths,
                    numAnnSteps,
                    t,
                    r0,
                    kappa,
                    theta,
                    sigma,
                    scheme,
                    seed):

    np.random.seed(seed)
    dt = 1.0 / numAnnSteps
    numSteps = int(t / dt)
    sigmasqrtdt = sigma * sqrt(dt)

    if scheme == FinVasicekNumericalScheme.NORMAL.value:
        ratePath = np.empty((numPaths, numSteps + 1))
        ratePath[:, 0] = r0
        for iPath in range(0, numPaths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                r += kappa * (theta - r) * dt + z[iStep - 1] * sigmasqrtdt
                ratePath[iPath, iStep] = r
    elif scheme == FinVasicekNumericalScheme.ANTITHETIC.value:
        ratePath = np.empty((2 * numPaths, numSteps + 1))
        ratePath[:, 0] = r0
        for iPath in range(0, numPaths):
            r1 = r0
            r2 = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                r1 = r1 + kappa * (theta - r1) * dt + \
                    z[iStep - 1] * sigmasqrtdt
                r2 = r2 + kappa * (theta - r2) * dt - \
                    z[iStep - 1] * sigmasqrtdt
                ratePath[iPath, iStep] = r1
                ratePath[iPath + numPaths, iStep] = r2
    return ratePath

###############################################################################


class FinCIRNumericalScheme(Enum):
    EULER = 1
    LOGNORMAL = 2
    MILSTEIN = 3
    KAHLJACKEL = 4
    EXACT = 5  # SAMPLES EXACT DISTRIBUTION

###############################################################################

@njit(float64[:, :](int64, int64, float64, float64, float64,
                    float64, float64, int64, int64), cache=True, fastmath=True)
def getCIRPaths(numPaths,
                numAnnSteps,
                t,
                r0,
                kappa,
                theta,
                sigma,
                scheme,
                seed):

    np.random.seed(seed)
    dt = 1.0 / numAnnSteps
    numSteps = int(t / dt)
    ratePath = np.empty(shape=(numPaths, numSteps + 1))
    ratePath[:, 0] = r0

    if scheme == FinCIRNumericalScheme.EULER.value:
        sigmasqrtdt = sigma * sqrt(dt)
        for iPath in range(0, numPaths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                rplus = max(r, 0.0)
                sqrtrplus = sqrt(rplus)
                r = r + kappa * (theta - rplus) * dt + \
                    sigmasqrtdt * z[iStep - 1] * sqrtrplus
                ratePath[iPath, iStep] = r

    elif scheme == FinCIRNumericalScheme.LOGNORMAL.value:
        x = exp(-kappa * dt)
        y = 1.0 - x
        for iPath in range(0, numPaths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                mean = x * r + theta * y
                var = sigma * sigma * y * (x * r + 0.50 * theta * y) / kappa
                sig = sqrt(log(1.0 + var / (mean * mean)))
                r = mean * exp(-0.5 * sig * sig + sig * z[iStep - 1])
                ratePath[iPath, iStep] = r

    elif scheme == FinCIRNumericalScheme.MILSTEIN.value:
        sigmasqrtdt = sigma * sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0
        for iPath in range(0, numPaths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                sqrtrplus = sqrt(max(r, 0.0))
                r = r + kappa * (theta - r) * dt + \
                    z[iStep - 1] * sigmasqrtdt * sqrtrplus
                r = r + sigma2dt * (z[iStep - 1]**2 - 1.0)
                ratePath[iPath, iStep] = r

    elif scheme == FinCIRNumericalScheme.KAHLJACKEL.value:
        bhat = theta - sigma * sigma / 4.0 / kappa
        sqrtdt = sqrt(dt)
        for iPath in range(0, numPaths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps))
            for iStep in range(1, numSteps + 1):
                beta = z[iStep - 1] / sqrtdt
                sqrtrplus = sqrt(max(r, 0.0))
                c = 1.0 + (sigma * beta - 2.0 * kappa *
                           sqrtrplus) * dt / 4.0 / sqrtrplus
                r = r + (kappa * (bhat - r) + sigma *
                         beta * sqrtrplus) * c * dt
                ratePath[iPath, iStep] = r

    return ratePath

###############################################################################
