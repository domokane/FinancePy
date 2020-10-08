##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
from scipy import integrate
from math import exp, log, pi
import numpy as np  # I USE NUMPY FOR EXP, LOG AND SQRT AS THEY HANDLE IMAGINARY PARTS

from ..finutils.FinGlobalVariables import gDaysInYear
from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinMath import norminvcdf
from ..finutils.FinError import FinError

##########################################################################
# Heston Process
# dS = rS dt + sqrt(V) * S * dz
# dV = kappa(theta-V) dt + sigma sqrt(V) dz
# corr(dV,dS) = rho dt
# Rewritten as
# dS = rS dt + sqrt(V) * S * (rhohat dz1 + rho dz2)
# dV = kappa(theta-V) dt + sigma sqrt(V) dz2
# where rhohat = sqrt(1-rho*rho)
###############################################################################
# TODO - DECIDE WHETHER TO OO MODEL
# TODO - NEEDS CHECKING FOR MC CONVERGENCE
###############################################################################

from enum import Enum


class FinHestonNumericalScheme(Enum):
    EULER = 1
    EULERLOG = 2
    QUADEXP = 3

###############################################################################


@njit(float64[:, :](float64, float64, float64, float64, float64, float64,
                    float64, float64, float64, float64, int64, int64, int64),
      cache=True, fastmath=True)
def getPaths(s0, r, q, v0, kappa, theta, sigma, rho, t, dt, numPaths,
             seed, scheme):

    np.random.seed(seed)
    numSteps = int(t / dt)
    sPaths = np.zeros(shape=(numPaths, numSteps))
    sPaths[:, 0] = s0
    sdt = np.sqrt(dt)
    rhohat = np.sqrt(1.0 - rho * rho)
    sigma2 = sigma * sigma

    if scheme == FinHestonNumericalScheme.EULER.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0, numPaths):
            s = s0
            v = v0
            for iStep in range(1, numSteps):
                z1 = np.random.normal(0.0, 1.0) * sdt
                z2 = np.random.normal(0.0, 1.0) * sdt
                zV = z1
                zS = rho * z1 + rhohat * z2
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                v += kappa * (theta - vplus) * dt + sigma * \
                    rtvplus * zV + 0.25 * sigma2 * (zV * zV - dt)
                s += (r - q) * s * dt + rtvplus * s * \
                    zS + 0.5 * s * vplus * (zV * zV - dt)
                sPaths[iPath, iStep] = s

    elif scheme == FinHestonNumericalScheme.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0, numPaths):
            x = log(s0)
            v = v0
            for iStep in range(1, numSteps):
                zV = np.random.normal(0.0, 1.0) * sdt
                zS = rho * zV + rhohat * np.random.normal(0.0, 1.0) * sdt
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                x += (r - q - 0.5 * vplus) * dt + rtvplus * zS
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
        mu = (r - q)
        c1 = sigma2 * Q * (1.0 - Q) / kappa
        c2 = theta * sigma2 * ((1.0 - Q)**2) / 2.0 / kappa

        for iPath in range(0, numPaths):
            x = log(s0)
            vn = v0
            for iStep in range(1, numSteps):
                zV = np.random.normal(0, 1)
                zS = rho * zV + rhohat * np.random.normal(0, 1)
                m = theta + (vn - theta) * Q
                m2 = m * m
                s2 = c1 * vn + c2
                psi = s2 / m2
                u = np.random.uniform(0.0, 1.0)

                if psi <= psic:
                    b2 = 2.0 / psi - 1.0 + \
                        np.sqrt((2.0 / psi) * (2.0 / psi - 1.0))
                    a = m / (1.0 + b2)
                    b = np.sqrt(b2)
                    zV = norminvcdf(u)
                    vnp = a * ((b + zV)**2)
                    d = (1.0 - 2.0 * A * a)
                    M = exp((A * b2 * a) / d) / np.sqrt(d)
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
                    np.sqrt(K3 * vn + K4 * vnp) * zS
                sPaths[iPath, iStep] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return sPaths

###############################################################################


class FinModelHeston():

    def __init__(self, v0, kappa, theta, sigma, rho):

        verbose = False

        if 2.0 * kappa * theta <= sigma and verbose:
            print("Feller condition not satisfied. Zero Variance possible")

        self._v0 = v0
        self._kappa = kappa
        self._theta = theta
        self._sigma = sigma
        self._rho = rho

###############################################################################

    def value_MC(self,
                 valueDate,
                 option,
                 stockPrice,
                 interestRate,
                 dividendYield,
                 numPaths,
                 numStepsPerYear,
                 seed,
                 scheme=FinHestonNumericalScheme.EULERLOG):

        tau = (option._expiryDate - valueDate) / gDaysInYear

        K = option._strikePrice
        dt = 1.0 / numStepsPerYear
        schemeValue = float(scheme.value)

        sPaths = getPaths(stockPrice,
                          interestRate,
                          dividendYield,
                          self._v0,
                          self._kappa,
                          self._theta,
                          self._sigma,
                          self._rho,
                          tau,
                          dt,
                          numPaths,
                          seed,
                          schemeValue)

        if option._optionType == FinOptionTypes.EUROPEAN_CALL:
            path_payoff = np.maximum(sPaths[:, -1] - K, 0.0)
        elif option._optionType == FinOptionTypes.EUROPEAN_PUT:
            path_payoff = np.maximum(K - sPaths[:, -1], 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(path_payoff)
        v = payoff * exp(-interestRate * tau)
        return v

###############################################################################

    def value_Lewis(self,
                    valueDate,
                    option,
                    stockPrice,
                    interestRate,
                    dividendYield):

        tau = (option._expiryDate - valueDate) / gDaysInYear

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        r = interestRate
        q = dividendYield
        S0 = stockPrice
        K = option._strikePrice
        F = S0 * exp((r - q) * tau)
        V = sigma * sigma

        def phi(k_in,):
            k = k_in + 0.5 * 1j
            b = kappa + 1j * rho * sigma * k
            d = np.sqrt(b**2 + V * k * (k - 1j))
            g = (b - d) / (b + d)
            T_m = (b - d) / V
            Q = np.exp(-d * tau)
            T = T_m * (1.0 - Q) / (1.0 - g * Q)
            W = kappa * theta * (tau * T_m - 2.0 *
                                 np.log((1.0 - g * Q) / (1.0 - g)) / V)
            phi = np.exp(W + v0 * T)
            return phi

        def phi_transform(x):
            def integrand(k): return 2.0 * np.real(np.exp(-1j * \
                          k * x) * phi(k)) / (k**2 + 1.0 / 4.0)
            return integrate.quad(integrand, 0, np.inf)[0]

        x = log(F / K)
        I1 = phi_transform(x) / (2.0 * pi)
        v1 = F * exp(-r * tau) - np.sqrt(K * F) * exp(-r * tau) * I1
#        v2 = S0 * exp(-q*tau) - K * exp(-r*tau) * I1
        return(v1)

###############################################################################

    def value_Lewis_Rouah(self,
                          valueDate,
                          option,
                          stockPrice,
                          interestRate,
                          dividendYield):

        tau = (option._expiryDate - valueDate) / gDaysInYear

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividendYield
        r = interestRate
        V = sigma * sigma

        def f(k_in):
            k = k_in + 0.5 * 1j
            b = (2.0 / V) * (1j * k * rho * sigma + kappa)
            e = np.sqrt(b**2 + 4.0 * k * (k - 1j) / V)
            g = (b - e) / 2.0
            h = (b - e) / (b + e)
            q = V * tau / 2.0
            Q = np.exp(-e * q)
            H = np.exp((2.0 * kappa * theta / V) * (q * g - np.log((1.0 -
                  h * Q) / (1.0 - h))) + v0 * g * (1.0 - Q) / (1.0 - h * Q))
            integrand = H * np.exp(-1j * k * X) / (k * k - 1j * k)
            return integrand.real

        S0 = stockPrice
        F = S0 * exp((r - q) * tau)
        K = option._strikePrice
        X = log(F / K)
        integral = integrate.quad(f, 0.0, np.inf)[0] * (1.0 / pi)
        v = S0 * exp(-q * tau) - K * exp(-r * tau) * integral
        return (v)

###############################################################################
# Taken from Nick Weber's VBA Finance book
###############################################################################

    def value_Weber(self,
                    valueDate,
                    option,
                    stockPrice,
                    interestRate,
                    dividendYield):

        tau = (option._expiryDate - valueDate) / gDaysInYear

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividendYield
        r = interestRate
        S0 = stockPrice
        K = option._strikePrice
        V = sigma**2

        def F(s, b):
            def integrand(u):
                beta = b - 1j * rho * sigma * u
                d = np.sqrt((beta**2) - V * u * (s * 1j - u))
                g = (beta - d) / (beta + d)
                Q = np.exp(-d * tau)
                B = (beta - d) * (1.0 - Q) / (1.0 - g * Q) / V
                A = kappa * ((beta - d) * tau - 2.0 *
                             np.log((1.0 - g * Q) / (1.0 - g))) / V
                v = np.exp(A * theta + B * v0 + 1j * u * \
                           np.log(S0 / (K * np.exp(-(r - q) * tau)))) / (u * 1j)
                return v.real

            area = 0.50 + (1.0 / pi) * integrate.quad(integrand, 0, np.inf)[0]
            return area

        v = S0 * exp(-q * tau) * F(1.0, kappa - rho * sigma) - \
            exp(-r * tau) * K * F(-1.0, kappa)
        return(v)

###############################################################################
# Gatheral book page 19 with definition of x given on page 16 and noting
# that the value C is a forward value and so needs to be discounted
###############################################################################

    def value_Gatheral(self,
                       valueDate,
                       option,
                       stockPrice,
                       interestRate,
                       dividendYield):

        tau = (option._expiryDate - valueDate) / gDaysInYear

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividendYield
        r = interestRate
        S0 = stockPrice
        K = option._strikePrice
        F = S0 * exp((r - q) * tau)
        x0 = log(F / K)

        def FF(j):
            def integrand(u):
                V = sigma * sigma
                A = -u * u / 2.0 - 1j * u / 2.0 + 1j * j * u
                B = kappa - rho * sigma * j - rho * sigma * 1j * u
                G = V / 2.0
                d = np.sqrt(B**2 - 4.0 * A * G)
                rplus = (B + d) / 2.0 / G
                rminus = (B - d) / 2.0 / G
                R = rminus / rplus
                Q = np.exp(-d * tau)
                D = rminus * (1.0 - Q) / (1.0 - R * Q)
                C = kappa * (rminus * tau - (2.0 / V) *
                             np.log((1.0 - R * Q) / (1.0 - R)))
                phi = np.exp(C * theta + D * v0 + 1j * u * x0) / (1j * u)
                return phi.real

            area = 0.50 + 1.0 / pi * integrate.quad(integrand, 0.0, np.inf)[0]
            return area

        v = S0 * exp(-q * tau) * FF(1) - K * exp(-r * tau) * FF(0)
        return(v)

###############################################################################
