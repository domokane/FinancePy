# Copyright (c) 2018, 2019, 2020 Dominic O'Kane

from numba import njit, float64, int64
from scipy import integrate
from math import exp, log, pi
import numpy as np

from ..utils.global_vars import G_DAYS_IN_YEARS
from ..utils.global_types import OptionTypes
from ..utils.math import norminvcdf
from ..utils.error import FinError

# Heston Process
# dS = rS dt + sqrt(v) * S * dz
# dV = kappa(theta-v) dt + sigma sqrt(v) dz
# corr(dV,dS) = rho dt
# Rewritten as
# dS = rS dt + sqrt(v) * S * (rhohat dz1 + rho dz2)
# dV = kappa(theta-v) dt + sigma sqrt(v) dz2
# where rhohat = sqrt(1-rho*rho)
# TODO - DECIDE WHETHER TO OO MODEL
# TODO - NEEDS CHECKING FOR MC CONVERGENCE

from enum import Enum

########################################################################################

class HestonNumericalScheme(Enum):

    euler = 1
    EULERLOG = 2
    QUADEXP = 3

########################################################################################

@njit(
    float64[:, :](
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64,
    ),
    cache=True,
    fastmath=True,
)
def get_paths(

    s0, r, q, v0, kappa, theta, sigma, rho, t, dt, num_paths, seed, scheme
):

    np.random.seed(seed)
    num_steps = int(t / dt)
    s_paths = np.zeros(shape=(num_paths, num_steps))
    s_paths[:, 0] = s0
    sdt = np.sqrt(dt)
    rhohat = np.sqrt(1.0 - rho * rho)
    sigma2 = sigma * sigma

    if scheme == HestonNumericalScheme.EULER.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            s = s0
            v = v0
            for i_step in range(1, num_steps):
                z1 = np.random.normal(0.0, 1.0) * sdt
                z2 = np.random.normal(0.0, 1.0) * sdt
                z_v = z1
                z_s = rho * z1 + rhohat * z2
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                v += (
                    kappa * (theta - vplus) * dt
                    + sigma * rtvplus * z_v
                    + 0.25 * sigma2 * (z_v * z_v - dt)
                )
                s += (
                    (r - q) * s * dt
                    + rtvplus * s * z_s
                    + 0.5 * s * vplus * (z_v * z_v - dt)
                )
                s_paths[i_path, i_step] = s

    elif scheme == HestonNumericalScheme.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            x = log(s0)
            v = v0
            for i_step in range(1, num_steps):
                z_v = np.random.normal(0.0, 1.0) * sdt
                z_s = rho * z_v + rhohat * np.random.normal(0.0, 1.0) * sdt
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                x += (r - q - 0.5 * vplus) * dt + rtvplus * z_s
                v += (
                    kappa * (theta - vplus) * dt
                    + sigma * rtvplus * z_v
                    + sigma2 * (z_v * z_v - dt) / 4.0
                )
                s_paths[i_path, i_step] = exp(x)

    elif scheme == HestonNumericalScheme.QUADEXP.value:
        # Due to Leif Andersen(2006)
        q = exp(-kappa * dt)
        psic = 1.50
        gamma1 = 0.50
        gamma2 = 0.50
        k0 = -rho * kappa * theta * dt / sigma
        k_1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho / sigma
        k_2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho / sigma
        k3 = gamma1 * dt * (1.0 - rho * rho)
        k4 = gamma2 * dt * (1.0 - rho * rho)
        a = k_2 + 0.5 * k4
        mu = r - q
        c1 = sigma2 * q * (1.0 - q) / kappa
        c2 = theta * sigma2 * ((1.0 - q) ** 2) / 2.0 / kappa

        for i_path in range(0, num_paths):
            x = log(s0)
            vn = v0
            for i_step in range(1, num_steps):
                z_v = np.random.normal(0, 1)
                z_s = rho * z_v + rhohat * np.random.normal(0, 1)
                m = theta + (vn - theta) * q
                m2 = m * m
                s2 = c1 * vn + c2
                psi = s2 / m2
                u = np.random.uniform(0.0, 1.0)

                if psi <= psic:
                    b2 = (
                        2.0 / psi
                        - 1.0
                        + np.sqrt((2.0 / psi) * (2.0 / psi - 1.0))
                    )
                    a = m / (1.0 + b2)
                    b = np.sqrt(b2)
                    z_v = norminvcdf(u)
                    vnp = a * ((b + z_v) ** 2)
                    d = 1.0 - 2.0 * a * a
                    m = exp((a * b2 * a) / d) / np.sqrt(d)
                    k0 = -log(m) - (k_1 + 0.5 * k3) * vn
                else:
                    p = (psi - 1.0) / (psi + 1.0)
                    beta = (1.0 - p) / m

                    if u <= p:
                        vnp = 0.0
                    else:
                        vnp = log((1.0 - p) / (1.0 - u)) / beta

                    m = p + beta * (1.0 - p) / (beta - a)
                    k0 = -log(m) - (k_1 + 0.5 * k3) * vn

                x += (
                    mu * dt
                    + k0
                    + (k_1 * vn + k_2 * vnp)
                    + np.sqrt(k3 * vn + k4 * vnp) * z_s
                )
                s_paths[i_path, i_step] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return s_paths

########################################################################################

class Heston:

    ####################################################################################

    def __init__(self, v0, kappa, theta, sigma, rho):

        verbose = False

        if 2.0 * kappa * theta <= sigma and verbose:
            print("Feller condition not satisfied. Zero Variance possible")

        self._v0 = v0
        self._kappa = kappa
        self._theta = theta
        self._sigma = sigma
        self._rho = rho

    ####################################################################################

    def value_mc(

        self,
        value_dt,
        option,
        stock_price,
        interest_rate,
        dividend_yield,
        num_paths,
        num_steps_per_year,
        seed,
        scheme=HestonNumericalScheme.EULERLOG,
    ):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        k = option.strike_price
        dt = 1.0 / num_steps_per_year
        scheme_value = float(scheme.value)

        s_paths = get_paths(
            stock_price,
            interest_rate,
            dividend_yield,
            self._v0,
            self._kappa,
            self._theta,
            self._sigma,
            self._rho,
            tau,
            dt,
            num_paths,
            seed,
            scheme_value,
        )

        if option.opt_type == OptionTypes.EUROPEAN_CALL:
            path_payoff = np.maximum(s_paths[:, -1] - k, 0.0)
        elif option.opt_type == OptionTypes.EUROPEAN_PUT:
            path_payoff = np.maximum(k - s_paths[:, -1], 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(path_payoff)
        v = payoff * exp(-interest_rate * tau)
        return v

    ####################################################################################

    def value_lewis(

        self, value_dt, option, stock_price, interest_rate, dividend_yield
    ):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        r = interest_rate
        q = dividend_yield
        s0 = stock_price
        k = option.strike_price
        f = s0 * exp((r - q) * tau)
        v = sigma * sigma

        ################################################################################

        def phi(

            k_in,
        ):
            k = k_in + 0.5 * 1j
            b = kappa + 1j * rho * sigma * k
            d = np.sqrt(b**2 + v * k * (k - 1j))
            g = (b - d) / (b + d)
            t_m = (b - d) / v
            q = np.exp(-d * tau)
            t = t_m * (1.0 - q) / (1.0 - g * q)
            w = (
                kappa
                * theta
                * (tau * t_m - 2.0 * np.log((1.0 - g * q) / (1.0 - g)) / v)
            )
            phi = np.exp(w + v0 * t)
            return phi

        ################################################################################

        def phi_transform(x):

            ############################################################################

            def integrand(k):

                return (
                    2.0
                    * np.real(np.exp(-1j * k * x) * phi(k))
                    / (k**2 + 1.0 / 4.0)
                )

            return integrate.quad(integrand, 0, np.inf)[0]

        x = log(f / k)
        i1 = phi_transform(x) / (2.0 * pi)
        v1 = f * exp(-r * tau) - np.sqrt(k * f) * exp(-r * tau) * i1
        #        v2 = s0 * exp(-q*tau) - k * exp(-r*tau) * i1
        return v1

    ####################################################################################

    def value_lewis_rouah(

        self, value_dt, option, stock_price, interest_rate, dividend_yield
    ):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        v = sigma * sigma

        ################################################################################

        def f(k_in):

            k = k_in + 0.5 * 1j
            b = (2.0 / v) * (1j * k * rho * sigma + kappa)
            e = np.sqrt(b**2 + 4.0 * k * (k - 1j) / v)
            g = (b - e) / 2.0
            h = (b - e) / (b + e)
            q = v * tau / 2.0
            q = np.exp(-e * q)
            h = np.exp(
                (2.0 * kappa * theta / v)
                * (q * g - np.log((1.0 - h * q) / (1.0 - h)))
                + v0 * g * (1.0 - q) / (1.0 - h * q)
            )
            integrand = h * np.exp(-1j * k * x) / (k * k - 1j * k)
            return integrand.real

        s0 = stock_price
        f = s0 * exp((r - q) * tau)
        k = option.strike_price
        x = log(f / k)
        integral = integrate.quad(f, 0.0, np.inf)[0] * (1.0 / pi)
        v = s0 * exp(-q * tau) - k * exp(-r * tau) * integral
        return v

    # Taken from Nick Weber's VBA Finance book

    ####################################################################################

    def value_weber(

        self, value_dt, option, stock_price, interest_rate, dividend_yield
    ):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        s0 = stock_price
        k = option.strike_price
        v = sigma**2

        ################################################################################

        def f(s, b):

            ############################################################################

            def integrand(u):

                beta = b - 1j * rho * sigma * u
                d = np.sqrt((beta**2) - v * u * (s * 1j - u))
                g = (beta - d) / (beta + d)
                q = np.exp(-d * tau)
                b = (beta - d) * (1.0 - q) / (1.0 - g * q) / v
                a = (
                    kappa
                    * (
                        (beta - d) * tau
                        - 2.0 * np.log((1.0 - g * q) / (1.0 - g))
                    )
                    / v
                )
                v = np.exp(
                    a * theta
                    + b * v0
                    + 1j * u * np.log(s0 / (k * np.exp(-(r - q) * tau)))
                ) / (u * 1j)
                return v.real

            area = 0.50 + (1.0 / pi) * integrate.quad(integrand, 0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * f(1.0, kappa - rho * sigma) - exp(
            -r * tau
        ) * k * f(-1.0, kappa)

        return v

    # Gatheral book page 19 with definition of x given on page 16 and noting
    # that the value c is a forward value and so needs to be discounted

    ####################################################################################

    def value_gatheral(

        self, value_dt, option, stock_price, interest_rate, dividend_yield
    ):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        s0 = stock_price
        k = option.strike_price
        f = s0 * exp((r - q) * tau)
        x0 = log(f / k)

        ################################################################################

        def ff(j):

            ############################################################################

            def integrand(u):

                v = sigma * sigma
                a = -u * u / 2.0 - 1j * u / 2.0 + 1j * j * u
                b = kappa - rho * sigma * j - rho * sigma * 1j * u
                g = v / 2.0
                d = np.sqrt(b**2 - 4.0 * a * g)
                rplus = (b + d) / 2.0 / g
                rminus = (b - d) / 2.0 / g
                r = rminus / rplus
                q = np.exp(-d * tau)
                d = rminus * (1.0 - q) / (1.0 - r * q)
                c = kappa * (
                    rminus * tau
                    - (2.0 / v) * np.log((1.0 - r * q) / (1.0 - r))
                )
                phi = np.exp(c * theta + d * v0 + 1j * u * x0) / (1j * u)
                return phi.real

            area = 0.50 + 1.0 / pi * integrate.quad(integrand, 0.0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * ff(1) - k * exp(-r * tau) * ff(0)
        return v

