##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
from math import exp, log, pi

from numba import njit, float64, int64
from scipy import integrate
import numpy as np

from ..utils.global_vars import G_DAYS_IN_YEARS
from ..utils.global_types import OptionTypes
from ..utils.math import norminvcdf
from ..utils.error import FinError

##########################################################################
# Heston Process
# dS = rS dt + sqrt(V) * S * dz
# dV = kappa(theta-V) dt + sigma sqrt(V) dz
# corr(dV,dS) = rho dt
# Rewritten as
# dS = rS dt + sqrt(V) * S * (rhohat dz1 + rho dz2)
# dV = kappa(theta-V) dt + sigma sqrt(V) dz2
# where rhohat = sqrt(1-rho*rho)
########################################################################################
# TODO - DECIDE WHETHER TO OO MODEL
# TODO - NEEDS CHECKING FOR MC CONVERGENCE
########################################################################################


class HestonNumericalScheme(Enum):
    EULER = 1
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
    parallel=False,
)
def get_paths(s0, r, q, v0, kappa, theta, sigma, rho, t, dt, num_paths, seed, scheme):

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
        qq = exp(-kappa * dt)
        psic = 1.50
        gamma1 = 0.50
        gamma2 = 0.50
        k_0 = -rho * kappa * theta * dt / sigma
        k_1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho / sigma
        k_2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho / sigma
        k_3 = gamma1 * dt * (1.0 - rho * rho)
        k_4 = gamma2 * dt * (1.0 - rho * rho)
        aa = k_2 + 0.5 * k_4
        mu = r - q
        c1 = sigma2 * qq * (1.0 - qq) / kappa
        c2 = theta * sigma2 * ((1.0 - qq) ** 2) / 2.0 / kappa

        for i_path in range(0, num_paths):
            x = log(s0)
            vn = v0
            for i_step in range(1, num_steps):
                z_v = np.random.normal(0, 1)
                z_s = rho * z_v + rhohat * np.random.normal(0, 1)
                m = theta + (vn - theta) * qq
                m2 = m * m
                s2 = c1 * vn + c2
                psi = s2 / m2
                u = np.random.uniform(0.0, 1.0)

                if psi <= psic:
                    b2 = 2.0 / psi - 1.0 + np.sqrt((2.0 / psi) * (2.0 / psi - 1.0))
                    a = m / (1.0 + b2)
                    b = np.sqrt(b2)
                    z_v = norminvcdf(u)
                    vnp = a * ((b + z_v) ** 2)
                    d = 1.0 - 2.0 * aa * a
                    m = exp((aa * b2 * a) / d) / np.sqrt(d)
                    k_0 = -log(m) - (k_1 + 0.5 * k_3) * vn
                else:
                    p = (psi - 1.0) / (psi + 1.0)
                    beta = (1.0 - p) / m

                    if u <= p:
                        vnp = 0.0
                    else:
                        vnp = log((1.0 - p) / (1.0 - u)) / beta

                    m = p + beta * (1.0 - p) / (beta - aa)
                    k_0 = -log(m) - (k_1 + 0.5 * k_3) * vn

                x += (
                    mu * dt
                    + k_0
                    + (k_1 * vn + k_2 * vnp)
                    + np.sqrt(k_3 * vn + k_4 * vnp) * z_s
                )
                s_paths[i_path, i_step] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return s_paths


########################################################################################


class Heston:

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

    def value_lewis(self, value_dt, option, stock_price, interest_rate, dividend_yield):

        tau = (option.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        rho = self._rho
        sigma = self._sigma
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        r = interest_rate
        q = dividend_yield
        s0 = stock_price
        kk = option.strike_price
        ff = s0 * exp((r - q) * tau)
        vv = sigma * sigma

        def phi(
            k_in,
        ):
            k = k_in + 0.5 * 1j
            b = kappa + 1j * rho * sigma * k
            d = np.sqrt(b**2 + vv * k * (k - 1j))
            g = (b - d) / (b + d)
            t_m = (b - d) / vv
            qq = np.exp(-d * tau)
            t = t_m * (1.0 - qq) / (1.0 - g * qq)
            ww = (
                kappa
                * theta
                * (tau * t_m - 2.0 * np.log((1.0 - g * qq) / (1.0 - g)) / vv)
            )
            phi = np.exp(ww + v0 * t)
            return phi

        def phi_transform(x):
            def integrand(k):
                return 2.0 * np.real(np.exp(-1j * k * x) * phi(k)) / (k**2 + 1.0 / 4.0)

            return integrate.quad(integrand, 0, np.inf)[0]

        x = log(ff / kk)
        i_1 = phi_transform(x) / (2.0 * pi)
        v1 = ff * exp(-r * tau) - np.sqrt(kk * ff) * exp(-r * tau) * i_1
        #        v2 = s0 * exp(-q*tau) - K * exp(-r*tau) * I1
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
        vv = sigma * sigma

        s0 = stock_price
        f = s0 * exp((r - q) * tau)
        k = option.strike_price
        x = log(f / k)

        def fn(k_in):
            k = k_in + 0.5 * 1j
            b = (2.0 / vv) * (1j * k * rho * sigma + kappa)
            e = np.sqrt(b**2 + 4.0 * k * (k - 1j) / vv)
            g = (b - e) / 2.0
            h = (b - e) / (b + e)
            q = vv * tau / 2.0
            qq = np.exp(-e * q)
            hh = np.exp(
                (2.0 * kappa * theta / vv)
                * (q * g - np.log((1.0 - h * qq) / (1.0 - h)))
                + v0 * g * (1.0 - qq) / (1.0 - h * qq)
            )
            integrand = hh * np.exp(-1j * k * x) / (k * k - 1j * k)
            return integrand.real

        integral = integrate.quad(fn, 0.0, np.inf)[0] * (1.0 / pi)
        v = s0 * exp(-q * tau) - k * exp(-r * tau) * integral
        return v

    ####################################################################################
    # Taken from Nick Weber's VBA Finance book
    ####################################################################################

    def value_weber(self, value_dt, option, stock_price, interest_rate, dividend_yield):

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
        vv = sigma**2

        def fn(s, b):
            def integrand(u):
                beta = b - 1j * rho * sigma * u
                d = np.sqrt((beta**2) - vv * u * (s * 1j - u))
                g = (beta - d) / (beta + d)
                qq = np.exp(-d * tau)
                bb = (beta - d) * (1.0 - qq) / (1.0 - g * qq) / vv
                aa = (
                    kappa
                    * ((beta - d) * tau - 2.0 * np.log((1.0 - g * qq) / (1.0 - g)))
                    / vv
                )
                v = np.exp(
                    aa * theta
                    + bb * v0
                    + 1j * u * np.log(s0 / (k * np.exp(-(r - q) * tau)))
                ) / (u * 1j)
                return v.real

            area = 0.50 + (1.0 / pi) * integrate.quad(integrand, 0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * fn(1.0, kappa - rho * sigma)
        v = v - exp(-r * tau) * k * fn(-1.0, kappa)

        return v

    ####################################@#@#############################################
    # Gatheral book page 19 with definition of x given on page 16 and noting
    # that the value C is a forward value and so needs to be discounted
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

        def ff(j):
            def integrand(u):
                vv = sigma * sigma
                aa = -u * u / 2.0 - 1j * u / 2.0 + 1j * j * u
                bb = kappa - rho * sigma * j - rho * sigma * 1j * u
                gg = vv / 2.0
                d = np.sqrt(bb**2 - 4.0 * aa * gg)
                rplus = (bb + d) / 2.0 / gg
                rminus = (bb - d) / 2.0 / gg
                rr = rminus / rplus
                qq = np.exp(-d * tau)
                dd = rminus * (1.0 - qq) / (1.0 - rr * qq)
                cc = kappa * (
                    rminus * tau - (2.0 / vv) * np.log((1.0 - rr * qq) / (1.0 - rr))
                )
                phi = np.exp(cc * theta + dd * v0 + 1j * u * x0) / (1j * u)
                return phi.real

            area = 0.50 + 1.0 / pi * integrate.quad(integrand, 0.0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * ff(1) - k * exp(-r * tau) * ff(0)
        return v


########################################################################################
