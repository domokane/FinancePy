# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from math import sqrt, exp, log
from enum import Enum

from numba import njit, float64, int64
import numpy as np

from ..utils.error import FinError
from ..utils.math import norminvcdf

########################################################################################

class ProcessTypes(Enum):

    gbm = 1
    cir = 2
    heston = 3
    VASICEK = 4
    cev = 5
    JUMP_DIFFUSION = 6

########################################################################################

class FinProcessSimulator:

    ####################################################################################

    def __init__(self):

        pass

    ####################################################################################

    def get_process(

        self, process_type, t, model_params, num_annual_steps, num_paths, seed
    ):

        if process_type == ProcessTypes.GBM:
            (stock_price, drift, volatility, scheme) = model_params
            paths = get_gbm_paths(
                num_paths,
                num_annual_steps,
                t,
                drift,
                stock_price,
                volatility,
                scheme.value,
                seed,
            )
            return paths

        elif process_type == ProcessTypes.HESTON:

            (stock_price, drift, v0, kappa, theta, sigma, rho, scheme) = model_params
            paths = get_heston_paths(
                num_paths,
                num_annual_steps,
                t,
                drift,
                stock_price,
                v0,
                kappa,
                theta,
                sigma,
                rho,
                scheme.value,
                seed,
            )
            return paths

        elif process_type == ProcessTypes.VASICEK:

            (r0, kappa, theta, sigma, scheme) = model_params
            paths = get_vasicek_paths(
                num_paths,
                num_annual_steps,
                t,
                r0,
                kappa,
                theta,
                sigma,
                scheme.value,
                seed,
            )
            return paths

        elif process_type == ProcessTypes.CIR:
            (r0, kappa, theta, sigma, scheme) = model_params
            paths = get_cir_paths(
                num_paths,
                num_annual_steps,
                t,
                r0,
                kappa,
                theta,
                sigma,
                scheme.value,
                seed,
            )
            return paths

        else:
            raise FinError("Unknown process" + str(process_type))

########################################################################################

class FinHestonNumericalScheme(Enum):

    euler = 1
    EULERLOG = 2
    QUADEXP = 3

########################################################################################

@njit(
    float64[:, :](
        int64,
        int64,
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
    ),
    cache=True,
    fastmath=True,
)
def get_heston_paths(

    num_paths,
    num_annual_steps,
    t,
    drift,
    s0,
    v0,
    kappa,
    theta,
    sigma,
    rho,
    scheme,
    seed,
):

    np.random.seed(seed)
    dt = 1.0 / num_annual_steps
    num_steps = int(t / dt)
    s_paths = np.empty(shape=(num_paths, num_steps + 1))
    s_paths[:, 0] = s0
    sdt = sqrt(dt)
    rhohat = sqrt(1.0 - rho * rho)
    sigma2 = sigma * sigma

    if scheme == FinHestonNumericalScheme.EULER.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            s = s0
            v = v0
            for i_step in range(1, num_steps + 1):
                z1 = np.random.normal(0.0, 1.0) * sdt
                z2 = np.random.normal(0.0, 1.0) * sdt
                z_v = z1
                z_s = rho * z1 + rhohat * z2
                v_plus = max(v, 0.0)
                rtv_plus = sqrt(v_plus)
                v += (
                    kappa * (theta - v_plus) * dt
                    + sigma * rtv_plus * z_v
                    + 0.25 * sigma2 * (z_v * z_v - dt)
                )
                s += (
                    drift * s * dt
                    + rtv_plus * s * z_s
                    + 0.5 * s * v_plus * (z_v * z_v - dt)
                )
                s_paths[i_path, i_step] = s

    elif scheme == FinHestonNumericalScheme.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            x = log(s0)
            v = v0
            for i_step in range(1, num_steps + 1):
                z_v = np.random.normal(0.0, 1.0) * sdt
                z_s = rho * z_v + rhohat * np.random.normal(0.0, 1.0) * sdt
                v_plus = max(v, 0.0)
                rtv_plus = sqrt(v_plus)
                x += (drift - 0.5 * v_plus) * dt + rtv_plus * z_s
                v += (
                    kappa * (theta - v_plus) * dt
                    + sigma * rtv_plus * z_v
                    + sigma2 * (z_v * z_v - dt) / 4.0
                )
                s_paths[i_path, i_step] = exp(x)

    elif scheme == FinHestonNumericalScheme.QUADEXP.value:
        # Due to Leif Andersen(2006)
        q = exp(-kappa * dt)
        psic = 1.50
        gamma1 = 0.50
        gamma2 = 0.50
        k_0 = -rho * kappa * theta * dt / sigma
        k_1 = gamma1 * dt * (kappa * rho / sigma - 0.5) - rho / sigma
        k_2 = gamma2 * dt * (kappa * rho / sigma - 0.5) + rho / sigma
        k_3 = gamma1 * dt * (1.0 - rho * rho)
        k4 = gamma2 * dt * (1.0 - rho * rho)
        a = k_2 + 0.5 * k4
        mu = drift
        c1 = sigma2 * q * (1.0 - q) / kappa
        c2 = theta * sigma2 * ((1.0 - q) ** 2) / 2.0 / kappa

        for i_path in range(0, num_paths):
            x = log(s0)
            vn = v0
            for i_step in range(1, num_steps + 1):
                z_v = np.random.normal(0, 1)
                z_s = rho * z_v + rhohat * np.random.normal(0, 1)
                m = theta + (vn - theta) * q
                m2 = m * m
                s2 = c1 * vn + c2
                psi = s2 / m2
                u = np.random.uniform(0.0, 1.0)

                if psi <= psic:
                    b2 = 2.0 / psi - 1.0 + sqrt((2.0 / psi) * (2.0 / psi - 1.0))
                    a = m / (1.0 + b2)
                    b = sqrt(b2)
                    z_v = norminvcdf(u)
                    vnp = a * ((b + z_v) ** 2)
                    d = 1.0 - 2.0 * a * a
                    m = exp((a * b2 * a) / d) / sqrt(d)
                    k_0 = -log(m) - (k_1 + 0.5 * k_3) * vn
                else:
                    p = (psi - 1.0) / (psi + 1.0)
                    beta = (1.0 - p) / m

                    if u <= p:
                        vnp = 0.0
                    else:
                        vnp = log((1.0 - p) / (1.0 - u)) / beta

                    m = p + beta * (1.0 - p) / (beta - a)
                    k_0 = -log(m) - (k_1 + 0.5 * k_3) * vn

                x += (
                    mu * dt
                    + k_0
                    + (k_1 * vn + k_2 * vnp)
                    + sqrt(k_3 * vn + k4 * vnp) * z_s
                )
                s_paths[i_path, i_step] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return s_paths

########################################################################################

class FinGBMNumericalScheme(Enum):

    normal = 1
    ANTITHETIC = 2

########################################################################################

@njit(
    float64[:, :](int64, int64, float64, float64, float64, float64, int64, int64),
    cache=True,
    fastmath=True,
)
def get_gbm_paths(num_paths, num_annual_steps, t, mu, stock_price, sigma, scheme, seed):

    np.random.seed(seed)
    dt = 1.0 / num_annual_steps
    num_time_steps = int(t / dt + 0.50)
    vsqrt_dt = sigma * sqrt(dt)
    m = exp((mu - sigma * sigma / 2.0) * dt)

    if scheme == FinGBMNumericalScheme.NORMAL.value:

        s_all = np.empty((num_paths, num_time_steps + 1))
        s_all[:, 0] = stock_price
        for it in range(1, num_time_steps + 1):
            g1_d = np.random.standard_normal((num_paths))
            for ip in range(0, num_paths):
                w = np.exp(g1_d[ip] * vsqrt_dt)
                s_all[ip, it] = s_all[ip, it - 1] * m * w

    elif scheme == FinGBMNumericalScheme.ANTITHETIC.value:

        s_all = np.empty((2 * num_paths, num_time_steps + 1))
        s_all[:, 0] = stock_price
        for it in range(1, num_time_steps + 1):
            g1_d = np.random.standard_normal((num_paths))
            for ip in range(0, num_paths):
                w = np.exp(g1_d[ip] * vsqrt_dt)
                s_all[ip, it] = s_all[ip, it - 1] * m * w
                s_all[ip + num_paths, it] = s_all[ip + num_paths, it - 1] * m / w

    else:

        raise FinError("Unknown FinGBMNumericalScheme")

    #    m = np.mean(s_all[:, -1])
    #    v = np.var(s_all[:, -1]/s_all[:, 0])
    #    print("gbm", num_paths, num_annual_steps, t, mu, stock_price, sigma, scheme, m,v)

    return s_all

########################################################################################

class FinVasicekNumericalScheme(Enum):

    normal = 1
    ANTITHETIC = 2

########################################################################################

@njit(
    float64[:, :](
        int64, int64, float64, float64, float64, float64, float64, int64, int64
    ),
    cache=True,
    fastmath=True,
)
def get_vasicek_paths(

    num_paths, num_annual_steps, t, r0, kappa, theta, sigma, scheme, seed
):

    np.random.seed(seed)
    dt = 1.0 / num_annual_steps
    num_steps = int(t / dt)
    sigma_sqrt_dt = sigma * sqrt(dt)

    if scheme == FinVasicekNumericalScheme.NORMAL.value:
        rate_path = np.empty((num_paths, num_steps + 1))
        rate_path[:, 0] = r0
        for i_path in range(0, num_paths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                r += kappa * (theta - r) * dt + z[i_step - 1] * sigma_sqrt_dt
                rate_path[i_path, i_step] = r
    elif scheme == FinVasicekNumericalScheme.ANTITHETIC.value:
        rate_path = np.empty((2 * num_paths, num_steps + 1))
        rate_path[:, 0] = r0
        for i_path in range(0, num_paths):
            r1 = r0
            r2 = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                r1 = r1 + kappa * (theta - r1) * dt + z[i_step - 1] * sigma_sqrt_dt
                r2 = r2 + kappa * (theta - r2) * dt - z[i_step - 1] * sigma_sqrt_dt
                rate_path[i_path, i_step] = r1
                rate_path[i_path + num_paths, i_step] = r2
    return rate_path

########################################################################################

class CIRNumericalScheme(Enum):

    euler = 1
    LOGNORMAL = 2
    MILSTEIN = 3
    KAHLJACKEL = 4
    exact = 5  # SAMPLES exact DISTRIBUTION

########################################################################################

@njit(
    float64[:, :](
        int64, int64, float64, float64, float64, float64, float64, int64, int64
    ),
    cache=True,
    fastmath=True,
)
def get_cir_paths(

    num_paths, num_annual_steps, t, r0, kappa, theta, sigma, scheme, seed
########################################################################################

):

    np.random.seed(seed)
    dt = 1.0 / num_annual_steps
    num_steps = int(t / dt)
    rate_path = np.empty(shape=(num_paths, num_steps + 1))
    rate_path[:, 0] = r0

    if scheme == CIRNumericalScheme.EULER.value:
        sigma_sqrt_dt = sigma * sqrt(dt)
        for i_path in range(0, num_paths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                rplus = max(r, 0.0)
                sqrt_rplus = sqrt(rplus)
                r = (
                    r
                    + kappa * (theta - rplus) * dt
                    + sigma_sqrt_dt * z[i_step - 1] * sqrt_rplus
                )
                rate_path[i_path, i_step] = r

    elif scheme == CIRNumericalScheme.LOGNORMAL.value:
        x = exp(-kappa * dt)
        y = 1.0 - x
        for i_path in range(0, num_paths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                mean = x * r + theta * y
                var = sigma * sigma * y * (x * r + 0.50 * theta * y) / kappa
                sig = sqrt(log(1.0 + var / (mean * mean)))
                r = mean * exp(-0.5 * sig * sig + sig * z[i_step - 1])
                rate_path[i_path, i_step] = r

    elif scheme == CIRNumericalScheme.MILSTEIN.value:
        sigma_sqrt_dt = sigma * sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0
        for i_path in range(0, num_paths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                sqrt_rplus = sqrt(max(r, 0.0))
                r = (
                    r
                    + kappa * (theta - r) * dt
                    + z[i_step - 1] * sigma_sqrt_dt * sqrt_rplus
                )
                r = r + sigma2dt * (z[i_step - 1] ** 2 - 1.0)
                rate_path[i_path, i_step] = r

    elif scheme == CIRNumericalScheme.KAHLJACKEL.value:
        bhat = theta - sigma * sigma / 4.0 / kappa
        sqrt_dt = sqrt(dt)
        for i_path in range(0, num_paths):
            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps))
            for i_step in range(1, num_steps + 1):
                beta = z[i_step - 1] / sqrt_dt
                sqrt_rplus = sqrt(max(r, 0.0))
                c = (
                    1.0
                    + (sigma * beta - 2.0 * kappa * sqrt_rplus) * dt / 4.0 / sqrt_rplus
                )
                r = r + (kappa * (bhat - r) + sigma * beta * sqrt_rplus) * c * dt
                rate_path[i_path, i_step] = r

    return rate_path

