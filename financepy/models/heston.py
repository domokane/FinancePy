##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, pi

from numba import njit, float64, int64
from scipy import integrate
import numpy as np

from ..utils.global_types import OptionTypes
from ..models.process_simulator import HestonNumericalSchemeTypes
from ..utils.math import norminvcdf
from ..utils.error import FinError
from ..models.black_scholes_analytic import implied_volatility


from enum import Enum


class HestonValueTypes(Enum):
    LEWIS = 0
    LEWIS_ROUAH = 1
    GATHERAL = 2
    WEBER = 3


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
def get_paths(
    s0: float,
    r: float,
    q: float,
    v0: float,
    kappa: float,
    theta: float,
    sigma: float,
    rho: float,
    t: float,
    dt: float,
    num_paths: int,
    seed: int,
    scheme: float,
) -> np.ndarray:

    np.random.seed(seed)

    num_steps = max(1, int(round(t / dt)))
    dt = t / num_steps

    s_paths = np.zeros((num_paths, num_steps + 1))
    s_paths[:, 0] = s0

    sdt = np.sqrt(dt)
    rhohat = np.sqrt(1.0 - rho * rho)
    sigma2 = sigma * sigma

    if scheme == HestonNumericalSchemeTypes.EULER.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            s = s0
            v = v0
            for i_step in range(1, num_steps + 1):
                z1 = np.random.normal(0.0, 1.0) * sdt
                z2 = np.random.normal(0.0, 1.0) * sdt
                z_v = z1
                z_s = rho * z1 + rhohat * z2
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                v += kappa * (theta - vplus) * dt + sigma * rtvplus * z_v + 0.25 * sigma2 * (z_v * z_v - dt)
                s += (r - q) * s * dt + rtvplus * s * z_s + 0.5 * s * vplus * (z_v * z_v - dt)
                s_paths[i_path, i_step] = s

    elif scheme == HestonNumericalSchemeTypes.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for i_path in range(0, num_paths):
            x = log(s0)
            v = v0
            for i_step in range(1, num_steps + 1):
                z_v = np.random.normal(0.0, 1.0) * sdt
                z_s = rho * z_v + rhohat * np.random.normal(0.0, 1.0) * sdt
                vplus = max(v, 0.0)
                rtvplus = np.sqrt(vplus)
                x += (r - q - 0.5 * vplus) * dt + rtvplus * z_s
                v += kappa * (theta - vplus) * dt + sigma * rtvplus * z_v + sigma2 * (z_v * z_v - dt) / 4.0
                s_paths[i_path, i_step] = exp(x)

    elif scheme == HestonNumericalSchemeTypes.QUADEXP.value:
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
            for i_step in range(1, num_steps + 1):
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

                x += mu * dt + k_0 + (k_1 * vn + k_2 * vnp) + np.sqrt(k_3 * vn + k_4 * vnp) * z_s
                s_paths[i_path, i_step] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown HestonNumericalSchme")

    return s_paths


########################################################################################


class Heston:

    def __init__(self, v0: float, kappa: float, theta: float, xi: float, rho: float):

        if v0 < 0.0:
            raise FinError("Initial variance must be non-negative.")

        if kappa <= 0.0:
            raise FinError("Mean-reversion speed must be positive.")

        if theta < 0.0:
            raise FinError("Long-run variance must be non-negative.")

        if xi <= 0.0:
            raise FinError("Volatility of variance must be positive.")

        if rho < -1.0 or rho > 1.0:
            raise FinError("Correlation must lie between -1 and 1.")

        self._v0 = v0
        self._kappa = kappa
        self._theta = theta
        self._xi = xi
        self._rho = rho

    ####################################################################################

    def value(
        self,
        stock_price,
        t_exp,
        strike,
        option_type,
        interest_rate,
        dividend_yield,
        method=HestonValueTypes.LEWIS,
    ):

        call_value = self.call_value(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
            method,
        )

        if option_type == OptionTypes.EUROPEAN_CALL.value:
            return call_value

        if option_type == OptionTypes.EUROPEAN_PUT.value:
            return call_value - stock_price * exp(-dividend_yield * t_exp) + strike * exp(-interest_rate * t_exp)

        raise FinError("Unsupported option type.")

    ####################################################################################

    def call_value(
        self,
        stock_price,
        t_exp,
        strike,
        interest_rate,
        dividend_yield,
        method=HestonValueTypes.LEWIS,
    ):

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive.")

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

        if method == HestonValueTypes.LEWIS:
            return self.value_call_lewis(
                t_exp,
                strike,
                stock_price,
                interest_rate,
                dividend_yield,
            )

        elif method == HestonValueTypes.LEWIS_ROUAH:
            return self.value_call_lewis_rouah(
                t_exp,
                strike,
                stock_price,
                interest_rate,
                dividend_yield,
            )

        elif method == HestonValueTypes.GATHERAL:
            return self.value_call_gatheral(
                t_exp,
                strike,
                stock_price,
                interest_rate,
                dividend_yield,
            )

        elif method == HestonValueTypes.WEBER:
            return self.value_call_weber(
                t_exp,
                strike,
                stock_price,
                interest_rate,
                dividend_yield,
            )

        raise FinError("Unknown Heston valuation method.")

    ####################################################################################

    def value_mc(
        self,
        stock_price: float,
        t_exp: float,
        strike: float,
        option_type: int,
        interest_rate: float,
        dividend_yield: float,
        num_paths: int,
        num_steps_per_year: int,
        seed: int,
        scheme=HestonNumericalSchemeTypes.EULERLOG,
    ):

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        if stock_price <= 0.0:
            raise FinError("Stock price must be positive.")

        if strike <= 0.0:
            raise FinError("Strike must be positive.")

        if num_paths <= 0:
            raise FinError("Number of paths must be positive.")

        if num_steps_per_year <= 0:
            raise FinError("Number of steps per year must be positive.")

        tau = t_exp
        k = strike
        dt = 1.0 / num_steps_per_year
        scheme_value = float(scheme.value)

        s_paths = get_paths(
            stock_price,
            interest_rate,
            dividend_yield,
            self._v0,
            self._kappa,
            self._theta,
            self._xi,
            self._rho,
            tau,
            dt,
            num_paths,
            seed,
            scheme_value,
        )

        if option_type == OptionTypes.EUROPEAN_CALL.value:
            path_payoff = np.maximum(s_paths[:, -1] - k, 0.0)
        elif option_type == OptionTypes.EUROPEAN_PUT.value:
            path_payoff = np.maximum(k - s_paths[:, -1], 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(path_payoff)
        v = payoff * exp(-interest_rate * tau)
        return v

    ####################################################################################

    def value_call_lewis(
        self,
        t_exp: float,
        strike: float,
        stock_price: float,
        interest_rate: float,
        dividend_yield: float,
    ) -> float:

        tau = t_exp

        rho = self._rho
        xi = self._xi
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        r = interest_rate
        q = dividend_yield
        s0 = stock_price
        kk = strike
        ff = s0 * exp((r - q) * tau)
        vv = xi * xi

        def phi(
            k_in,
        ):
            k = k_in + 0.5 * 1j
            b = kappa + 1j * rho * xi * k
            d = np.sqrt(b**2 + vv * k * (k - 1j))
            g = (b - d) / (b + d)
            t_m = (b - d) / vv
            qq = np.exp(-d * tau)
            t = t_m * (1.0 - qq) / (1.0 - g * qq)
            ww = kappa * theta * (tau * t_m - 2.0 * np.log((1.0 - g * qq) / (1.0 - g)) / vv)
            phi = np.exp(ww + v0 * t)
            return phi

        def phi_transform(x):
            def integrand(k):
                return 2.0 * np.real(np.exp(-1j * k * x) * phi(k)) / (k**2 + 1.0 / 4.0)

            return integrate.quad(integrand, 0, np.inf)[0]

        x = log(ff / kk)
        i_1 = phi_transform(x) / (2.0 * pi)
        v1 = ff * exp(-r * tau) - np.sqrt(kk * ff) * exp(-r * tau) * i_1
        return v1

    ####################################################################################

    def value_call_lewis_rouah(
        self,
        t_exp: float,
        strike: float,
        stock_price: float,
        interest_rate: float,
        dividend_yield: float,
    ) -> float:

        tau = t_exp

        rho = self._rho
        xi = self._xi
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        vv = xi * xi

        s0 = stock_price
        f = s0 * exp((r - q) * tau)
        k = strike
        x = log(f / k)

        def fn(k_in):
            k = k_in + 0.5 * 1j
            b = (2.0 / vv) * (1j * k * rho * xi + kappa)
            e = np.sqrt(b**2 + 4.0 * k * (k - 1j) / vv)
            g = (b - e) / 2.0
            h = (b - e) / (b + e)
            q = vv * tau / 2.0
            qq = np.exp(-e * q)
            hh = np.exp(
                (2.0 * kappa * theta / vv) * (q * g - np.log((1.0 - h * qq) / (1.0 - h)))
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

    def value_call_weber(
        self,
        t_exp: float,
        strike: float,
        stock_price: float,
        interest_rate: float,
        dividend_yield: float,
    ) -> float:

        tau = t_exp

        rho = self._rho
        xi = self._xi
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        s0 = stock_price
        k = strike
        vv = xi**2

        def fn(s, b):
            def integrand(u):
                beta = b - 1j * rho * xi * u
                d = np.sqrt((beta**2) - vv * u * (s * 1j - u))
                g = (beta - d) / (beta + d)
                qq = np.exp(-d * tau)
                bb = (beta - d) * (1.0 - qq) / (1.0 - g * qq) / vv
                aa = kappa * ((beta - d) * tau - 2.0 * np.log((1.0 - g * qq) / (1.0 - g))) / vv
                v = np.exp(aa * theta + bb * v0 + 1j * u * np.log(s0 / (k * np.exp(-(r - q) * tau)))) / (u * 1j)
                return v.real

            area = 0.50 + (1.0 / pi) * integrate.quad(integrand, 0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * fn(1.0, kappa - rho * xi)
        v = v - exp(-r * tau) * k * fn(-1.0, kappa)

        return v

    #################################### @#@#############################################
    # Gatheral book page 19 with definition of x given on page 16 and noting
    # that the value C is a forward value and so needs to be discounted
    ####################################################################################

    def value_call_gatheral(
        self,
        t_exp: float,
        strike: float,
        stock_price: float,
        interest_rate: float,
        dividend_yield: float,
    ) -> float:

        tau = t_exp
        rho = self._rho
        xi = self._xi
        v0 = self._v0
        kappa = self._kappa
        theta = self._theta

        q = dividend_yield
        r = interest_rate
        s0 = stock_price
        k = strike
        f = s0 * exp((r - q) * tau)
        x0 = log(f / k)

        def ff(j):
            def integrand(u):
                vv = xi * xi
                aa = -u * u / 2.0 - 1j * u / 2.0 + 1j * j * u
                bb = kappa - rho * xi * j - rho * xi * 1j * u
                gg = vv / 2.0
                d = np.sqrt(bb**2 - 4.0 * aa * gg)
                rplus = (bb + d) / 2.0 / gg
                rminus = (bb - d) / 2.0 / gg
                rr = rminus / rplus
                qq = np.exp(-d * tau)
                dd = rminus * (1.0 - qq) / (1.0 - rr * qq)
                cc = kappa * (rminus * tau - (2.0 / vv) * np.log((1.0 - rr * qq) / (1.0 - rr)))
                phi = np.exp(cc * theta + dd * v0 + 1j * u * x0) / (1j * u)
                return phi.real

            area = 0.50 + 1.0 / pi * integrate.quad(integrand, 0.0, np.inf)[0]
            return area

        v = s0 * exp(-q * tau) * ff(1) - k * exp(-r * tau) * ff(0)
        return v

    ####################################################################################

    def feller_condition(self) -> bool:
        return 2.0 * self._kappa * self._theta >= self._xi**2

    ####################################################################################

    def implied_volatility(
        self,
        stock_price,
        t_exp,
        strike,
        interest_rate,
        dividend_yield,
        method=HestonValueTypes.LEWIS,
    ):

        price = self.call_value(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
            method,
        )

        return implied_volatility(
            stock_price,
            t_exp,
            strike,
            interest_rate,
            dividend_yield,
            price,
            OptionTypes.EUROPEAN_CALL.value,
        )

    ####################################################################################

    def volatility_smile(
        self,
        t_exp,
        strikes,
        stock_price,
        interest_rate,
        dividend_yield,
        method=HestonValueTypes.LEWIS,
    ):

        strikes = np.asarray(strikes, dtype=float)

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if np.any(strikes <= 0.0):
            raise FinError("Strikes must be positive.")

        if t_exp <= 0.0:
            raise FinError("Time to expiry must be positive.")

        implied_vols = np.empty(len(strikes))

        for i, strike in enumerate(strikes):

            implied_vols[i] = self.implied_volatility(
                stock_price,
                t_exp,
                strike,
                interest_rate,
                dividend_yield,
                method,
            )

        return implied_vols


########################################################################################
