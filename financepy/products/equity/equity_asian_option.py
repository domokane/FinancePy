##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

import numpy as np
from numba import njit

# TODO: Add perturbatory risk using the analytical methods !!
# TODO: Add Sobol to Monte Carlo

from ...utils.math import covar
from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.error import FinError

from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve

from ...utils.math import normcdf


########################################################################################


class AsianOptionValuationMethods(Enum):
    GEOMETRIC = (1,)
    TURNBULL_WAKEMAN = (2,)
    CURRAN = 3


########################################################################################


error_str = "In averaging period so need to enter accrued average."

########################################################################################
# An Asian option on an arithmetic average and strike K has a payoff
# Max(SA(T)-K,0) where SA is the arithmetic average
# We define three dates
# - Valuation date for which we want the price
# - Start Averaging Date for when the averaging starts
# - Expiry date for when the payoff is made and the option expires
#
# In the model we have
# tv = is the time now
# t0 = time to the start averaging date in years
# t = time to the expiry date in years
# tau = length of averaging period in years at the start of the option
#
# We can be before the start of the averaging period in which case t0 > 0
# We can be after the start of the averaging period in which case we set t0=0
# and we note that t <= tau
#
# If we are in the averaging period then we need to know the accrued average
# I call this AA and the new average is now given by the accrued average
# The option payoff is now Max( (AA x (tau-t) + SA(t0) x (t-t0))/tau - K,0)
# This simplifies to
#
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#
########################################################################################


@njit(cache=True, fastmath=True)
def _value_mc_numba(
    t0,
    t,
    tau,
    k,
    n,
    opt_type,
    stock_price,
    interest_rate,
    dividend_yield,
    volatility,
    num_paths,
    seed,
    accrued_average,
):

    # Start pricing here
    np.random.seed(seed)
    multiplier = 1.0

    if t0 < 0.0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    mu = interest_rate - dividend_yield
    v2 = volatility**2
    dt = (t - t0) / n

    payoff_a = 0.0

    for _ in range(0, num_paths):

        # evolve stock price to start of averaging period
        g = np.random.standard_normal(1)
        s_1 = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 + g[0] * np.sqrt(t0) * volatility
        )
        s_2 = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 - g[0] * np.sqrt(t0) * volatility
        )

        # enter averaging period
        s_1_arithmetic = 0.0
        s_2_arithmetic = 0.0

        g = np.random.standard_normal(n)

        for obs in range(0, n):

            s_1 = s_1 * np.exp((mu - v2 / 2.0) * dt + g[obs] * np.sqrt(dt) * volatility)

            s_2 = s_2 * np.exp((mu - v2 / 2.0) * dt - g[obs] * np.sqrt(dt) * volatility)

            s_1_arithmetic += s_1
            s_2_arithmetic += s_2

        s_1_arithmetic /= n
        s_2_arithmetic /= n

        if opt_type == OptionTypes.EUROPEAN_CALL:
            payoff_a += max(s_1_arithmetic - k, 0.0)
            payoff_a += max(s_2_arithmetic - k, 0.0)
        elif opt_type == OptionTypes.EUROPEAN_PUT:
            payoff_a += max(k - s_1_arithmetic, 0.0)
            payoff_a += max(k - s_2_arithmetic, 0.0)
        else:
            return None

    v_a = payoff_a * np.exp(-interest_rate * t) / num_paths / 2.0
    v_a = v_a * multiplier
    return v_a


########################################################################################


@njit(cache=True, fastmath=True)
def _value_mc_fast_numba(
    t0: float,
    t: float,
    tau: float,
    k: float,
    n: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths: int,
    seed: int,
    accrued_average: float,
):

    np.random.seed(seed)
    mu = interest_rate - dividend_yield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interest_rate
    num_paths = int(num_paths)

    multiplier = 1.0

    if t0 < 0.0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    # g = np.random.normal(0.0, 1.0, size=(num_paths))

    gg = np.empty(num_paths, np.float64)
    for i_path in range(0, num_paths):
        rv = np.random.normal()
        gg[i_path] = rv

    s_1 = np.empty(num_paths)
    s_2 = np.empty(num_paths)

    for ip in range(0, num_paths):
        s_1[ip] = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 + gg[ip] * np.sqrt(t0) * volatility
        )
        s_2[ip] = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 - gg[ip] * np.sqrt(t0) * volatility
        )

    s_1_arithmetic = np.zeros(num_paths)
    s_2_arithmetic = np.zeros(num_paths)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths)

        for ip in range(0, num_paths):
            s_1[ip] = s_1[ip] * np.exp(
                (mu - v2 / 2.0) * dt + g[ip] * np.sqrt(dt) * volatility
            )
            s_2[ip] = s_2[ip] * np.exp(
                (mu - v2 / 2.0) * dt - g[ip] * np.sqrt(dt) * volatility
            )

        for ip in range(0, num_paths):
            s_1_arithmetic[ip] += s_1[ip] / n
            s_2_arithmetic[ip] += s_2[ip] / n

    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
    else:
        raise FinError("Unknown option type.")

    payoff_a = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    v_a = multiplier * payoff_a * np.exp(-r * t) / 2.0
    return v_a


########################################################################################


@njit(cache=True, fastmath=True)
def _value_mc_fast_cv_numba(
    t0,
    t,
    tau,
    k,
    n,
    opt_type,
    stock_price,
    interest_rate,
    dividend_yield,
    volatility,
    num_paths,
    seed,
    accrued_average,
    v_g_exact,
):

    np.random.seed(seed)
    mu = interest_rate - dividend_yield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interest_rate

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    g = np.random.normal(0.0, 1.0, size=num_paths)

    s_1 = np.empty(num_paths)
    s_2 = np.empty(num_paths)

    for ip in range(0, num_paths):
        s_1[ip] = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 + g[ip] * np.sqrt(t0) * volatility
        )
        s_2[ip] = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 - g[ip] * np.sqrt(t0) * volatility
        )

    s_1_arithmetic = np.zeros(num_paths)
    s_2_arithmetic = np.zeros(num_paths)
    ln_s_1_geometric = np.zeros(num_paths)
    ln_s_2_geometric = np.zeros(num_paths)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths)
        for ip in range(0, num_paths):
            s_1[ip] = s_1[ip] * np.exp(
                (mu - v2 / 2.0) * dt + g[ip] * np.sqrt(dt) * volatility
            )
            s_2[ip] = s_2[ip] * np.exp(
                (mu - v2 / 2.0) * dt - g[ip] * np.sqrt(dt) * volatility
            )

        for ip in range(0, num_paths):
            s_1_arithmetic[ip] += s_1[ip]
            s_2_arithmetic[ip] += s_2[ip]
            ln_s_1_geometric[ip] += np.log(s_1[ip])
            ln_s_2_geometric[ip] += np.log(s_2[ip])

    s_1_geometric = np.empty(num_paths)
    s_2_geometric = np.empty(num_paths)

    for ip in range(0, num_paths):
        s_1_arithmetic[ip] /= n
        s_1_geometric[ip] = np.exp(ln_s_1_geometric[ip] / n)
        s_2_arithmetic[ip] /= n
        s_2_geometric[ip] = np.exp(ln_s_2_geometric[ip] / n)

    if opt_type == OptionTypes.EUROPEAN_CALL:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_g_1 = np.maximum(s_1_geometric - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
        payoff_g_2 = np.maximum(s_2_geometric - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_g_1 = np.maximum(k - s_1_geometric, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
        payoff_g_2 = np.maximum(k - s_2_geometric, 0.0)
    else:
        raise FinError("Unknown Option Type")

    payoff_a = np.concatenate((payoff_a_1, payoff_a_2), axis=0)
    payoff_g = np.concatenate((payoff_g_1, payoff_g_2), axis=0)

    # Now we do the control variate adjustment
    m = covar(payoff_a, payoff_a)

    if np.abs(m[1][1]) < 1e-10:
        lam = 0.0
    else:
        lam = m[0][1] / m[1][1]

    payoff_a_mean = np.mean(payoff_a)
    payoff_g_mean = np.mean(payoff_g)

    v_a = payoff_a_mean * np.exp(-r * t) * multiplier
    v_g = payoff_g_mean * np.exp(-r * t) * multiplier

    epsilon = v_g_exact - v_g
    v_a_cv = v_a + lam * epsilon

    return v_a_cv


########################################################################################


class EquityAsianOption:
    """Class for an Equity Asian Option. This is an option with a final payoff
    linked to the averaging of the stock price over some specified period
    before the option expires. The valuation is done for both an arithmetic and
    a geometric average but the former can only be done either using an
    analytical approximation of the arithmetic average distribution or by using
    Monte-Carlo simulation."""

    def __init__(
        self,
        start_averaging_dt: Date,
        expiry_dt: Date,
        strike_price: float,
        opt_type: OptionTypes,
        num_obs: int = 100,
    ):
        """Create an EquityAsian option object which takes a start date for
        the averaging, an expiry date, a strike price, an option type and a
        number of observations."""

        check_argument_types(self.__init__, locals())

        if start_averaging_dt > expiry_dt:
            raise FinError("Averaging starts after expiry date")

        self.start_averaging_date = start_averaging_dt
        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.opt_type = opt_type
        self.num_observations = num_obs

    ####################################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        method: AsianOptionValuationMethods,
        accrued_average: float = None,
    ):
        """Calculate the value of an Asian option using one of the specified
        analytical approximations for an average rate option. These are the
        three enumerated values in the enum AsianOptionValuationMethods. The
        choices of approximation are (i) GEOMETRIC - the average is a geometric
        one as in paper by Kenna and Worst (1990), (ii) TURNBULL_WAKEMAN -
        this is a value based on an edgeworth expansion of the moments of the
        arithmetic average, and (iii) CURRAN - another approximative approach
        by Curran based on conditioning on the geometric mean price. Just
        choose the corresponding enumerated value to switch between these
        different approaches.

        Note that the accrued average is only required if the value date is
        inside the averaging period for the option."""

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date"
            )

        if method == AsianOptionValuationMethods.GEOMETRIC:
            v = self._value_geometric(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                accrued_average,
            )

        elif method == AsianOptionValuationMethods.TURNBULL_WAKEMAN:
            v = self._value_turnbull_wakeman(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                accrued_average,
            )

        elif method == AsianOptionValuationMethods.CURRAN:
            v = self._value_curran(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                accrued_average,
            )
        else:
            raise FinError("Unknown valuation model")

        return v

    ####################################################################################

    def _value_geometric(
        self,
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        accrued_average,
    ):
        """This option valuation is based on paper by Kemna and Vorst 1990. It
        calculates the Geometric Asian option price which is a lower bound on
        the Arithmetic option price. This should not be used as a valuation
        model for the Arithmetic Average option but can be used as a control
        variate for other approaches."""

        if value_dt > self.expiry_dt:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        k = self.strike_price
        n = self.num_observations
        s0 = stock_price

        multiplier = 1.0

        if t0 < 0:  # we are in the averaging period

            if accrued_average is None:
                raise FinError(error_str)

            # we adjust the strike to account for the accrued coupon
            k = (k * tau + accrued_average * t0) / t_exp
            # the number of options is rescaled also
            multiplier = t_exp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled
            n = n * t_exp / tau

        sig_sq = volatility**2
        mean_geo = (r - q - sig_sq / 2.0) * (t0 + (t_exp - t0) / 2.0)
        var_geo = sig_sq * (t0 + (t_exp - t0) * (2 * n - 1) / (6 * n))
        eg = s0 * np.exp(mean_geo + var_geo / 2.0)

        if np.abs(var_geo) < 1e-10:
            raise FinError("Asian option geometric variance is zero.")

        d1 = (mean_geo + np.log(s0 / k) + var_geo) / np.sqrt(var_geo)
        d2 = d1 - np.sqrt(var_geo)

        # the Geometric price is the lower bound
        call_g = np.exp(-r * t_exp) * (eg * normcdf(d1) - k * normcdf(d2))

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            v = call_g
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            put_g = call_g - (eg - k) * np.exp(-r * t_exp)
            v = put_g
        else:
            raise FinError("Unknown option type " + str(self.opt_type))

        v = v * multiplier
        return v

    ####################################################################################

    def _value_curran(
        self,
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        accrued_average,
    ):
        """Valuation of an Asian option using the result by Vorst."""

        if value_dt > self.expiry_dt:
            raise FinError("Value date after option expiry date.")

        # the years to the start of the averaging period
        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        multiplier = 1.0

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        s0 = stock_price
        b = r - q
        sigma2 = volatility**2
        k = self.strike_price

        n = self.num_observations

        if t0 < 0:  # we are in the averaging period

            if accrued_average is None:
                raise FinError(error_str)

            # we adjust the strike to account for the accrued coupon
            k = (k * tau + accrued_average * t0) / t_exp
            # the number of options is rescaled also
            multiplier = t_exp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * t_exp / tau + 0.5) + 1

        h = (t_exp - t0) / (n - 1)
        u = (1.0 - np.exp(b * h * n)) / (1.0 - np.exp(b * h))
        w = (1.0 - np.exp((2 * b + sigma2) * h * n)) / (
            1.0 - np.exp((2 * b + sigma2) * h)
        )

        fa = (s0 / n) * np.exp(b * t0) * u
        ea2 = (s0 * s0 / n / n) * np.exp((2.0 * b + sigma2) * t0)
        ea2 = ea2 * (w + 2.0 / (1.0 - np.exp((b + sigma2) * h)) * (u - w))
        sigma_aa = np.sqrt((np.log(ea2) - 2.0 * np.log(fa)) / t_exp)

        d1 = (np.log(fa / k) + sigma_aa * sigma_aa * t_exp / 2.0) / (
            sigma_aa * np.sqrt(t_exp)
        )
        d2 = d1 - sigma_aa * np.sqrt(t_exp)

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            v = np.exp(-r * t_exp) * (fa * normcdf(d1) - k * normcdf(d2))
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            v = np.exp(-r * t_exp) * (k * normcdf(-d2) - fa * normcdf(-d1))
        else:
            return None

        v = v * multiplier
        return v

    ####################################################################################

    def _value_turnbull_wakeman(
        self,
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        accrued_average,
    ):
        """Asian option valuation based on paper by Turnbull and Wakeman 1991
        which uses the edgeworth expansion to find the first two moments of the
        arithmetic average."""

        if value_dt > self.expiry_dt:
            raise FinError("Value date after option expiry date.")

        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        k = self.strike_price
        multiplier = 1.0
        n = self.num_observations

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        if t0 < 0:  # we are in the averaging period

            if accrued_average is None:
                raise FinError(error_str)

            # we adjust the strike to account for the accrued coupon
            k = (k * tau + accrued_average * t0) / t_exp
            # the number of options is rescaled also
            multiplier = t_exp / tau
            # there is no pre-averaging time
            t0 = 0.0
            # the number of observations is scaled and floored at 1
            n = int(n * t_exp / tau + 0.5) + 1

        # need to handle this
        b = r - q
        sigma2 = volatility**2
        a1 = b + sigma2
        a2 = 2 * b + sigma2
        s0 = stock_price

        dt = t_exp - t0

        if b == 0:
            m1 = 1.0
            m2 = 2.0 * np.exp(sigma2 * t_exp) - 2.0 * np.exp(sigma2 * t0) * (
                1.0 + sigma2 * dt
            )
            m2 = m2 / sigma2 / sigma2 / dt / dt
        else:
            m1 = s0 * (np.exp(b * t_exp) - np.exp(b * t0)) / (b * dt)
            m2 = np.exp(a2 * t_exp) / a1 / a2 / dt / dt + (
                np.exp(a2 * t0) / b / dt / dt
            ) * (1.0 / a2 - np.exp(b * dt) / a1)
            m2 = 2.0 * m2 * s0 * s0

        f0 = m1
        sigma2 = (1.0 / t_exp) * np.log(m2 / m1 / m1)
        sigma = np.sqrt(sigma2)

        d1 = (np.log(f0 / k) + sigma2 * t_exp / 2) / sigma / np.sqrt(t_exp)
        d2 = d1 - sigma * np.sqrt(t_exp)

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            call = np.exp(-r * t_exp) * (f0 * normcdf(d1) - k * normcdf(d2))
            v = call
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            put = np.exp(-r * t_exp) * (k * normcdf(-d2) - f0 * normcdf(-d1))
            v = put
        else:
            return None

        v = v * multiplier
        return v

    ####################################################################################

    def _value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_paths: int,
        seed: int,
        accrued_average: float,
    ):
        """Monte Carlo valuation of the Asian Average option using standard
        Monte Carlo code enhanced by Numba. I have discontinued the use of this
        as it is both slow and has limited variance reduction."""

        # Basic validation
        if value_dt > self.expiry_dt:
            raise FinError("Value date after option expiry date.")

        if value_dt > self.start_averaging_date and accrued_average is None:
            raise FinError(error_str)

        # the years to the start of the averaging period
        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        k = self.strike_price
        n = self.num_observations

        v = _value_mc_numba(
            t0,
            t_exp,
            tau,
            k,
            n,
            self.opt_type,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
        )

        return v

    ####################################################################################

    def _value_mc_fast(
        self,
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,  # Yield
        model,  # Model
        num_paths,  # Numpaths integer
        seed,
        accrued_average,
    ):
        """Monte Carlo valuation of the Asian Average option. This method uses
        a lot of Numpy vectorisation. It is also helped by Numba."""

        # the years to the start of the averaging period
        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        k = self.strike_price
        n = self.num_observations

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        v = _value_mc_fast_numba(
            t0,
            t_exp,
            tau,
            k,
            n,
            self.opt_type.value,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
        )

        return v

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_paths: int,
        seed: int,
        accrued_average: float,
    ):
        """Monte Carlo valuation of the Asian Average option using a control
        variate method that improves accuracy and reduces the variance of the
        price. This uses Numpy and Numba. This is the standard MC pricer."""

        # the years to the start of the averaging period
        t0 = (self.start_averaging_date - value_dt) / G_DAYS_IN_YEARS
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        tau = (self.expiry_dt - self.start_averaging_date) / G_DAYS_IN_YEARS

        k = self.strike_price
        n = self.num_observations

        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        volatility = model.volatility

        # For control variate we price a Geometric average option exactly
        v_g_exact = self._value_geometric(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            accrued_average,
        )

        v = _value_mc_fast_cv_numba(
            t0,
            t_exp,
            tau,
            k,
            n,
            self.opt_type,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
            v_g_exact,
        )

        return v

    ####################################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START AVERAGING DATE", self.start_averaging_date)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("NUM OBSERVATIONS", self.num_observations, "")
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
