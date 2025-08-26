##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
from math import exp, log, sqrt
import numpy as np

from ...utils.error import FinError
from ...utils.math import normcdf
from ...utils.global_vars import G_DAYS_IN_YEARS
from ...products.fx.fx_option import FXOption
from ...models.process_simulator import FinProcessSimulator
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date


########################################################################################


class FinFXBarrierTypes(Enum):
    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8


########################################################################################


class FXBarrierOption(FXOption):

    def __init__(
        self,
        expiry_dt: Date,
        strike_fx_rate: float,  # 1 unit of foreign in domestic
        currency_pair: str,  # FORDOM
        opt_type: FinFXBarrierTypes,
        barrier_level: float,
        num_obs_per_year: int,
        notional: float,
        notional_currency: str,
    ):
        """Create FX Barrier option product. This is an option that cancels if
        the FX rate crosses a barrier during the life of the option."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.strike_fx_rate = float(strike_fx_rate)
        self.currency_pair = currency_pair
        self.barrier_level = float(barrier_level)
        self.num_obs_per_year = int(num_obs_per_year)
        self.opt_type = opt_type
        self.notional = notional
        self.notional_currency = notional_currency

    ##########################################################################

    def value(self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model):
        """Value FX Barrier Option using Black-Scholes model with closed-form
        analytical models."""

        # This prices the option using the formulae given in the paper
        # by Clewlow, Llanos and Strickland December 1994 which can be found at
        # https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as option value date"
            )

        if foreign_curve.value_dt != value_dt:
            raise FinError("Foreign Curve valuation date not same as option value date")

        k = self.strike_fx_rate
        s0 = spot_fx_rate
        h = self.barrier_level

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        ln_s0_k = log(float(s0) / k)
        sqrt_t = sqrt(t)

        dq = foreign_curve.df_t(t)
        df = domestic_curve.df_t(t)
        r_d = -log(df) / t
        rf = -log(dq) / t

        volatility = model.volatility
        sigma_root_t = volatility * sqrt_t
        v2 = volatility * volatility
        mu = r_d - rf
        d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / sigma_root_t
        d2 = (ln_s0_k + (mu - v2 / 2.0) * t) / sigma_root_t

        c = s0 * dq * normcdf(d1) - k * df * normcdf(d2)
        p = k * df * normcdf(-d2) - s0 * dq * normcdf(-d1)
        #        print("CALL:",c,"PUT:",p)

        if self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL and s0 <= h:
            return 0.0
        if self.opt_type == FinFXBarrierTypes.UP_AND_OUT_CALL and s0 >= h:
            return 0.0
        if self.opt_type == FinFXBarrierTypes.UP_AND_OUT_PUT and s0 >= h:
            return 0.0
        if self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT and s0 <= h:
            return 0.0
        if self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_CALL and s0 <= h:
            return c
        if self.opt_type == FinFXBarrierTypes.UP_AND_IN_CALL and s0 >= h:
            return c
        if self.opt_type == FinFXBarrierTypes.UP_AND_IN_PUT and s0 >= h:
            return p
        if self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_PUT and s0 <= h:
            return p

        num_observations = t * self.num_obs_per_year

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        if self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
        else:
            raise FinError("Unknown barrier option type." + str(self.opt_type))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        ll = (mu + v2 / 2.0) / v2
        y = log(h * h / (s0 * k)) / sigma_root_t + ll * sigma_root_t
        x1 = log(s0 / h) / sigma_root_t + ll * sigma_root_t
        y1 = log(h / s0) / sigma_root_t + ll * sigma_root_t
        h_over_s = h / s0

        if self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            if h >= k:
                c_do = (
                    s0 * dq * normcdf(x1)
                    - k * df * normcdf(x1 - sigma_root_t)
                    - s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
                    + k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * normcdf(y1 - sigma_root_t)
                )
                price = c_do
            else:
                c_di = s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
                    h_over_s, 2.0 * ll - 2.0
                ) * normcdf(y - sigma_root_t)
                price = c - c_di
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            if h <= k:
                c_di = s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
                    h_over_s, 2.0 * ll - 2.0
                ) * normcdf(y - sigma_root_t)
                price = c_di
            else:
                c_do = (
                    s0 * dq * normcdf(x1)
                    - k * df * normcdf(x1 - sigma_root_t)
                    - s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
                    + k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * normcdf(y1 - sigma_root_t)
                )
                price = c - c_do
        elif self.opt_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            if h >= k:
                c_ui = (
                    s0 * dq * normcdf(x1)
                    - k * df * normcdf(x1 - sigma_root_t)
                    - s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
                    + k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * (normcdf(-y + sigma_root_t) - normcdf(-y1 + sigma_root_t))
                )
                price = c_ui
            else:
                price = c
        elif self.opt_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            if h > k:
                c_ui = (
                    s0 * dq * normcdf(x1)
                    - k * df * normcdf(x1 - sigma_root_t)
                    - s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
                    + k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * (normcdf(-y + sigma_root_t) - normcdf(-y1 + sigma_root_t))
                )
                price = c - c_ui
            else:
                price = 0.0
        elif self.opt_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            if h > k:
                p_ui = -s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
                    h_over_s, 2.0 * ll - 2.0
                ) * normcdf(-y + sigma_root_t)
                price = p_ui
            else:
                p_uo = (
                    -s0 * dq * normcdf(-x1)
                    + k * df * normcdf(-x1 + sigma_root_t)
                    + s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
                    - k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * normcdf(-y1 + sigma_root_t)
                )
                price = p - p_uo
        elif self.opt_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            if h >= k:
                p_ui = -s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
                    h_over_s, 2.0 * ll - 2.0
                ) * normcdf(-y + sigma_root_t)
                price = p - p_ui
            else:
                p_uo = (
                    -s0 * dq * normcdf(-x1)
                    + k * df * normcdf(-x1 + sigma_root_t)
                    + s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
                    - k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * normcdf(-y1 + sigma_root_t)
                )
                price = p_uo
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            if h >= k:
                price = 0.0
            else:
                p_di = (
                    -s0 * dq * normcdf(-x1)
                    + k * df * normcdf(-x1 + sigma_root_t)
                    + s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
                    - k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * (normcdf(y - sigma_root_t) - normcdf(y1 - sigma_root_t))
                )
                price = p - p_di
        elif self.opt_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            if h >= k:
                price = p
            else:
                p_di = (
                    -s0 * dq * normcdf(-x1)
                    + k * df * normcdf(-x1 + sigma_root_t)
                    + s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
                    - k
                    * df
                    * pow(h_over_s, 2.0 * ll - 2.0)
                    * (normcdf(y - sigma_root_t) - normcdf(y1 - sigma_root_t))
                )
                price = p_di
        else:
            raise FinError("Unknown barrier option type." + str(self.opt_type))

        return price

    ###########################################################################

    def value_mc(
        self,
        value_dt,
        spot_fx_rate,
        dom_interest_rate,
        process_type,
        model_params,
        num_ann_steps=552,
        num_paths=5000,
        seed=4242,
    ):
        """Value the FX Barrier Option using Monte Carlo."""

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        num_time_steps = int(t * num_ann_steps)
        k = self.strike_fx_rate
        b = self.barrier_level
        s0 = spot_fx_rate
        opt_type = self.opt_type

        process = FinProcessSimulator()

        r_d = dom_interest_rate

        #######################################################################

        if opt_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL and s0 <= b:
            return 0.0
        elif opt_type == FinFXBarrierTypes.UP_AND_OUT_CALL and s0 >= b:
            return 0.0
        elif opt_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT and s0 <= b:
            return 0.0
        elif opt_type == FinFXBarrierTypes.UP_AND_OUT_PUT and s0 >= b:
            return 0.0

        #######################################################################

        simple_call = False
        simple_put = False

        if opt_type == FinFXBarrierTypes.DOWN_AND_IN_CALL and s0 <= b:
            simple_call = True
        elif opt_type == FinFXBarrierTypes.UP_AND_IN_CALL and s0 >= b:
            simple_call = True
        elif opt_type == FinFXBarrierTypes.UP_AND_IN_PUT and s0 >= b:
            simple_put = True
        elif opt_type == FinFXBarrierTypes.DOWN_AND_IN_PUT and s0 <= b:
            simple_put = True

        if simple_put or simple_call:
            s_all = process.get_process(
                process_type, t, model_params, 1, num_paths, seed
            )

            if simple_call:
                s_t = s_all[:, -1]
                c = (np.maximum(s_t - k, 0.0)).mean()
                c = c * exp(-r_d * t)
                return c

            if simple_put:
                s_t = s_all[:, -1]
                p = (np.maximum(k - s_t, 0.0)).mean()
                p = p * exp(-r_d * t)
                return p

        # Otherwise get full set of paths
        s_all = process.get_process(
            process_type, t, model_params, num_time_steps, num_paths, seed
        )

        (num_paths, num_time_steps) = s_all.shape

        if opt_type in (
            FinFXBarrierTypes.DOWN_AND_IN_CALL,
            FinFXBarrierTypes.DOWN_AND_OUT_CALL,
            FinFXBarrierTypes.DOWN_AND_IN_PUT,
            FinFXBarrierTypes.DOWN_AND_OUT_PUT,
        ):

            barrier_crossed_from_above = [False] * num_paths

            for p in range(0, num_paths):
                barrier_crossed_from_above[p] = np.any(s_all[p] <= b)

        if opt_type in (
            FinFXBarrierTypes.UP_AND_IN_CALL,
            FinFXBarrierTypes.UP_AND_OUT_CALL,
            FinFXBarrierTypes.UP_AND_IN_PUT,
            FinFXBarrierTypes.UP_AND_OUT_PUT,
        ):

            barrier_crossed_from_below = [False] * num_paths
            for p in range(0, num_paths):
                barrier_crossed_from_below[p] = np.any(s_all[p] >= b)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if opt_type == FinFXBarrierTypes.DOWN_AND_OUT_CALL:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * (
                ones - barrier_crossed_from_above
            )
        elif opt_type == FinFXBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_above
        elif opt_type == FinFXBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_below
        elif opt_type == FinFXBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(s_all[:, -1] - k, 0.0) * (
                ones - barrier_crossed_from_below
            )
        elif opt_type == FinFXBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_below
        elif opt_type == FinFXBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * (
                ones - barrier_crossed_from_below
            )
        elif opt_type == FinFXBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * (
                ones - barrier_crossed_from_above
            )
        elif opt_type == FinFXBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_above
        else:
            raise FinError("Unknown barrier option type." + str(self.opt_type))

        v = payoff.mean() * exp(-r_d * t)

        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_level)
        s += label_to_string("NUM OBSERVATIONS", self.num_obs_per_year)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("NOTIONAL CURRENCY", self.notional_currency, "")
        return s

    ###########################################################################

    def _print(self):
        """Print a lis_t of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        print(self)


########################################################################################
