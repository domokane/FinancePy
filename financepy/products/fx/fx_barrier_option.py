##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.global_types import BarrierTypes
from ...products.fx.fx_option import FXOption
from ...models.process_simulator import ProcessTypes
from ...utils.global_types import GBMNumericalSchemeTypes
from ...models.barrier_option_mc import value_barrier_option_mc
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...models.fx_barrier_model import fx_barrier_value
from ...market.curves.discount_curve import DiscountCurve
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_strike_price
from ...utils.helpers import option_years

########################################################################################


class FXBarrierOption(FXOption):

    def __init__(
        self,
        expiry_dt: Date,
        strike_fx_rate: float,  # 1 unit of foreign in domestic
        currency_pair: str,  # FORDOM
        barrier_type: BarrierTypes,
        barrier_level: float,
        num_obs_per_year: int,
        notional_currency: str,
        notional: float = 1.0,
    ) -> None:
        """Create FX Barrier option product. This is an option that cancels if
        the FX rate crosses a barrier during the life of the option."""

        check_argument_types(self.__init__, locals())

        check_strike_price(strike_fx_rate)

        self.expiry_dt = expiry_dt
        self.strike_fx_rate = float(strike_fx_rate)
        self.currency_pair = currency_pair
        self.barrier_type = barrier_type
        self.barrier_level = float(barrier_level)
        self.num_obs_per_year = int(num_obs_per_year)
        self.notional = notional
        self.notional_currency = notional_currency

    ##########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: float | np.ndarray,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ) -> float | np.ndarray:
        """Value an FX barrier option for a scalar or array of spot rates."""

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date")

        t_exp = option_years(value_dt, self.expiry_dt)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)
        check_stock_price(spot_fx_rate)

        domestic_df = domestic_curve.df_t(t_exp)
        foreign_df = foreign_curve.df_t(t_exp)

        values = fx_barrier_value(
            spot_fx_rate,
            self.strike_fx_rate,
            self.barrier_level,
            t_exp,
            domestic_df,
            foreign_df,
            model.volatility,
            self.num_obs_per_year,
            self.barrier_type.value,
        )

        values = values * self.notional

        if isinstance(spot_fx_rate, float):
            return values
        else:
            return np.array(values)

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float | np.ndarray,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_obs_per_year: int=252,
        num_paths: int = 10000,
        seed: int = 42,
    ):
        """Value the FX Barrier Option using Monte Carlo."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        check_stock_price(spot_fx_rate)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        mu = r_d - r_f

        scheme = GBMNumericalSchemeTypes.ANTITHETIC

        model_params = (spot_fx_rate, mu, model.volatility, scheme)

        process_type = ProcessTypes.GBM_PROCESS

        value = value_barrier_option_mc(
            t_exp,
            self.strike_fx_rate,
            self.barrier_type,
            self.barrier_level,
            spot_fx_rate,
            r_d,
            process_type,
            model_params,
            num_obs_per_year,
            num_paths,
            seed,
        )

        value = value * self.notional

        return value

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("BARRIER TYPE", self.barrier_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_level)
        s += label_to_string("NUM OBSERVATIONS", self.num_obs_per_year)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("NOTIONAL CURRENCY", self.notional_currency, "")
        return s

    ###########################################################################

    # def value_mc(
    #     self,
    #     value_dt,
    #     spot_fx_rate,
    #     dom_interest_rate,
    #     process_type,
    #     model_params,
    #     num_ann_steps=552,
    #     num_paths=5000,
    #     seed=4242,
    # ):
    #     """Value the FX Barrier Option using Monte Carlo."""

    #     t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEAR
    #     num_time_steps = int(t * num_ann_steps)
    #     k = self.strike_fx_rate
    #     b = self.barrier_level
    #     s0 = spot_fx_rate
    #     barrier_type = self.barrier_type
    #     process = ProcessSimulator()

    #     r_d = dom_interest_rate

    #     #######################################################################

    #     if barrier_type == BarrierTypes.DOWN_AND_OUT_CALL and s0 <= b:
    #         return 0.0
    #     elif barrier_type == BarrierTypes.UP_AND_OUT_CALL and s0 >= b:
    #         return 0.0
    #     elif barrier_type == BarrierTypes.DOWN_AND_OUT_PUT and s0 <= b:
    #         return 0.0
    #     elif barrier_type == BarrierTypes.UP_AND_OUT_PUT and s0 >= b:
    #         return 0.0

    #     #######################################################################

    #     simple_call = False
    #     simple_put = False

    #     if barrier_type == BarrierTypes.DOWN_AND_IN_CALL and s0 <= b:
    #         simple_call = True
    #     elif barrier_type == BarrierTypes.UP_AND_IN_CALL and s0 >= b:
    #         simple_call = True
    #     elif barrier_type == BarrierTypes.UP_AND_IN_PUT and s0 >= b:
    #         simple_put = True
    #     elif barrier_type == BarrierTypes.DOWN_AND_IN_PUT and s0 <= b:
    #         simple_put = True

    #     if simple_put or simple_call:
    #         s_all = process.get_process(process_type, t, model_params, 1, num_paths, seed)

    #         if simple_call:
    #             s_t = s_all[:, -1]
    #             c = (np.maximum(s_t - k, 0.0)).mean()
    #             c = c * exp(-r_d * t)
    #             return c

    #         if simple_put:
    #             s_t = s_all[:, -1]
    #             p = (np.maximum(k - s_t, 0.0)).mean()
    #             p = p * exp(-r_d * t)
    #             return p

    #     # Otherwise get full set of paths
    #     s_all = process.get_process(process_type, t, model_params, num_time_steps, num_paths, seed)

    #     num_paths, num_time_steps = s_all.shape

    #     if barrier_type in (
    #         BarrierTypes.DOWN_AND_IN_CALL,
    #         BarrierTypes.DOWN_AND_OUT_CALL,
    #         BarrierTypes.DOWN_AND_IN_PUT,
    #         BarrierTypes.DOWN_AND_OUT_PUT,
    #     ):

    #         barrier_crossed_from_above = [False] * num_paths

    #         for p in nb.prange(num_paths):
    #             barrier_crossed_from_above[p] = np.any(s_all[p] <= b)

    #     if barrier_type in (
    #         BarrierTypes.UP_AND_IN_CALL,
    #         BarrierTypes.UP_AND_OUT_CALL,
    #         BarrierTypes.UP_AND_IN_PUT,
    #         BarrierTypes.UP_AND_OUT_PUT,
    #     ):

    #         barrier_crossed_from_below = [False] * num_paths
    #         for p in nb.prange(num_paths):
    #             barrier_crossed_from_below[p] = np.any(s_all[p] >= b)

    #     payoff = np.zeros(num_paths)
    #     ones = np.ones(num_paths)

    #     if barrier_type == BarrierTypes.DOWN_AND_OUT_CALL:
    #         payoff = np.maximum(s_all[:, -1] - k, 0.0) * (ones - barrier_crossed_from_above)
    #     elif barrier_type == BarrierTypes.DOWN_AND_IN_CALL:
    #         payoff = np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_above
    #     elif barrier_type == BarrierTypes.UP_AND_IN_CALL:
    #         payoff = np.maximum(s_all[:, -1] - k, 0.0) * barrier_crossed_from_below
    #     elif barrier_type == BarrierTypes.UP_AND_OUT_CALL:
    #         payoff = np.maximum(s_all[:, -1] - k, 0.0) * (ones - barrier_crossed_from_below)
    #     elif barrier_type == BarrierTypes.UP_AND_IN_PUT:
    #         payoff = np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_below
    #     elif barrier_type == BarrierTypes.UP_AND_OUT_PUT:
    #         payoff = np.maximum(k - s_all[:, -1], 0.0) * (ones - barrier_crossed_from_below)
    #     elif barrier_type == BarrierTypes.DOWN_AND_OUT_PUT:
    #         payoff = np.maximum(k - s_all[:, -1], 0.0) * (ones - barrier_crossed_from_above)
    #     elif barrier_type == BarrierTypes.DOWN_AND_IN_PUT:
    #         payoff = np.maximum(k - s_all[:, -1], 0.0) * barrier_crossed_from_above
    #     else:
    #         raise FinError("Unknown barrier option type." + str(self.barrier_type))

    #     v = payoff.mean() * exp(-r_d * t)

    #     return v

    ###########################################################################

    # def value_old(self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model):
    #     """Value FX Barrier Option using Black-Scholes model with closed-form
    #     analytical models."""

    #     # This prices the option using the formulae given in the paper
    #     # by Clewlow, Llanos and Strickland December 1994 which can be found at
    #     # https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf

    #     if isinstance(value_dt, Date) is False:
    #         raise FinError("Valuation date is not a Date")

    #     t_exp = option_years(value_dt, self.expiry_dt)
    #     check_curve_dt(value_dt, domestic_curve)
    #     check_curve_dt(value_dt, foreign_curve)
    #     check_stock_price(spot_fx_rate)

    #     if domestic_curve.anchor_dt != value_dt:
    #         raise FinError("Domestic Curve valuation date not same as option value date")

    #     if foreign_curve.anchor_dt != value_dt:
    #         raise FinError("Foreign Curve valuation date not same as option value date")

    #     k = self.strike_fx_rate
    #     s0 = spot_fx_rate
    #     h = self.barrier_level

    #     ln_s0_k = np.log(s0 / k)
    #     sqrt_t = np.sqrt(t_exp)

    #     dq = foreign_curve.df_t(t)
    #     df = domestic_curve.df_t(t)
    #     r_d = -np.log(df) / t
    #     rf = -np.log(dq) / t

    #     volatility = model.volatility
    #     sigma_root_t = volatility * sqrt_t
    #     v2 = volatility * volatility
    #     mu = r_d - rf
    #     d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / sigma_root_t
    #     d2 = (ln_s0_k + (mu - v2 / 2.0) * t) / sigma_root_t

    #     c = s0 * dq * normcdf(d1) - k * df * normcdf(d2)
    #     p = k * df * normcdf(-d2) - s0 * dq * normcdf(-d1)
    #     #        print("CALL:",c,"PUT:",p)

    #     if self.barrier_type == BarrierTypes.DOWN_AND_OUT_CALL and s0 <= h:
    #         return 0.0
    #     if self.barrier_type == BarrierTypes.UP_AND_OUT_CALL and s0 >= h:
    #         return 0.0
    #     if self.barrier_type == BarrierTypes.UP_AND_OUT_PUT and s0 >= h:
    #         return 0.0
    #     if self.barrier_type == BarrierTypes.DOWN_AND_OUT_PUT and s0 <= h:
    #         return 0.0
    #     if self.barrier_type == BarrierTypes.DOWN_AND_IN_CALL and s0 <= h:
    #         return c
    #     if self.barrier_type == BarrierTypes.UP_AND_IN_CALL and s0 >= h:
    #         return c
    #     if self.barrier_type == BarrierTypes.UP_AND_IN_PUT and s0 >= h:
    #         return p
    #     if self.barrier_type == BarrierTypes.DOWN_AND_IN_PUT and s0 <= h:
    #         return p

    #     num_observations = t * self.num_obs_per_year

    #     # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
    #     # Adjusts the barrier for discrete and not continuous observations
    #     h_adj = h
    #     if self.barrier_type == BarrierTypes.DOWN_AND_OUT_CALL:
    #         h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_IN_CALL:
    #         h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.UP_AND_IN_CALL:
    #         h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.UP_AND_OUT_CALL:
    #         h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.UP_AND_IN_PUT:
    #         h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.UP_AND_OUT_PUT:
    #         h_adj = h * exp(0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_OUT_PUT:
    #         h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_IN_PUT:
    #         h_adj = h * exp(-0.5826 * volatility * sqrt(t / num_observations))
    #     else:
    #         raise FinError("Unknown barrier option type." + str(self.barrier_type))

    #     h = h_adj

    #     if abs(volatility) < 1e-5:
    #         volatility = 1e-5

    #     ll = (mu + v2 / 2.0) / v2
    #     y = np.log(h * h / (s0 * k)) / sigma_root_t + ll * sigma_root_t
    #     x1 = np.log(s0 / h) / sigma_root_t + ll * sigma_root_t
    #     y1 = np.log(h / s0) / sigma_root_t + ll * sigma_root_t
    #     h_over_s = h / s0

    #     if self.barrier_type == BarrierTypes.DOWN_AND_OUT_CALL:
    #         if h >= k:
    #             c_do = (
    #                 s0 * dq * normcdf(x1)
    #                 - k * df * normcdf(x1 - sigma_root_t)
    #                 - s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
    #                 + k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(y1 - sigma_root_t)
    #             )
    #             price = c_do
    #         else:
    #             c_di = s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
    #                 h_over_s, 2.0 * ll - 2.0
    #             ) * normcdf(y - sigma_root_t)
    #             price = c - c_di
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_IN_CALL:
    #         if h <= k:
    #             c_di = s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
    #                 h_over_s, 2.0 * ll - 2.0
    #             ) * normcdf(y - sigma_root_t)
    #             price = c_di
    #         else:
    #             c_do = (
    #                 s0 * dq * normcdf(x1)
    #                 - k * df * normcdf(x1 - sigma_root_t)
    #                 - s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
    #                 + k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(y1 - sigma_root_t)
    #             )
    #             price = c - c_do
    #     elif self.barrier_type == BarrierTypes.UP_AND_IN_CALL:
    #         if h >= k:
    #             c_ui = (
    #                 s0 * dq * normcdf(x1)
    #                 - k * df * normcdf(x1 - sigma_root_t)
    #                 - s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
    #                 + k
    #                 * df
    #                 * pow(h_over_s, 2.0 * ll - 2.0)
    #                 * (normcdf(-y + sigma_root_t) - normcdf(-y1 + sigma_root_t))
    #             )
    #             price = c_ui
    #         else:
    #             price = c
    #     elif self.barrier_type == BarrierTypes.UP_AND_OUT_CALL:
    #         if h > k:
    #             c_ui = (
    #                 s0 * dq * normcdf(x1)
    #                 - k * df * normcdf(x1 - sigma_root_t)
    #                 - s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
    #                 + k
    #                 * df
    #                 * pow(h_over_s, 2.0 * ll - 2.0)
    #                 * (normcdf(-y + sigma_root_t) - normcdf(-y1 + sigma_root_t))
    #             )
    #             price = c - c_ui
    #         else:
    #             price = 0.0
    #     elif self.barrier_type == BarrierTypes.UP_AND_IN_PUT:
    #         if h > k:
    #             p_ui = -s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
    #                 h_over_s, 2.0 * ll - 2.0
    #             ) * normcdf(-y + sigma_root_t)
    #             price = p_ui
    #         else:
    #             p_uo = (
    #                 -s0 * dq * normcdf(-x1)
    #                 + k * df * normcdf(-x1 + sigma_root_t)
    #                 + s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
    #                 - k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(-y1 + sigma_root_t)
    #             )
    #             price = p - p_uo
    #     elif self.barrier_type == BarrierTypes.UP_AND_OUT_PUT:
    #         if h >= k:
    #             p_ui = -s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
    #                 h_over_s, 2.0 * ll - 2.0
    #             ) * normcdf(-y + sigma_root_t)
    #             price = p - p_ui
    #         else:
    #             p_uo = (
    #                 -s0 * dq * normcdf(-x1)
    #                 + k * df * normcdf(-x1 + sigma_root_t)
    #                 + s0 * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
    #                 - k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(-y1 + sigma_root_t)
    #             )
    #             price = p_uo
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_OUT_PUT:
    #         if h >= k:
    #             price = 0.0
    #         else:
    #             p_di = (
    #                 -s0 * dq * normcdf(-x1)
    #                 + k * df * normcdf(-x1 + sigma_root_t)
    #                 + s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
    #                 - k * df * pow(h_over_s, 2.0 * ll - 2.0) * (normcdf(y - sigma_root_t) - normcdf(y1 - sigma_root_t))
    #             )
    #             price = p - p_di
    #     elif self.barrier_type == BarrierTypes.DOWN_AND_IN_PUT:
    #         if h >= k:
    #             price = p
    #         else:
    #             p_di = (
    #                 -s0 * dq * normcdf(-x1)
    #                 + k * df * normcdf(-x1 + sigma_root_t)
    #                 + s0 * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
    #                 - k * df * pow(h_over_s, 2.0 * ll - 2.0) * (normcdf(y - sigma_root_t) - normcdf(y1 - sigma_root_t))
    #             )
    #             price = p_di
    #     else:
    #         raise FinError("Unknown barrier option type." + str(self.barrier_type))

    #     return price
