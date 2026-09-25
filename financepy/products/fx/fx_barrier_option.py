##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.global_types import BarrierTypes
from ...products.fx.fx_option import FXOption
from ...models.process_simulator import ProcessTypes
from ...models.process_simulator import GBMNumericalSchemeTypes
from ...models.barrier_option_mc import value_barrier_option_mc
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...models.barrier_option_model import barrier_option_value
from ...market.curves.discount_curve import DiscountCurve
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_strike_price
from ...utils.helpers import option_years
from ...models.model import Model

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

        values = barrier_option_value(
            self.strike_fx_rate,
            self.barrier_level,
            t_exp,
            spot_fx_rate,
            domestic_df,
            foreign_df,
            model.volatility,
            self.barrier_type.value,
            self.num_obs_per_year,
        )

        values = values * self.notional

        if np.isscalar(spot_fx_rate):
            return float(values)

        return np.asarray(values)

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float | np.ndarray,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_obs_per_year: int = 252,
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
