##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import List, Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import check_argument_types, label_to_string, option_years
from ...utils.check_values import (
    check_curve_dt,
    check_stock_price,
    check_strike_price,
    check_shapes,
)
from ...utils.mc_result import MCResult
from ...market.curves.discount_curve import DiscountCurve
from ...models.model import Model
from ...models.black_scholes import BlackScholes
from ...models.bs_power_option import power_option_value, power_option_delta
from .fx_option import FXOption

########################################################################################


class FXPowerOption(FXOption):
    """Class for managing European FX power options (calls and puts).

    Under Black-Scholes / Garman-Kohlhagen dynamics, the exchange rate S_t
    (domestic currency per unit of foreign currency) follows a geometric
    Brownian motion with domestic risk-free rate r_d and foreign risk-free rate r_f:

        dS_t = (r_d - r_f) S_t dt + sigma S_t dW_t

    At expiry T, the payoff in domestic currency is:
        Call: max(S_T**n - K, 0) * notional
        Put:  max(K - S_T**n, 0) * notional

    where n > 0 is the power parameter and K is the strike price.
    """

    def __init__(
        self,
        expiry_dt: Union[Date, List[Date]],
        strike_fx_rate: Union[float, np.ndarray],
        opt_type: OptionTypes,
        power: float = 1.0,
        notional: float = 1.0,
        currency_pair: Union[str, None] = None,
    ) -> None:
        """Create European FX power option."""
        check_argument_types(self.__init__, locals())

        check_strike_price(strike_fx_rate)

        if not np.isfinite(power) or power <= 0.0:
            raise FinError("Power must be strictly positive.")

        if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            raise FinError(
                f"Invalid option type: {opt_type}. Must be EUROPEAN_CALL or EUROPEAN_PUT."
            )

        if not np.isfinite(notional) or notional <= 0.0:
            raise FinError("Notional must be strictly positive.")

        if currency_pair is not None and len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.expiry_dt = expiry_dt
        self.strike_fx_rate = strike_fx_rate
        self.opt_type = opt_type
        self.opt_type_value = opt_type.value
        self.power = float(power)
        self.notional = float(notional)
        self.currency_pair = currency_pair

    ###########################################################################

    def intrinsic(
        self,
        value_dt: Date,
        spot_fx_rate: Union[float, np.ndarray],
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
    ):
        """Calculate the discounted intrinsic value of the FX power option."""
        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(spot_fx_rate)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_fx_rate
        s0 = spot_fx_rate
        n = self.power

        check_shapes(s0, k, t_exp)

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            phi = 1.0
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            phi = -1.0
        else:
            raise FinError("Unknown option type")

        s_t = s0 * np.exp((r_d - r_f) * t_exp)

        intrinsic_value = (
            np.exp(-r_d * t_exp)
            * np.maximum(phi * (s_t**n - k), 0.0)
        )

        return intrinsic_value * self.notional

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: Union[float, np.ndarray],
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """Value the European FX power option using the Black-Scholes / Garman-Kohlhagen model."""
        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(spot_fx_rate)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_fx_rate
        s0 = spot_fx_rate
        n = self.power

        check_shapes(s0, k, t_exp)

        if isinstance(model, BlackScholes):
            v = model.volatility
            if isinstance(s0, (list, np.ndarray)):
                s0_arr = np.asarray(s0, dtype=float)
                values = np.empty_like(s0_arr)
                for idx, s_val in np.ndenumerate(s0_arr):
                    values[idx] = power_option_value(
                        float(s_val), t_exp, k, r_d, r_f, v, self.opt_type_value, n
                    )
                return values * self.notional
            else:
                value = power_option_value(
                    float(s0), t_exp, k, r_d, r_f, v, self.opt_type_value, n
                )
                return value * self.notional
        else:
            raise FinError("Unknown Model Type")

    ###########################################################################

    def delta(
        self,
        value_dt: Date,
        spot_fx_rate: Union[float, np.ndarray],
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """Analytical spot delta for a European FX power call or put."""
        if not isinstance(model, BlackScholes):
            return super().delta(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(spot_fx_rate)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)
        check_shapes(spot_fx_rate, self.strike_fx_rate, t_exp)

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        if isinstance(spot_fx_rate, (list, np.ndarray)):
            s_arr = np.asarray(spot_fx_rate, dtype=float)
            deltas = np.empty_like(s_arr)
            for idx, s_val in np.ndenumerate(s_arr):
                deltas[idx] = power_option_delta(
                    float(s_val),
                    t_exp,
                    self.strike_fx_rate,
                    r_d,
                    r_f,
                    model.volatility,
                    self.opt_type_value,
                    self.power,
                )
            return deltas * self.notional
        else:
            d = power_option_delta(
                float(spot_fx_rate),
                t_exp,
                self.strike_fx_rate,
                r_d,
                r_f,
                model.volatility,
                self.opt_type_value,
                self.power,
            )
            return d * self.notional

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
    ) -> MCResult:
        """Value the FX power option using Monte Carlo simulation under Garman-Kohlhagen dynamics.

        Parameters
        ----------
        value_dt : Date
            Valuation date.
        spot_fx_rate : float
            Current spot exchange rate (DOM per FOR).
        domestic_curve : DiscountCurve
            Domestic risk-free discount curve.
        foreign_curve : DiscountCurve
            Foreign risk-free discount curve.
        model : Model
            Black-Scholes model.
        num_paths : int
            Number of Monte Carlo paths.
        seed : int
            Random number seed.

        Returns
        -------
        MCResult
            Monte Carlo value and standard error.
        """
        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(spot_fx_rate)
        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        if not isinstance(model, BlackScholes):
            raise FinError("Unknown Model Type")

        if self.power <= 0.0:
            raise FinError("Power must be positive")

        if num_paths <= 0:
            raise FinError("Number of paths must be positive")

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        s0 = float(spot_fx_rate)
        k = float(self.strike_fx_rate)
        n = self.power
        v = model.volatility

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            phi = 1.0
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            phi = -1.0
        else:
            raise FinError("Unknown option type")

        rng = np.random.default_rng(seed)
        num_pairs = (num_paths + 1) // 2
        z = rng.standard_normal(num_pairs)
        z = np.concatenate((z, -z))[:num_paths]

        mu = r_d - r_f - 0.5 * v * v
        s_t = s0 * np.exp(mu * t_exp + v * np.sqrt(t_exp) * z)

        payoff = np.maximum(phi * (s_t**n - k), 0.0)
        discounted_payoff = np.exp(-r_d * t_exp) * payoff * self.notional

        value = float(np.mean(discounted_payoff))
        std_error = float(np.std(discounted_payoff, ddof=1) / np.sqrt(num_paths))

        return MCResult(value, std_error)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("POWER", self.power)
        s += label_to_string("NOTIONAL", self.notional)
        if self.currency_pair is not None:
            s += label_to_string("CURRENCY PAIR", self.currency_pair, "")
        else:
            s = s.rstrip("\n")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
