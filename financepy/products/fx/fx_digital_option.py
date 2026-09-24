##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np


from ...utils.math import normcdf_vect  # normcdf_prime_vect

from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.error import FinError

from ...utils.date import Date

from ...utils.helpers import check_argument_types
from ...utils.global_types import OptionTypes
from ...utils.check_values import check_curve_dt

from ...models.black_scholes import BlackScholes
from ...models.model import Model
from ...market.curves.discount_curve import DiscountCurve

########################################################################################


class FXDigitalOption:
    """FX Digital Option"""

    def __init__(
        self,
        expiry_dt: Date,
        strike_fx_rate: Union[float, np.ndarray],
        currency_pair: str,  # FORDOM
        opt_type: Union[OptionTypes, list],
        notional: float,
        prem_currency: str,
        spot_days: int = 0,
    ) -> None:
        """Create the FX Digital Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency."""

        check_argument_types(self.__init__, locals())

        delivery_dt = expiry_dt.add_weekdays(spot_days)

        if delivery_dt < expiry_dt:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.expiry_dt = expiry_dt
        self.delivery_dt = delivery_dt

        if np.any(strike_fx_rate <= 0.0):
            raise FinError("Strike must be greater than zero.")

        self.strike_fx_rate = strike_fx_rate

        self.currency_pair = currency_pair
        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        if prem_currency not in [self.dom_name, self.for_name]:
            raise FinError("Notional currency not in currency pair.")

        self.prem_currency = prem_currency

        self.notional = notional

        if opt_type not in [OptionTypes.DIGITAL_CALL, OptionTypes.DIGITAL_PUT]:
            raise FinError("Unknown Digital Option Type:" + str(opt_type))

        self.opt_type = opt_type
        self.spot_days = spot_days

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: float,  # 1 unit of foreign in domestic
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """Valuation of a digital option using Black-Scholes model. This
        allows for 4 cases - first upper barriers that when crossed pay out
        cash (calls) and lower barriers than when crossed from above cause a
        cash payout (puts) PLUS the fact that the cash payment can be in
        domestic or foreign currency."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        if isinstance(value_dt, Date):
            spot_dt = value_dt.add_weekdays(self.spot_days)
            t_del = (self.delivery_dt - spot_dt) / G_DAYS_IN_YEAR
            t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEAR
        else:
            t_del = value_dt
            t_exp = t_del

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("spot_fx_rate must be greater than zero.")

        if np.any(t_del < 0.0):
            raise FinError("Option time to maturity is less than zero.")

        t_del = np.maximum(t_del, 1e-10)

        # TODO RESOLVE t_del versus TEXP
        dom_df = domestic_curve.df_t(t_del)
        for_df = foreign_curve.df_t(t_del)

        if not isinstance(model, BlackScholes):
            raise FinError("Model must be BlackScholes.")

        f = spot_fx_rate * for_df / dom_df
        k = self.strike_fx_rate
        v = model.volatility
        vol_sqrt_t = v * np.sqrt(t_exp)

        d1 = (np.log(f / k) + 0.5 * (v**2) * t_exp) / vol_sqrt_t
        d2 = d1 - vol_sqrt_t

        if self.opt_type == OptionTypes.DIGITAL_CALL and self.for_name == self.prem_currency:
            v = spot_fx_rate * for_df * normcdf_vect(d1)

        elif self.opt_type == OptionTypes.DIGITAL_PUT and self.for_name == self.prem_currency:
            v = spot_fx_rate * for_df * normcdf_vect(-d1)

        elif self.opt_type == OptionTypes.DIGITAL_CALL and self.dom_name == self.prem_currency:
            v = dom_df * normcdf_vect(d2)

        elif self.opt_type == OptionTypes.DIGITAL_PUT and self.dom_name == self.prem_currency:
            v = dom_df * normcdf_vect(-d2)

        else:
            raise FinError("Unknown option type")

        v = v * self.notional

        return v


########################################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
    ):
        """Value the FX digital option using Monte Carlo simulation."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        if spot_fx_rate <= 0.0:
            raise FinError("spot_fx_rate must be greater than zero.")

        if num_paths <= 0:
            raise FinError("num_paths must be positive")

        if not isinstance(model, BlackScholes):
            raise FinError("Model must be BlackScholes.")

        # ------------------------------------------------------------------
        # Times
        # ------------------------------------------------------------------

        spot_dt = value_dt.add_weekdays(self.spot_days)

        t_del = (
            self.delivery_dt - spot_dt
        ) / G_DAYS_IN_YEAR

        t_exp = (
            self.expiry_dt - value_dt
        ) / G_DAYS_IN_YEAR

        if t_exp < 0.0:
            raise FinError("Option time to expiry is less than zero.")

        if t_del < 0.0:
            raise FinError("Option time to delivery is less than zero.")

        # ------------------------------------------------------------------
        # Market parameters
        # ------------------------------------------------------------------

        t = max(t_exp, 1.0e-10)

        dom_df = domestic_curve.df_t(t)
        for_df = foreign_curve.df_t(t)

        r_d = -np.log(dom_df) / t
        r_f = -np.log(for_df) / t

        mu = r_d - r_f

        v = model.volatility

        # ------------------------------------------------------------------
        # Simulate terminal FX rate
        # ------------------------------------------------------------------

        # Use Antithetic variables
        # A local RNG makes the Monte Carlo result reproducible for a given seed.
        rng = np.random.default_rng(seed)
        g = rng.standard_normal(num_paths)

        s = spot_fx_rate * np.exp((mu - v * v / 2.0) * t)
        m = np.exp(g * np.sqrt(t) * v)

        s_1 = s * m
        s_2 = s / m

        payoff_1 = None
        payoff_2 = None
        indicator_1 = None
        indicator_2 = None

        # ------------------------------------------------------------------
        # Digital event
        # ------------------------------------------------------------------

        if self.opt_type == OptionTypes.DIGITAL_CALL:
            indicator_1 = s_1 > self.strike_fx_rate
            indicator_2 = s_2 > self.strike_fx_rate

        elif self.opt_type == OptionTypes.DIGITAL_PUT:
            indicator_1 = s_1 < self.strike_fx_rate
            indicator_2 = s_2 < self.strike_fx_rate

        else:
            raise FinError(
                "Unknown Digital Option Type: "
                + str(self.opt_type)
            )

        # ------------------------------------------------------------------
        # Payoff
        # ------------------------------------------------------------------

        if self.prem_currency == self.dom_name:

            # Domestic cash-or-nothing digital:
            #
            #     payoff = N * 1(S_T > K)
            #
            # Analytic:
            #
            #     N * DF_dom * N(d2)

            payoff_1 = self.notional * indicator_1
            payoff_2 = self.notional * indicator_2

        elif self.prem_currency == self.for_name:

            # Foreign-currency digital:
            #
            #     payoff in foreign currency = N * 1(S_T > K)
            #
            # Converted into domestic currency at expiry:
            #
            #     payoff = N * S_T * 1(S_T > K)
            #
            # Analytic:
            #
            #     N * S_0 * DF_for * N(d1)

            payoff_1 = self.notional * s_1 * indicator_1
            payoff_2 = self.notional * s_2 * indicator_2

        else:
            raise FinError(
                "Notional currency not in currency pair."
            )

        # ------------------------------------------------------------------
        # Discount domestic-currency payoff
        # ------------------------------------------------------------------

        payoff = np.mean(payoff_1) + np.mean(payoff_2)
        value = dom_df * payoff / 2.0

        return value
