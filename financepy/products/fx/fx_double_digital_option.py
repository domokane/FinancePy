##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.math import n_vect  # n_prime_vect

from ...utils.global_vars import g_days_in_year
from ...utils.error import FinError

# from ...products.equity.EquityOption import FinOption
from ...utils.date import Date

# from ...products.fx.FinFXModelTypes import FinFXModel
from ...models.black_scholes import BlackScholes
from ...utils.helpers import check_argument_types


class FXDoubleDigitalOption:

    def __init__(
        self,
        expiry_dt: Date,
        upper_strike: (float, np.ndarray),
        lower_strike: (float, np.ndarray),
        currency_pair: str,  # FORDOM
        notional: float,
        prem_currency: str,
        spot_days: int = 0,
    ):
        """Create the FX Double Digital Option object. Inputs include
        expiry date, upper strike, lower strike, currency pair,
        option type notional and the currency of the notional.
        An adjustment for spot days is enabled. All currency rates
        must be entered in the price in domestic currency of one unit
        of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency."""

        check_argument_types(self.__init__, locals())

        delivery_dt = expiry_dt.add_weekdays(spot_days)

        if delivery_dt < expiry_dt:
            raise FinError("Delivery date must be on or after expiry date.")

        self.expiry_dt = expiry_dt
        self.delivery_dt = delivery_dt

        if np.any(upper_strike < 0.0):
            raise FinError("Negative upper strike.")

        if np.any(lower_strike < 0.0):
            raise FinError("Negative lower strike.")

        self.upper_strike = upper_strike
        self.lower_strike = lower_strike

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.currency_pair = currency_pair
        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        if prem_currency != self.dom_name and prem_currency != self.for_name:
            raise FinError("Notional currency not in currency pair.")

        self.prem_currency = prem_currency

        self.notional = notional

        self.spot_days = spot_days

    ###########################################################################

    def value(
        self,
        value_dt,
        spot_fx_rate,  # 1 unit of foreign in domestic
        domestic_curve,
        foreign_curve,
        model,
    ):
        """Valuation of a double digital option using Black-Scholes model.
        The option pays out the notional in the premium currency if the
        fx rate is between the upper and lower strike at maturity. The
        valuation is equivalent to the valuation of the difference of
        the value of two digital puts, one with the upper and the other
        with the lower strike"""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as valuation date"
            )

        if foreign_curve.value_dt != value_dt:
            raise FinError(
                "Foreign Curve valuation date not same as valuation date"
            )

        if isinstance(value_dt, Date):
            spot_dt = value_dt.add_weekdays(self.spot_days)
            t_del = (self.delivery_dt - spot_dt) / g_days_in_year
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
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

        r_d = -np.log(dom_df) / t_del
        r_f = -np.log(for_df) / t_del

        s0 = spot_fx_rate
        K1 = self.lower_strike
        K2 = self.upper_strike

        if isinstance(model, BlackScholes):

            volatility = model.volatility
            ln_s0_k1 = np.log(s0 / K1)
            ln_s0_k2 = np.log(s0 / K2)
            den = volatility * np.sqrt(t_exp)
            v2 = volatility * volatility
            mu = r_d - r_f
            lower_d2 = (ln_s0_k1 + (mu - v2 / 2.0) * t_del) / den
            upper_d2 = (ln_s0_k2 + (mu - v2 / 2.0) * t_del) / den

            if self.prem_currency == self.for_name:
                lower_digital = s0 * np.exp(-r_f * t_del) * n_vect(-lower_d2)
                upper_digital = s0 * np.exp(-r_f * t_del) * n_vect(-upper_d2)
            elif self.prem_currency == self.dom_name:
                lower_digital = np.exp(-r_f * t_del) * n_vect(-lower_d2)
                upper_digital = np.exp(-r_f * t_del) * n_vect(-upper_d2)

            v = (upper_digital - lower_digital) * self.notional

        return v


###############################################################################
