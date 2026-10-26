##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
from numba import njit
from numpy.typing import ArrayLike

from ...utils.date import Date
from ...utils.math import normcdf_prime
from ...utils.global_vars import G_SMALL
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.global_types import FXDeltaMethodTypes
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_currency_pair

from ...models.equity_crr_tree import crr_tree_val_avg
from ...models.sabr import vol_function_sabr
from ...models.sabr import SABR
from ...models.black_scholes import BlackScholes
from ...models.black import Black

from ...models.bs_european_analytic import european_value, delta

from ...utils.helpers import check_argument_types, label_to_string

from ...utils.math import normcdf

from ...models.model import Model
from ...market.curves.discount_curve import DiscountCurve

from ...utils.helpers import option_years
from ...utils.fx_conventions import fx_spot_date

########################################################################################
# TODO: Refactor code to use FinBlackScholesAnalytic
########################################################################################


def f(volatility: float, *args):
    """This is the objective function used in the determination of the FX
    Option implied volatility which is computed in the class below."""

    self = args[0]
    value_dt = args[1]
    spot_fx_rate = args[2]
    domestic_curve = args[3]
    foreign_curve = args[4]
    price = args[5]

    model = BlackScholes(volatility)

    vdf = self.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)["v"]

    obj_fn = vdf - price

    return obj_fn


########################################################################################


def fvega(volatility: float, *args):
    """This is the derivative of the objective function with respect to the
    option volatility. It is used to speed up the determination of the FX
    Option implied volatility which is computed in the class below."""

    self = args[0]
    value_dt = args[1]
    spot_fx_rate = args[2]
    domestic_curve = args[3]
    foreign_curve = args[4]
#    price = args[5]
    currency_pair = args[6]

    model = BlackScholes(volatility)

    fprime = self.vega(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)

    return fprime


########################################################################################


@njit(fastmath=True, cache=True)
def fast_delta(s, t, k, r_d, r_f, vol, delta_type_value, opt_type_value):
    """Calculation of the FX Option delta. Used in the determination of
    the volatility sur_face. Avoids discount curve interpolation so it
    should be slightly faster than the full calculation of delta."""

    pips_spot_delta = delta(s, t, k, r_d, r_f, vol, opt_type_value)

    if delta_type_value == FXDeltaMethodTypes.SPOT_DELTA.value:
        return pips_spot_delta

    elif delta_type_value == FXDeltaMethodTypes.FORWARD_DELTA.value:

        pips_fwd_delta = pips_spot_delta * np.exp(r_f * t)
        return pips_fwd_delta

    elif delta_type_value == FXDeltaMethodTypes.SPOT_DELTA_PREM_ADJ.value:

        vpctf = european_value(s, t, k, r_d, r_f, vol, opt_type_value) / s
        pct_spot_delta_prem_adj = pips_spot_delta - vpctf
        return pct_spot_delta_prem_adj

    elif delta_type_value == FXDeltaMethodTypes.FORWARD_DELTA_PREM_ADJ.value:

        vpctf = european_value(s, t, k, r_d, r_f, vol, opt_type_value) / s
        pct_fwd_delta_prem_adj = np.exp(r_f * t) * (pips_spot_delta - vpctf)
        return pct_fwd_delta_prem_adj

    else:
        raise FinError("Unknown FinFXDeltaMethod")


########################################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EURUSD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
########################################################################################


class FXVanillaOption:
    """This is a class for an FX Option trade. It permits the user to
    calculate the price of an FX Option trade which can be expressed in a
    number of ways depending on the investor or hedger's currency. It aslo
    allows the calculation of the option's delta in a number of forms as
    well as the various Greek risk sensitivies."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_fx_rate: float,  # DOM per unit of FOR
        currency_pair: str,
        opt_type: OptionTypes,
        notional: float,
        prem_currency: str
    ) -> None:
        """Create the FX Vanilla Option object. Inputs include expiry date,
        strike, currency pair, OPTION_TYPE (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency."""

        check_argument_types(self.__init__, locals())

        check_currency_pair(currency_pair)

        self.currency_pair = currency_pair

        """ The FX rate the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1)"""

        self.expiry_dt = expiry_dt

        # Need to check if this is correct
        self.delivery_dt = fx_spot_date(currency_pair, expiry_dt)

        self.strike_fx_rate = strike_fx_rate

        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        if prem_currency != self.dom_name and prem_currency != self.for_name:
            raise FinError("Premium currency not in currency pair.")

        self.prem_currency = prem_currency
        self.notional = notional

        if (
            opt_type != OptionTypes.EUROPEAN_CALL
            and opt_type != OptionTypes.EUROPEAN_PUT
            and opt_type != OptionTypes.AMERICAN_CALL
            and opt_type != OptionTypes.AMERICAN_PUT
        ):
            raise FinError("Unknown OPTION_TYPE:" + str(opt_type))

        self.opt_type = opt_type

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: float | np.ndarray,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        forward_fx_rate: float | None = None,
    ):
        """Value the option in domestic currency.

        BlackScholes and SABR use spot FX. A Black model prices European
        options using Black-76 on an FX forward. For Black, the forward is
        inferred from the spot and discount curves unless forward_fx_rate
        is supplied explicitly (domestic currency per foreign unit).
        """

        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        spot_dt = fx_spot_date(self.currency_pair, value_dt)

        # Volatility horizon: valuation date to option expiry
        t_exp = option_years(value_dt, self.expiry_dt)

        s0 = spot_fx_rate
        k = self.strike_fx_rate

        # Discount factors to spot and delivery dates
        df_d_spot = domestic_curve.df(spot_dt)
        df_f_spot = foreign_curve.df(spot_dt)

        df_d_del = domestic_curve.df(self.delivery_dt)
        df_f_del = foreign_curve.df(self.delivery_dt)

        r_d = -np.log(df_d_del / df_d_spot) / t_exp
        r_f = -np.log(df_f_del / df_f_spot) / t_exp

        # Forward FX rate for the delivery date
        f0t = s0 * (df_f_del / df_f_spot) / (df_d_del / df_d_spot)

        vdf = None

        if isinstance(model, Black):
            if self.opt_type not in (
                OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT
            ):
                raise FinError("Black-76 FX pricing supports European options only.")

            domestic_df = df_d_del

            forward = (
                f0t if forward_fx_rate is None else forward_fx_rate
            )

            if np.any(np.asarray(forward) <= 0.0):
                raise FinError("FX forward rate must be positive.")
            vdf = model.value(forward, k, t_exp, domestic_df, self.opt_type)

        elif isinstance(model, BlackScholes) or isinstance(model, SABR):
            if forward_fx_rate is not None:
                raise FinError("Explicit FX forward requires a Black model.")

            if isinstance(model, BlackScholes):
                volatility = model.volatility
            elif isinstance(model, SABR):

                params_list = np.array([model.alpha, model.beta, model.rho, model.nu])

                volatility = vol_function_sabr(params_list, f0t, k, t_exp)

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            v = np.maximum(volatility, 1e-10)

            if self.opt_type == OptionTypes.EUROPEAN_CALL:

                vdf = df_d_del * european_value(
                    f0t, t_exp, k, 0.0, 0.0, v,
                    OptionTypes.EUROPEAN_CALL.value
                )

            elif self.opt_type == OptionTypes.EUROPEAN_PUT:

                vdf = df_d_del * european_value(
                    f0t, t_exp, k, 0.0, 0.0, v,
                    OptionTypes.EUROPEAN_PUT.value
                )

            elif self.opt_type == OptionTypes.AMERICAN_CALL:

                num_steps_per_year = 100

                vdf = df_d_spot * crr_tree_val_avg(
                    s0,
                    r_d,
                    r_f,
                    volatility,
                    num_steps_per_year,
                    t_exp,
                    OptionTypes.AMERICAN_CALL.value,
                    k,
                    True
                )["value"]

            elif self.opt_type == OptionTypes.AMERICAN_PUT:

                num_steps_per_year = 100

                vdf = df_d_spot * crr_tree_val_avg(
                    s0,
                    r_d,
                    r_f,
                    volatility,
                    num_steps_per_year,
                    t_exp,
                    OptionTypes.AMERICAN_PUT.value,
                    k,
                    True
                )["value"]
            else:
                raise FinError("Unknown OPTION_TYPE")

        # The option value v is in domestic currency terms but the value of
        # the option may be quoted in either currency terms and so we calculate
        # these

        if self.prem_currency == self.dom_name:
            notional_dom = self.notional
            notional_for = self.notional / self.strike_fx_rate
        elif self.prem_currency == self.for_name:
            notional_dom = self.notional * self.strike_fx_rate
            notional_for = self.notional
        else:
            raise FinError("Invalid notional currency.")

        s0 = np.asarray(s0, dtype=float)

        pips_dom = vdf
        pips_for = vdf / (s0 * self.strike_fx_rate)

        cash_dom = vdf * notional_dom / self.strike_fx_rate
        cash_for = vdf * notional_for / s0

        pct_dom = vdf / self.strike_fx_rate
        pct_for = vdf / s0

        return {
            "v": vdf,
            "cash_dom": cash_dom,
            "cash_for": cash_for,
            "pips_dom": pips_dom,
            "pips_for": pips_for,
            "pct_dom": pct_dom,
            "pct_for": pct_for,
            "not_dom": notional_dom,
            "not_for": notional_for,
            "ccy_dom": self.dom_name,
            "ccy_for": self.for_name,
        }

    ###########################################################################

    def delta_bump(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        ccy1_discount_curve: DiscountCurve,
        ccy2_discount_curve: DiscountCurve,
        model: Model,
    ):
        """Calculation of the FX option delta by bumping the spot FX rate by
        1 cent of its value. This gives the FX spot delta. For speed we prefer
        to use the analytical calculation of the derivative given below."""

        check_curve_dt(value_dt, ccy1_discount_curve, ccy2_discount_curve)

        bump_size = 0.0001 * spot_fx_rate

        v = self.value(value_dt,
                       spot_fx_rate,
                       ccy1_discount_curve,
                       ccy2_discount_curve,
                       model)

        v_bumped = self.value(
            value_dt,
            spot_fx_rate + bump_size,
            ccy1_discount_curve,
            ccy2_discount_curve,
            model,
        )

        if isinstance(v_bumped, dict):
            delta = (v_bumped["v"] - v["v"]) / bump_size
        else:
            delta = (v_bumped - v) / bump_size

        return delta

    ###########################################################################

    def delta(self, value_dt: Date,
              spot_fx_rate: float,
              domestic_curve: DiscountCurve,
              foreign_curve: DiscountCurve,
              model: Model):
        """Calculation of the FX Option delta. There are several definitions
        of delta and so we are required to return a dictionary of values. The
        definitions can be found on Page 44 of Foreign Exchange Option Pricing
        by Iain Clark, published by Wiley Finance."""

        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        s0 = spot_fx_rate
        t_exp = option_years(value_dt, self.expiry_dt)

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        k = self.strike_fx_rate

        pips_spot_delta = None
        pips_fwd_delta = None
        pct_spot_delta_prem_adj = None
        pct_fwd_delta_prem_adj = None

        if isinstance(model, BlackScholes):

            v = model.volatility

            if np.any(v < 0.0):
                raise FinError("Volatility should not be negative.")

            v = np.maximum(v, G_SMALL)

            pips_spot_delta = delta(s0, t_exp, k, r_d, r_f, v, self.opt_type.value)
            pips_fwd_delta = pips_spot_delta * np.exp(r_f * t_exp)
            vpctf = european_value(s0, t_exp, k, r_d, r_f, v, self.opt_type.value) / s0
            pct_spot_delta_prem_adj = pips_spot_delta - vpctf
            pct_fwd_delta_prem_adj = np.exp(r_f * t_exp) * (pips_spot_delta - vpctf)

        else:
            raise FinError("Unknown Model Type")

        return {
            "pips_spot_delta": pips_spot_delta,
            "pips_fwd_delta": pips_fwd_delta,
            "pct_spot_delta_prem_adj": pct_spot_delta_prem_adj,
            "pct_fwd_delta_prem_adj": pct_fwd_delta_prem_adj,
        }

    ###########################################################################

    def fast_delta(self, t, s, r_d, r_f, vol):
        """Calculation of the FX Option delta. Used in the determination of
        the volatility sur_face. Avoids discount curve interpolation so it
        should be slightly faster than the full calculation of delta."""

        k = self.strike_fx_rate

        pips_spot_delta = delta(s, t, k, r_d, r_f, vol, self.opt_type.value)
        pips_fwd_delta = pips_spot_delta * np.exp(r_f * t)

        vpctf = european_value(s, t, k, r_d, r_f, vol, self.opt_type.value) / s

        pct_spot_delta_prem_adj = pips_spot_delta - vpctf
        pct_fwd_delta_prem_adj = np.exp(r_f * t) * (pips_spot_delta - vpctf)

        return {
            "pips_spot_delta": pips_spot_delta,
            "pips_fwd_delta": pips_fwd_delta,
            "pct_spot_delta_prem_adj": pct_spot_delta_prem_adj,
            "pct_fwd_delta_prem_adj": pct_fwd_delta_prem_adj,
        }

    ###########################################################################

    def gamma(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """This function calculates the FX Option Gamma using spot delta."""
        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        t_exp = option_years(value_dt, self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        k = self.strike_fx_rate
        s0 = spot_fx_rate

        if isinstance(model, BlackScholes):

            volatility = model.volatility

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            ln_s0_k = np.log(s0 / k)
            sqrt_t = np.sqrt(t_exp)
            den = volatility * sqrt_t
            mu = r_d - r_f
            v2 = volatility * volatility
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t_exp) / den
            gamma = np.exp(-r_f * t_exp) * normcdf_prime(d1)
            gamma = gamma / s0 / den
        else:
            raise FinError("Unknown Model Type")

        return gamma

    ###########################################################################

    def vega(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """This function calculates the FX Option Vega using the spot delta."""

        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        s0 = spot_fx_rate

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_fx_rate
        t_exp = option_years(value_dt, self.expiry_dt)

        if isinstance(model, BlackScholes):

            volatility = model.volatility

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            ln_s0_k = np.log(s0 / k)
            sqrt_t = np.sqrt(t_exp)
            den = volatility * sqrt_t
            mu = r_d - r_f
            v2 = volatility * volatility
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t_exp) / den
            vega = s0 * sqrt_t * np.exp(-r_f * t_exp) * normcdf_prime(d1)
        else:
            raise FinError("Unknown Model type")

        return vega

    ###########################################################################

    def theta(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """This function calculates the time decay of the FX option."""

        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        s0 = spot_fx_rate

        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)
        k = self.strike_fx_rate
        t_exp = option_years(value_dt, self.expiry_dt)

        if isinstance(model, BlackScholes):

            vol = model.volatility

            if np.any(vol < 0.0):
                raise FinError("Volatility should not be negative.")

            vol = np.maximum(vol, 1e-10)

            ln_s0_k = np.log(s0 / k)
            sqrt_t = np.sqrt(t_exp)
            den = vol * sqrt_t
            mu = r_d - r_f
            v2 = vol * vol
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t_exp) / den
            d2 = (ln_s0_k + (mu - v2 / 2.0) * t_exp) / den

            if self.opt_type == OptionTypes.EUROPEAN_CALL:
                v = -s0 * np.exp(-r_f * t_exp) * normcdf_prime(d1) * vol / 2.0 / sqrt_t
                v = v + r_f * s0 * np.exp(-r_f * t_exp) * normcdf(d1)
                v = v - r_d * k * np.exp(-r_d * t_exp) * normcdf(d2)
            elif self.opt_type == OptionTypes.EUROPEAN_PUT:
                v = -s0 * np.exp(-r_f * t_exp) * normcdf_prime(d1) * vol / 2.0 / sqrt_t
                v = v + r_d * k * np.exp(-r_d * t_exp) * normcdf(-d2)
                v = v - r_f * s0 * np.exp(-r_f * t_exp) * normcdf(-d1)
            else:
                raise FinError("Unknown OPTION_TYPE")

        else:
            raise FinError("Unknown Model Type")

        return v

    ###########################################################################

    def implied_volatility(self,
                           value_dt: Date,
                           spot_fx_rate: float,
                           domestic_curve: DiscountCurve,
                           foreign_curve: DiscountCurve,
                           model: Model,
                           price: float):
        """This function determines the implied volatility of an FX option
        given a price and the other option details. It uses a one-dimensional
        Newton root search algorith to determine the implied volatility."""

        s0 = spot_fx_rate

        argtuple = (
            self,
            value_dt,
            s0,
            domestic_curve,
            foreign_curve,
            price,
            self.currency_pair
        )

        sigma = optimize.newton(
            f,
            x0=0.2,
            fprime=fvega,
            args=argtuple,
            tol=1e-6,
            maxiter=50,
            fprime2=None,
        )
        return sigma

    ###########################################################################

    def value_mc(self,
                 value_dt: Date,
                 spot_fx_rate: float,  # DOM per unit of FOR
                 domestic_curve: DiscountCurve,
                 foreign_curve: DiscountCurve,
                 model: Model,
                 num_paths: int = 10000,
                 seed=4242,
                 ):
        """Calculate the value of an FX Option using Monte Carlo methods.
        This function can be used to validate the risk measures calculated
        above or used as the starting code for a model exotic FX product that
        cannot be priced analytically. This function uses Numpy vectorisation
        for speed of execution."""

        check_curve_dt(value_dt, domestic_curve, foreign_curve)

        spot_dt = fx_spot_date(self.currency_pair, value_dt)

        df_d_spot = domestic_curve.df(spot_dt)
        df_f_spot = foreign_curve.df(spot_dt)
        df_d_del = domestic_curve.df(self.delivery_dt)
        df_f_del = foreign_curve.df(self.delivery_dt)

        f0t = spot_fx_rate * (df_f_del / df_f_spot) / (df_d_del / df_d_spot)

        k = self.strike_fx_rate
        t_exp = option_years(value_dt, self.expiry_dt)

        if isinstance(model, BlackScholes):
            volatility = model.volatility
        else:
            raise FinError("Model Type invalid")

        np.random.seed(seed)

        mu = np.log(f0t / spot_fx_rate) / t_exp

        v2 = volatility**2
        k = self.strike_fx_rate
        sqrt_dt = np.sqrt(t_exp)
        s0 = spot_fx_rate

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, num_paths))
        s = s0 * np.exp((mu - v2 / 2.0) * t_exp)
        m = np.exp(g * sqrt_dt * volatility)
        s_1 = s * m
        s_2 = s / m

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            payoff_a_1 = np.maximum(s_1 - k, 0.0)
            payoff_a_2 = np.maximum(s_2 - k, 0.0)
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            payoff_a_1 = np.maximum(k - s_1, 0.0)
            payoff_a_2 = np.maximum(k - s_2, 0.0)
        else:
            raise FinError("Unknown OPTION_TYPE.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * df_d_del / 2.0
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("DELIVERY DATE", self.delivery_dt)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("PREMIUM CCY", self.prem_currency)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
