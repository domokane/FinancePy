##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
from numba import njit

from ...utils.date import Date
from ...utils.math import nprime
from ...utils.global_vars import g_days_in_year, g_small
from ...utils.error import FinError
from ...utils.global_types import OptionTypes

# from ...products.fx.FinFXModelTypes import FinFXModel
# from ...products.fx.FinFXModelTypes import FinFXModelBlackScholes
# from ...products.fx.FinFXModelTypes import FinFXModelSABR
from ...products.fx.fx_mkt_conventions import FinFXDeltaMethod

from ...models.equity_crr_tree import crr_tree_val_avg
from ...models.sabr import vol_function_sabr
from ...models.sabr import SABR
from ...models.black_scholes import BlackScholes

from ...models.black_scholes_analytic import bs_value, bs_delta

from ...utils.helpers import check_argument_types, label_to_string

from ...utils.math import N


###############################################################################
# TODO: Refactor code to use FinBlackScholesAnalytic
###############################################################################


def f(volatility, *args):
    """This is the objective function used in the determination of the FX
    Option implied volatility which is computed in the class below."""

    self = args[0]
    value_dt = args[1]
    spot_fx_rate = args[2]
    domestic_curve = args[3]
    foreign_curve = args[4]
    price = args[5]

    model = BlackScholes(volatility)

    vdf = self.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    obj_fn = vdf - price

    return obj_fn


###############################################################################


def fvega(volatility, *args):
    """This is the derivative of the objective function with respect to the
    option volatility. It is used to speed up the determination of the FX
    Option implied volatility which is computed in the class below."""

    self = args[0]
    value_dt = args[1]
    spot_fx_rate = args[2]
    domestic_curve = args[3]
    foreign_curve = args[4]

    model = BlackScholes(volatility)

    fprime = self.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    return fprime


###############################################################################


@njit(fastmath=True, cache=True)
def fast_delta(s, t, k, rd, rf, vol, deltaTypeValue, option_type_value):
    """Calculation of the FX Option delta. Used in the determination of
    the volatility surface. Avoids discount curve interpolation so it
    should be slightly faster than the full calculation of delta."""

    pips_spot_delta = bs_delta(s, t, k, rd, rf, vol, option_type_value)

    if deltaTypeValue == FinFXDeltaMethod.SPOT_DELTA.value:
        return pips_spot_delta
    elif deltaTypeValue == FinFXDeltaMethod.FORWARD_DELTA.value:
        pips_fwd_delta = pips_spot_delta * np.exp(rf * t)
        return pips_fwd_delta
    elif deltaTypeValue == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:
        vpctf = bs_value(s, t, k, rd, rf, vol, option_type_value) / s
        pct_spot_delta_prem_adj = pips_spot_delta - vpctf
        return pct_spot_delta_prem_adj
    elif deltaTypeValue == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:
        vpctf = bs_value(s, t, k, rd, rf, vol, option_type_value) / s
        pct_fwd_delta_prem_adj = np.exp(rf * t) * (pips_spot_delta - vpctf)
        return pct_fwd_delta_prem_adj
    else:
        raise FinError("Unknown FinFXDeltaMethod")


###############################################################################

# def g(K, *args):
#     """ This is the objective function used in the determination of the FX
#     Option implied strike which is computed in the class below. """

#     self = args[0]
#     value_dt = args[1]
#     stock_price = args[2]
#     domestic_curve = args[3]
#     foreign_curve = args[4]
#     delta = args[5]
#     deltaType = args[6]
#     volatility = args[7]

#     model = FinFXModelBlackScholes(volatility)

#     self.strike_fx_rate = K

#     deltaDict = self.delta(value_dt,
#                            stock_price,
#                            domestic_curve,
#                            foreign_curve,
#                            model)

#     delta_out = deltaDict[deltaType]
#     obj_fn = delta - delta_out
#     return obj_fn

# ## THIS IS A HOPEFULLY FASTER VERSION WHICH AVOIDS CALLING DF

# def g2(K, *args):
#     """ This is the objective function used in the determination of the FX
#     Option implied strike which is computed in the class below. """

#     self = args[0]
#     value_dt = args[1]
#     stock_price = args[2]
#     dom_df = args[3]
#     for_df = args[4]
#     delta = args[5]
#     deltaType = args[6]
#     volatility = args[7]

#     self.strike_fx_rate = K

#     deltaDict = self.fast_delta(value_dt,
#                                stock_price,
#                                dom_df,
#                                for_df,
#                                volatility)

#     delta_out = deltaDict[deltaType]
#     obj_fn = delta - delta_out
#     return obj_fn


###############################################################################

###############################################################################
# ALL CCY RATES MUST BE IN NUM UNITS OF DOMESTIC PER UNIT OF FOREIGN CURRENCY
# SO EURUSD = 1.30 MEANS 1.30 DOLLARS PER EURO SO DOLLAR IS THE DOMESTIC AND
# EUR IS THE FOREIGN CURRENCY
###############################################################################


class FXVanillaOption:
    """This is a class for an FX Option trade. It permits the user to
    calculate the price of an FX Option trade which can be expressed in a
    number of ways depending on the investor or hedger's currency. It aslo
    allows the calculation of the option's delta in a number of forms as
    well as the various Greek risk sensitivies."""

    def __init__(
        self,
        expiry_dt: Date,
        # 1 unit of foreign in domestic
        strike_fx_rate: (float, np.ndarray),
        currency_pair: str,  # FORDOM
        option_type: (OptionTypes, list),
        notional: float,
        prem_currency: str,
        spot_days: int = 0,
    ):
        """Create the FX Vanilla Option object. Inputs include expiry date,
        strike, currency pair, option type (call or put), notional and the
        currency of the notional. And adjustment for spot days is enabled. All
        currency rates must be entered in the price in domestic currency of
        one unit of foreign. And the currency pair should be in the form FORDOM
        where FOR is the foreign currency pair currency code and DOM is the
        same for the domestic currency."""

        check_argument_types(self.__init__, locals())

        delivery_dt = expiry_dt.add_weekdays(spot_days)

        """ The FX rate the price in domestic currency ccy2 of a single unit
        of the foreign currency which is ccy1. For example EURUSD of 1.3 is the
        price in USD (CCY2) of 1 unit of EUR (CCY1)"""

        if delivery_dt < expiry_dt:
            raise FinError("Delivery date must be on or after expiry date.")

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.expiry_dt = expiry_dt
        self.delivery_dt = delivery_dt

        if np.any(strike_fx_rate < 0.0):
            raise FinError("Negative strike.")

        self.strike_fx_rate = strike_fx_rate

        self.currency_pair = currency_pair
        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        if prem_currency != self.dom_name and prem_currency != self.for_name:
            raise FinError("Premium currency not in currency pair.")

        self.prem_currency = prem_currency

        self.notional = notional

        if (
            option_type != OptionTypes.EUROPEAN_CALL
            and option_type != OptionTypes.EUROPEAN_PUT
            and option_type != OptionTypes.AMERICAN_CALL
            and option_type != OptionTypes.AMERICAN_PUT
        ):
            raise FinError("Unknown Option Type:" + option_type)

        self.option_type = option_type
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
        """This function calculates the value of the option using a specified
        model with the resulting value being in domestic i.e. ccy2 terms.
        Recall that Domestic = CCY2 and Foreign = CCY1 and FX rate is in
        price in domestic of one unit of foreign currency."""

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

        if t_del < 0.0:
            raise FinError("Time to expiry must be positive.")

        t_del = np.maximum(t_del, 1e-10)

        # TODO RESOLVE t_del versus TEXP
        dom_df = domestic_curve.df_t(t_del)
        for_df = foreign_curve.df_t(t_del)

        r_d = -np.log(dom_df) / t_del
        r_f = -np.log(for_df) / t_del

        s0 = spot_fx_rate
        k = self.strike_fx_rate
        f0t = s0 * np.exp((r_d - r_f) * t_del)

        if isinstance(model, BlackScholes) or isinstance(model, SABR):

            if isinstance(model, BlackScholes):
                volatility = model.volatility
            elif isinstance(model, SABR):

                params_list = np.array(
                    [model.alpha, model.beta, model.rho, model.nu]
                )

                volatility = vol_function_sabr(params_list, f0t, k, t_del)

            if np.any(volatility < 0.0):
                raise FinError("Volatility should not be negative.")

            v = np.maximum(volatility, 1e-10)

            if self.option_type == OptionTypes.EUROPEAN_CALL:

                vdf = bs_value(
                    s0, t_exp, k, r_d, r_f, v, OptionTypes.EUROPEAN_CALL.value
                )

            elif self.option_type == OptionTypes.EUROPEAN_PUT:

                vdf = bs_value(
                    s0, t_exp, k, r_d, r_f, v, OptionTypes.EUROPEAN_PUT.value
                )

            elif self.option_type == OptionTypes.AMERICAN_CALL:

                num_steps_per_year = 100

                vdf = crr_tree_val_avg(
                    s0,
                    r_d,
                    r_f,
                    volatility,
                    num_steps_per_year,
                    t_exp,
                    OptionTypes.AMERICAN_CALL.value,
                    k,
                )["value"]

            elif self.option_type == OptionTypes.AMERICAN_PUT:

                num_steps_per_year = 100

                vdf = crr_tree_val_avg(
                    s0,
                    r_d,
                    r_f,
                    volatility,
                    num_steps_per_year,
                    t_exp,
                    OptionTypes.AMERICAN_PUT.value,
                    k,
                )["value"]
            else:
                raise FinError("Unknown option type")

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

        vdf = vdf
        pips_dom = vdf
        pips_for = vdf / (spot_fx_rate * self.strike_fx_rate)

        cash_dom = vdf * notional_dom / self.strike_fx_rate
        cash_for = vdf * notional_for / spot_fx_rate

        pct_dom = vdf / self.strike_fx_rate
        pct_for = vdf / spot_fx_rate

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
        value_dt,
        spot_fx_rate,
        ccy1DiscountCurve,
        ccy2DiscountCurve,
        model,
    ):
        """Calculation of the FX option delta by bumping the spot FX rate by
        1 cent of its value. This gives the FX spot delta. For speed we prefer
        to use the analytical calculation of the derivative given below."""

        bump = 0.0001 * spot_fx_rate

        v = self.value(
            value_dt, spot_fx_rate, ccy1DiscountCurve, ccy2DiscountCurve, model
        )

        v_bumped = self.value(
            value_dt,
            spot_fx_rate + bump,
            ccy1DiscountCurve,
            ccy2DiscountCurve,
            model,
        )

        if type(v_bumped) is dict:
            delta = (v_bumped["value"] - v["value"]) / bump
        else:
            delta = (v_bumped - v) / bump

        return delta

    ###########################################################################

    def delta(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculation of the FX Option delta. There are several definitions
        of delta and so we are required to return a dictionary of values. The
        definitions can be found on Page 44 of Foreign Exchange Option Pricing
        by Iain Clark, published by Wiley Finance."""

        if isinstance(value_dt, Date):
            spot_dt = value_dt.add_weekdays(self.spot_days)
            t_del = (self.delivery_dt - spot_dt) / g_days_in_year
            t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t_del = value_dt
            t_exp = t_del

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if np.any(t_del < 0.0):
            raise FinError("Time to expiry must be positive.")

        t_del = np.maximum(t_del, 1e-10)

        dom_df = domestic_curve.df_t(t_del)
        r_d = -np.log(dom_df) / t_del

        for_df = foreign_curve.df_t(t_del)
        r_f = -np.log(for_df) / t_del

        s0 = spot_fx_rate
        K = self.strike_fx_rate

        if isinstance(model, BlackScholes):

            v = model.volatility

            if np.any(v < 0.0):
                raise FinError("Volatility should not be negative.")

            v = np.maximum(v, g_small)

            pips_spot_delta = bs_delta(
                s0, t_exp, K, r_d, r_f, v, self.option_type.value
            )
            pips_fwd_delta = pips_spot_delta * np.exp(r_f * t_del)
            vpctf = (
                bs_value(s0, t_exp, K, r_d, r_f, v, self.option_type.value)
                / s0
            )
            pct_spot_delta_prem_adj = pips_spot_delta - vpctf
            pct_fwd_delta_prem_adj = np.exp(r_f * t_del) * (
                pips_spot_delta - vpctf
            )

        return {
            "pips_spot_delta": pips_spot_delta,
            "pips_fwd_delta": pips_fwd_delta,
            "pct_spot_delta_prem_adj": pct_spot_delta_prem_adj,
            "pct_fwd_delta_prem_adj": pct_fwd_delta_prem_adj,
        }

    ###########################################################################

    def fast_delta(self, t, s, rd, rf, vol):
        """Calculation of the FX Option delta. Used in the determination of
        the volatility surface. Avoids discount curve interpolation so it
        should be slightly faster than the full calculation of delta."""

        #        spot_dt = value_dt.add_weekdays(self.spot_days)
        #        t_del = (self.delivery_dt - value_dt) / g_days_in_year
        #        t_del = np.maximum(t_del, g_small)

        #        r_d = -np.log(dom_df)/t_del
        #         r_f = -np.log(for_df)/t_del
        k = self.strike_fx_rate

        #        print("FAST DELTA IN OPTION CLASS", s,t,k,rd,rf,vol)

        pips_spot_delta = bs_delta(
            s, t, k, rd, rf, vol, self.option_type.value
        )
        pips_fwd_delta = pips_spot_delta * np.exp(rf * t)

        vpctf = bs_value(s, t, k, rd, rf, vol, self.option_type.value) / s

        pct_spot_delta_prem_adj = pips_spot_delta - vpctf
        pct_fwd_delta_prem_adj = np.exp(rf * t) * (pips_spot_delta - vpctf)

        return {
            "pips_spot_delta": pips_spot_delta,
            "pips_fwd_delta": pips_fwd_delta,
            "pct_spot_delta_prem_adj": pct_spot_delta_prem_adj,
            "pct_fwd_delta_prem_adj": pct_fwd_delta_prem_adj,
        }

    ###########################################################################

    def gamma(
        self,
        value_dt,
        spot_fx_rate,  # value of a unit of foreign in domestic currency
        domestic_curve,
        foreign_curve,
        model,
    ):
        """This function calculates the FX Option Gamma using spot delta."""

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("FX Rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        dom_df = domestic_curve.df_t(t)
        r_d = -np.log(dom_df) / t

        for_df = foreign_curve.df_t(t)
        r_f = -np.log(for_df) / t

        K = self.strike_fx_rate
        s0 = spot_fx_rate

        if isinstance(model, BlackScholes):

            volatility = model.volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            ln_s0_k = np.log(s0 / K)
            sqrt_t = np.sqrt(t)
            den = volatility * sqrt_t
            mu = r_d - r_f
            v2 = volatility * volatility
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / den
            gamma = np.exp(-r_f * t) * nprime(d1)
            gamma = gamma / s0 / den
        else:
            raise FinError("Unknown Model Type")

        return gamma

    ###########################################################################

    def vega(
        self,
        value_dt,
        spot_fx_rate,  # value of a unit of foreign in domestic currency
        domestic_curve,
        foreign_curve,
        model,
    ):
        """This function calculates the FX Option Vega using the spot delta."""

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        dom_df = domestic_curve.df_t(t)
        r_d = -np.log(dom_df) / t

        for_df = foreign_curve.df_t(t)
        r_f = -np.log(for_df) / t

        K = self.strike_fx_rate
        s0 = spot_fx_rate

        if isinstance(model, BlackScholes):

            volatility = model.volatility

            if np.any(volatility) < 0.0:
                raise FinError("Volatility should not be negative.")

            volatility = np.maximum(volatility, 1e-10)

            ln_s0_k = np.log(s0 / K)
            sqrt_t = np.sqrt(t)
            den = volatility * sqrt_t
            mu = r_d - r_f
            v2 = volatility * volatility
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / den
            vega = s0 * sqrt_t * np.exp(-r_f * t) * nprime(d1)
        else:
            raise FinError("Unknown Model type")

        return vega

    ###########################################################################

    def theta(
        self,
        value_dt,
        spot_fx_rate,  # value of a unit of foreign in domestic currency
        domestic_curve,
        foreign_curve,
        model,
    ):
        """This function calculates the time decay of the FX option."""

        if isinstance(value_dt, Date):
            t = (self.expiry_dt - value_dt) / g_days_in_year
        else:
            t = value_dt

        if np.any(spot_fx_rate <= 0.0):
            raise FinError("Spot FX Rate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        dom_df = domestic_curve.df_t(t)
        r_d = -np.log(dom_df) / t

        for_df = foreign_curve.df_t(t)
        r_f = -np.log(for_df) / t

        K = self.strike_fx_rate
        s0 = spot_fx_rate

        if isinstance(model, BlackScholes):

            vol = model.volatility

            if np.any(vol) < 0.0:
                raise FinError("Volatility should not be negative.")

            vol = np.maximum(vol, 1e-10)

            ln_s0_k = np.log(s0 / K)
            sqrt_t = np.sqrt(t)
            den = vol * sqrt_t
            mu = r_d - r_f
            v2 = vol * vol
            d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / den
            d2 = (ln_s0_k + (mu - v2 / 2.0) * t) / den

            if self.option_type == OptionTypes.EUROPEAN_CALL:
                v = -s0 * np.exp(-r_f * t) * nprime(d1) * vol / 2.0 / sqrt_t
                v = v + r_f * s0 * np.exp(-r_f * t) * N(d1)
                v = v - r_d * K * np.exp(-r_d * t) * N(d2)
            elif self.option_type == OptionTypes.EUROPEAN_PUT:
                v = -s0 * np.exp(-r_f * t) * nprime(d1) * vol / 2.0 / sqrt_t
                v = v + r_d * K * np.exp(-r_d * t) * N(-d2)
                v = v - r_f * s0 * np.exp(-r_f * t) * N(-d1)
            else:
                raise FinError("Unknown option type")

        else:
            raise FinError("Unknown Model Type")

        return v

    ###########################################################################

    def implied_volatility(
        self, value_dt, stock_price, discount_curve, dividend_curve, price
    ):
        """This function determines the implied volatility of an FX option
        given a price and the other option details. It uses a one-dimensional
        Newton root search algorith to determine the implied volatility."""

        argtuple = (
            self,
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            price,
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

    def value_mc(
        self,
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
        num_paths=10000,
        seed=4242,
    ):
        """Calculate the value of an FX Option using Monte Carlo methods.
        This function can be used to validate the risk measures calculated
        above or used as the starting code for a model exotic FX product that
        cannot be priced analytically. This function uses Numpy vectorisation
        for speed of execution."""

        if isinstance(model, BlackScholes):
            volatility = model.volatility
        else:
            raise FinError("Model Type invalid")

        np.random.seed(seed)
        t = (self.expiry_dt - value_dt) / g_days_in_year

        dom_Df = domestic_curve.df(self.expiry_dt)
        for_Df = foreign_curve.df(self.expiry_dt)

        r_d = -np.log(dom_Df) / t
        r_f = -np.log(for_Df) / t

        mu = r_d - r_f
        v2 = volatility**2
        K = self.strike_fx_rate
        sqrt_dt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, num_paths))
        s = spot_fx_rate * np.exp((mu - v2 / 2.0) * t)
        m = np.exp(g * sqrt_dt * volatility)
        s_1 = s * m
        s_2 = s / m

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            payoff_a_1 = np.maximum(s_1 - K, 0.0)
            payoff_a_2 = np.maximum(s_2 - K, 0.0)
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            payoff_a_1 = np.maximum(K - s_1, 0.0)
            payoff_a_2 = np.maximum(K - s_2, 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * np.exp(-r_d * t) / 2.0
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("CURRENCY PAIR", self.currency_pair)
        s += label_to_string("PREMIUM CCY", self.prem_currency)
        s += label_to_string("STRIKE FX RATE", self.strike_fx_rate)
        s += label_to_string("OPTION TYPE", self.option_type)
        s += label_to_string("SPOT DAYS", self.spot_days)
        s += label_to_string("NOTIONAL", self.notional, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
