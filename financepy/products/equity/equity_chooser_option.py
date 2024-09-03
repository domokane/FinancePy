##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from scipy import optimize

from ...utils.math import M
from ...utils.global_vars import g_days_in_year
from ...utils.global_vars import g_small
from ...utils.error import FinError

from ...products.equity.equity_option import EquityOption
from ...utils.global_types import OptionTypes
from ...market.curves.discount_curve_flat import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...models.black_scholes_analytic import bs_value

from scipy.stats import norm

N = norm.cdf


###############################################################################
# TODO: Vectorise pricer
# TODO: NUMBA ??
# TODO: Monte Carlo pricer
###############################################################################


def _f(ss, *args):
    """ Complex chooser option solve for critical stock price that makes the
    forward starting call and put options have the same price on the chooser
    date. """

    t = args[0]
    tc = args[1]
    tp = args[2]
    rtc = args[3]
    rtp = args[4]
    kc = args[5]
    kp = args[6]
    v = args[7]
    q = args[8]

    v_call = bs_value(ss, tc - t, kc, rtc, q, v,
                      OptionTypes.EUROPEAN_CALL.value)
    v_put = bs_value(ss, tp - t, kp, rtp, q, v,
                     OptionTypes.EUROPEAN_PUT.value)

    v = v_call - v_put
    return v

###############################################################################


class EquityChooserOption(EquityOption):
    """ A EquityChooserOption is an option which allows the holder to
    either enter into a call or a put option on a later expiry date, with both
    strikes potentially different and both expiry dates potentially different.
    This is known as a complex chooser. All the option details are set at trade
    initiation. """

    def __init__(self,
                 choose_dt: Date,
                 call_expiry_dt: Date,
                 put_expiry_dt: Date,
                 call_strike_price: float,
                 put_strike_price: float):
        """ Create the EquityChooserOption by passing in the chooser date
        and then the put and call expiry dates as well as the corresponding put
        and call strike prices. """

        check_argument_types(self.__init__, locals())

        if choose_dt > call_expiry_dt:
            raise FinError("Expiry date must precede call option expiry date")

        if choose_dt > put_expiry_dt:
            raise FinError("Expiry date must precede put option expiry date")

        self.chooseDate = choose_dt
        self.call_expiry_dt = call_expiry_dt
        self.put_expiry_dt = put_expiry_dt
        self.call_strike = float(call_strike_price)
        self.put_strike = float(put_strike_price)

    ###########################################################################

    def value(self,
              value_dt: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Value the complex chooser option using an approach by Rubinstein
        (1991). See also Haug page 129 for complex chooser options. """

        if value_dt > self.chooseDate:
            raise FinError("Value date after choose date.")

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.call_expiry_dt:
            raise FinError("Valuation date after call expiry date.")

        if value_dt > self.put_expiry_dt:
            raise FinError("Valuation date after put expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date")

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date")

        DEBUG_MODE = False

        t = (self.chooseDate - value_dt) / g_days_in_year
        tc = (self.call_expiry_dt - value_dt) / g_days_in_year
        tp = (self.put_expiry_dt - value_dt) / g_days_in_year

        rt = discount_curve.cc_rate(self.chooseDate)
        rtc = discount_curve.cc_rate(self.call_expiry_dt)
        rtp = discount_curve.cc_rate(self.put_expiry_dt)

        q = dividend_curve.cc_rate(self.chooseDate)

        t = max(t, g_small)
        tc = max(tc, g_small)
        tp = max(tp, g_small)

        v = model.volatility
        v = max(v, g_small)

        s0 = stock_price
        xc = self.call_strike
        xp = self.put_strike
        bt = rt - q
        btc = rtc - q
        btp = rtp - q

        argtuple = (t, tc, tp, rtc, rtp, xc, xp, v, q)
        if DEBUG_MODE:
            print("args", argtuple)

        istar = optimize.newton(_f, x0=s0, args=argtuple, tol=1e-8,
                                maxiter=50, fprime2=None)

        if DEBUG_MODE:
            print("istar", istar)

        d1 = (np.log(s0 / istar) + (bt + v * v / 2) * t) / v / np.sqrt(t)
        d2 = d1 - v * np.sqrt(t)

        if DEBUG_MODE:
            print("d1", d1)
            print("d2", d2)

        y1 = (np.log(s0 / xc) + (btc + v * v / 2) * tc) / v / np.sqrt(tc)
        y2 = (np.log(s0 / xp) + (btp + v * v / 2) * tp) / v / np.sqrt(tp)

        if DEBUG_MODE:
            print("y1", y1)
            print("y2", y2)

        rho1 = np.sqrt(t / tc)
        rho2 = np.sqrt(t / tp)

        if DEBUG_MODE:
            print("rho1", rho1)
            print("rho2", rho2)

        w = s0 * np.exp(-q * tc) * M(d1, y1, rho1)
        w = w - xc * np.exp(-rtc * tc) * M(d2, y1 - v * np.sqrt(tc), rho1)
        w = w - s0 * np.exp(-q * tp) * M(-d1, -y2, rho2)
        w = w + xp * np.exp(-rtp * tp) * M(-d2, -y2 + v * np.sqrt(tp), rho2)
        return w

    ###########################################################################

    def value_mc(self,
                 value_dt: Date,
                 stock_price: float,
                 discount_curve: DiscountCurve,
                 dividend_curve: DiscountCurve,
                 model,
                 num_paths: int = 10000,
                 seed: int = 4242):
        """ Value the complex chooser option Monte Carlo. """

        dft = discount_curve.df(self.chooseDate)
        dftc = discount_curve.df(self.call_expiry_dt)
        dftp = discount_curve.df(self.put_expiry_dt)

        t = (self.chooseDate - value_dt) / g_days_in_year
        tc = (self.call_expiry_dt - value_dt) / g_days_in_year
        tp = (self.put_expiry_dt - value_dt) / g_days_in_year

        rt = -np.log(dft) / t
        rtc = -np.log(dftc) / tc
        rtp = -np.log(dftp) / tp

        t = max(t, 1e-6)
        tc = max(tc, 1e-6)
        tp = max(tp, 1e-6)

        v = model.volatility
        v = max(v, 1e-6)

        # SHOULD THIS CARE ABOUT TERM STRUCTURE OF Q
        dq = dividend_curve.df(self.chooseDate)
        q = -np.log(dq) / t

        #        q = dividend_yield
        kc = self.call_strike
        kp = self.put_strike

        np.random.seed(seed)
        sqrt_dt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, num_paths))
        s = stock_price * np.exp((rt - q - v * v / 2.0) * t)
        m = np.exp(g * sqrt_dt * v)

        s_1 = s * m
        s_2 = s / m

        v_call_1 = bs_value(s_1, tc - t, kc, rtc, q, v,
                            OptionTypes.EUROPEAN_CALL.value)
        v_put_1 = bs_value(s_1, tp - t, kp, rtp, q, v,
                           OptionTypes.EUROPEAN_PUT.value)

        v_call_2 = bs_value(s_2, tc - t, kc, rtc, q, v,
                            OptionTypes.EUROPEAN_CALL.value)
        v_put_2 = bs_value(s_2, tp - t, kp, rtp, q, v,
                           OptionTypes.EUROPEAN_PUT.value)

        payoff_1 = np.maximum(v_call_1, v_put_1)
        payoff_2 = np.maximum(v_call_2, v_put_2)

        payoff = np.mean(payoff_1) + np.mean(payoff_2)
        v = payoff * dft / 2.0
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("CHOOSER DATE", self.chooseDate)
        s += label_to_string("CALL EXPIRY DATE", self.call_expiry_dt)
        s += label_to_string("CALL STRIKE PRICE", self.call_strike)
        s += label_to_string("PUT EXPIRY DATE", self.put_expiry_dt)
        s += label_to_string("PUT STRIKE PRICE", self.put_strike, "")
        return s

    ###########################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
