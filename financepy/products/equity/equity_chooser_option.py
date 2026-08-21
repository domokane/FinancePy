##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from scipy import optimize

from ...utils.math import M
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.global_vars import G_SMALL
from ...utils.error import FinError

from ...products.equity.equity_option import EquityOption
from ...utils.global_types import OptionTypes
from ...market.curves.flat_discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...models.black_scholes_analytic import european_value

DEBUG_MODE = False

########################################################################################
# TODO: Vectorise pricer
# TODO: NUMBA ??
# TODO: Monte Carlo pricer
########################################################################################


def _f(ss, *args):
    """Complex chooser option solve for critical stock price that makes the
    forward starting call and put options have the same price on the chooser
    date."""

    t = args[0]
    tc = args[1]
    tp = args[2]
    rtc = args[3]
    rtp = args[4]
    kc = args[5]
    kp = args[6]
    v = args[7]
    qc = args[8]
    qp = args[9]

    call_int = OptionTypes.EUROPEAN_CALL.value
    put_int = OptionTypes.EUROPEAN_PUT.value

    if tc == t:
        v_call = np.maximum(ss - kc, 0.0)
    else:

        v_call = european_value(ss, tc - t, kc, rtc, qc, v, call_int)

    if tp == t:
        v_put = np.maximum(kp - ss, 0.0)
    else:
        v_put = european_value(ss, tp - t, kp, rtp, qp, v, put_int)

    v = v_call - v_put
    return v


########################################################################################


class EquityChooserOption(EquityOption):
    """A EquityChooserOption is an option which allows the holder to
    either enter into a call or a put option on a later expiry date, with both
    strikes potentially different and both expiry dates potentially different.
    This is known as a complex chooser. All the option details are set at trade
    initiation."""

    def __init__(
        self,
        choose_dt: Date,
        call_expiry_dt: Date,
        put_expiry_dt: Date,
        call_strike_price: float,
        put_strike_price: float,
    ):
        """Create the EquityChooserOption by passing in the chooser date
        and then the put and call expiry dates as well as the corresponding put
        and call strike prices."""

        check_argument_types(self.__init__, locals())

        if choose_dt > call_expiry_dt:
            raise FinError("Expiry date must precede call option expiry date")

        if choose_dt > put_expiry_dt:
            raise FinError("Expiry date must precede put option expiry date")

        self.choose_dt = choose_dt
        self.call_expiry_dt = call_expiry_dt
        self.put_expiry_dt = put_expiry_dt
        self.call_strike = float(call_strike_price)
        self.put_strike = float(put_strike_price)

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Value the complex chooser option using an approach by Rubinstein
        (1991). See also Haug page 129 for complex chooser options."""

        if value_dt > self.choose_dt:
            raise FinError("Value date after choose date.")

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.call_expiry_dt:
            raise FinError("Valuation date after call expiry date.")

        if value_dt > self.put_expiry_dt:
            raise FinError("Valuation date after put expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError("Discount Curve valuation date not same as option value date")

        if dividend_curve.value_dt != value_dt:
            raise FinError("Dividend Curve valuation date not same as option value date")

        if value_dt == self.choose_dt:
            v = self.value_dt_on_choose_dt(value_dt, stock_price, discount_curve, dividend_curve, model)
            return v

        t_choose = (self.choose_dt - value_dt) / G_DAYS_IN_YEAR
        t_call = (self.call_expiry_dt - value_dt) / G_DAYS_IN_YEAR
        t_put = (self.put_expiry_dt - value_dt) / G_DAYS_IN_YEAR

        df_t = discount_curve.df(self.choose_dt)
        df_c = discount_curve.df(self.call_expiry_dt)
        df_p = discount_curve.df(self.put_expiry_dt)

        dq_t = dividend_curve.df(self.choose_dt)
        dq_c = dividend_curve.df(self.call_expiry_dt)
        dq_p = dividend_curve.df(self.put_expiry_dt)

        if t_call > t_choose:
            rfc = -np.log(df_c / df_t) / (t_call - t_choose)
            qfc = -np.log(dq_c / dq_t) / (t_call - t_choose)
        else:
            rfc = 0.0
            qfc = 0.0

        if t_put > t_choose:
            rfp = -np.log(df_p / df_t) / (t_put - t_choose)
            qfp = -np.log(dq_p / dq_t) / (t_put - t_choose)
        else:
            rfp = 0.0
            qfp = 0.0

        vol = model.volatility
        vol = max(vol, G_SMALL)
        vol2 = vol * vol

        scalar_input = np.isscalar(stock_price)
        s0 = np.atleast_1d(np.asarray(stock_price, dtype=float))

        xc = self.call_strike
        xp = self.put_strike

        argtuple = (t_choose, t_call, t_put, rfc, rfp, xc, xp, vol, qfc, qfp)
        if DEBUG_MODE:
            print("args", argtuple)

        x_init = 0.5 * (xc + xp)
        istar = optimize.newton(_f, x0=x_init, args=argtuple, tol=1e-8, maxiter=50)

        if DEBUG_MODE:
            print("istar", istar)

        sqrt_tc = np.sqrt(t_call)
        sqrt_tp = np.sqrt(t_put)
        sqrt_t = np.sqrt(t_choose)

        d1 = (np.log(s0 / istar) + np.log(dq_t / df_t) + 0.5 * vol2 * t_choose) / vol / sqrt_t
        d2 = d1 - vol * sqrt_t

        if DEBUG_MODE:
            print("d1", d1)
            print("d2", d2)

        y1 = (np.log(s0 / xc) + np.log(dq_c / df_c) + 0.5 * vol2 * t_call) / vol / sqrt_tc
        y2 = (np.log(s0 / xp) + np.log(dq_p / df_p) + 0.5 * vol2 * t_put) / vol / sqrt_tp

        if DEBUG_MODE:
            print("y1", y1)
            print("y2", y2)

        rho1 = sqrt_t / sqrt_tc
        rho2 = sqrt_t / sqrt_tp

        if DEBUG_MODE:
            print("rho1", rho1)
            print("rho2", rho2)

        if 1 == 1:
            w = s0 * dq_c * M(d1, y1, rho1)
            w = w - xc * df_c * M(d2, y1 - vol * sqrt_tc, rho1)
            w = w - s0 * dq_p * M(-d1, -y2, rho2)
            w = w + xp * df_p * M(-d2, -y2 + vol * sqrt_tp, rho2)
        else:
            m1 = np.array([M(a, b, rho1) for a, b in zip(d1, y1)])
            m2 = np.array([M(a, b, rho1) for a, b in zip(d2, y1 - vol * sqrt_tc)])
            m3 = np.array([M(a, b, rho2) for a, b in zip(-d1, -y2)])
            m4 = np.array([M(a, b, rho2) for a, b in zip(-d2, -y2 + vol * sqrt_tp)])

            w = s0 * dq_c * m1
            w -= xc * df_c * m2
            w -= s0 * dq_p * m3
            w += xp * df_p * m4

        if scalar_input:
            return w[0]

        return w

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_paths: int = 10000,
        seed: int = 4242,
    ):
        """Value the complex chooser option Monte Carlo."""

        if value_dt == self.choose_dt:
            v = self.value_dt_on_choose_dt(value_dt, stock_price, discount_curve, dividend_curve, model)
            return v

        t = (self.choose_dt - value_dt) / G_DAYS_IN_YEAR
        t_c = (self.call_expiry_dt - value_dt) / G_DAYS_IN_YEAR
        t_p = (self.put_expiry_dt - value_dt) / G_DAYS_IN_YEAR

        vol = model.volatility
        vol = max(vol, 1e-6)

        df_t = discount_curve.df(self.choose_dt)
        df_c = discount_curve.df(self.call_expiry_dt)
        df_p = discount_curve.df(self.put_expiry_dt)

        dq_t = dividend_curve.df(self.choose_dt)
        dq_c = dividend_curve.df(self.call_expiry_dt)
        dq_p = dividend_curve.df(self.put_expiry_dt)

        if t_c > t:
            rfc = -np.log(df_c / df_t) / (t_c - t)
            qfc = -np.log(dq_c / dq_t) / (t_c - t)
        else:
            rfc = 0.0
            qfc = 0.0

        if t_p > t:
            rfp = -np.log(df_p / df_t) / (t_p - t)
            qfp = -np.log(dq_p / dq_t) / (t_p - t)
        else:
            rfp = 0.0
            qfp = 0.0

        kc = self.call_strike
        kp = self.put_strike

        rng = np.random.default_rng(seed)
        g = rng.normal(0.0, 1.0, size=(1, num_paths))
        sqrt_dt = np.sqrt(t)

        forward_growth = dq_t / df_t
        s = stock_price * forward_growth * np.exp(-0.5 * vol * vol * t)
        m = np.exp(g * sqrt_dt * vol)

        s_1 = s * m
        s_2 = s / m

        if t_c == t:
            v_call_1 = np.maximum(s_1 - kc, 0.0)
            v_call_2 = np.maximum(s_2 - kc, 0.0)
        else:
            v_call_1 = european_value(s_1, t_c - t, kc, rfc, qfc, vol, OptionTypes.EUROPEAN_CALL.value)
            v_call_2 = european_value(s_2, t_c - t, kc, rfc, qfc, vol, OptionTypes.EUROPEAN_CALL.value)

        if t_p == t:
            v_put_1 = np.maximum(kp - s_1, 0.0)
            v_put_2 = np.maximum(kp - s_2, 0.0)
        else:
            v_put_1 = european_value(s_1, t_p - t, kp, rfp, qfp, vol, OptionTypes.EUROPEAN_PUT.value)
            v_put_2 = european_value(s_2, t_p - t, kp, rfp, qfp, vol, OptionTypes.EUROPEAN_PUT.value)

        payoff_1 = np.maximum(v_call_1, v_put_1)
        payoff_2 = np.maximum(v_call_2, v_put_2)

        payoff = np.mean(payoff_1) + np.mean(payoff_2)
        v = payoff * df_t / 2.0
        return v

    ###########################################################################

    def value_dt_on_choose_dt(self, value_dt, stock_price, discount_curve, dividend_curve, model):

        stock_price = np.asarray(stock_price, dtype=float)

        vol = max(model.volatility, G_SMALL)

        t_call = (self.call_expiry_dt - value_dt) / G_DAYS_IN_YEAR
        t_put = (self.put_expiry_dt - value_dt) / G_DAYS_IN_YEAR

        if t_call == 0.0:
            call_value = np.maximum(stock_price - self.call_strike, 0.0)
        else:
            df_c = discount_curve.df(self.call_expiry_dt)
            dq_c = dividend_curve.df(self.call_expiry_dt)

            rc = -np.log(df_c) / t_call
            qc = -np.log(dq_c) / t_call

            call_value = european_value(
                stock_price,
                t_call,
                self.call_strike,
                rc,
                qc,
                vol,
                OptionTypes.EUROPEAN_CALL.value,
            )

        if t_put == 0.0:
            put_value = np.maximum(self.put_strike - stock_price, 0.0)
        else:
            df_p = discount_curve.df(self.put_expiry_dt)
            dq_p = dividend_curve.df(self.put_expiry_dt)

            rp = -np.log(df_p) / t_put
            qp = -np.log(dq_p) / t_put

            put_value = european_value(
                stock_price,
                t_put,
                self.put_strike,
                rp,
                qp,
                vol,
                OptionTypes.EUROPEAN_PUT.value,
            )

        return np.maximum(call_value, put_value)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("CHOOSER DATE", self.choose_dt)
        s += label_to_string("CALL EXPIRY DATE", self.call_expiry_dt)
        s += label_to_string("CALL STRIKE PRICE", self.call_strike)
        s += label_to_string("PUT EXPIRY DATE", self.put_expiry_dt)
        s += label_to_string("PUT STRIKE PRICE", self.put_strike, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
