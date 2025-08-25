##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import ONE_MILLION
from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.global_types import OptionTypes
from ...models.black_scholes import BlackScholes
from ...utils.helpers import label_to_string, check_argument_types

from .equity_vanilla_option import EquityVanillaOption

########################################################################################


class EquityVarianceSwap:
    """Class for managing an equity variance swap contract."""

    def __init__(
        self,
        start_dt: Date,
        maturity_dt_or_tenor: Union[Date, str],
        strike_variance: float,
        notional: float = ONE_MILLION,
        pay_strike_flag: bool = True,
    ):
        """Create variance swap contract."""

        check_argument_types(self.__init__, locals())

        if isinstance(maturity_dt_or_tenor, Date):
            maturity_dt = maturity_dt_or_tenor
        else:
            maturity_dt = start_dt.add_tenor(maturity_dt_or_tenor)

        if start_dt >= maturity_dt:
            raise FinError("Start date after or same as maturity date")

        self.start_dt = start_dt
        self.maturity_dt = maturity_dt
        self.strike_variance = strike_variance
        self.notional = notional
        self.pay_strike_flag = pay_strike_flag

        # Replication portfolio is stored
        self.num_put_options = 0
        self.num_call_options = 0
        self.put_wts = []
        self.put_strikes = []
        self.call_wts = []
        self.call_strikes = []

    ###########################################################################

    def value(self, value_dt, realised_var, fair_strike_var, libor_curve):
        """Calculate the value of the variance swap based on the realised
        volatility to the valuation date, the forward looking implied
        volatility to the maturity date using the libor discount curve."""

        t1 = (value_dt - self.start_dt) / G_DAYS_IN_YEARS
        t2 = (self.maturity_dt - self.start_dt) / G_DAYS_IN_YEARS

        expected_variance = t1 * realised_var / t2
        expected_variance += (t2 - t1) * fair_strike_var / t2

        payoff = expected_variance - self.strike_variance

        df = libor_curve.df(self.maturity_dt)
        v = payoff * self.notional * df
        return v

    ###########################################################################

    def fair_strike_approx(self, value_dt, fwd_stock_price, strikes, volatilities):
        """This is an approximation of the fair strike variance by Demeterfi
        et al. (1999) which assumes that sigma(K) = sigma(F) - b(K-F)/F where
        F is the forward stock price and sigma(F) is the ATM forward vol."""

        f = fwd_stock_price

        # TODO Linear interpolation - to be revisited
        atm_vol = np.interp(f, strikes, volatilities)
        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEARS

        # Calculate the slope of the volatility curve by taking the end
        # points in the volatilities and strikes to calculate the gradient

        dvol = volatilities[-1] - volatilities[0]
        dk = strikes[-1] - strikes[0]
        b = f * dvol / dk
        var = (atm_vol**2) * np.sqrt(1.0 + 3.0 * t_mat * (b**2))
        return var

    ###########################################################################

    def fair_strike(
        self,
        value_dt,
        stock_price,
        dividend_curve,
        volatility_curve,
        num_call_options,
        num_put_options,
        strike_spacing,
        discount_curve,
        use_forward=True,
    ):
        """Calculate the implied variance according to the volatility surface
        using a static replication methodology with a specially weighted
        portfolio of put and call options across a range of strikes using the
        approximate method set out by Demeterfi et al. 1999."""

        self.num_put_options = num_put_options
        self.num_call_options = num_call_options

        call_type = OptionTypes.EUROPEAN_CALL
        put_type = OptionTypes.EUROPEAN_PUT

        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEARS

        df = discount_curve.df_t(t_mat)
        r = -np.log(df) / t_mat

        dq = dividend_curve.df_t(t_mat)
        q = -np.log(dq) / t_mat

        s0 = stock_price
        g = np.exp((r - q) * t_mat)
        fwd = stock_price * g

        # This fixes the centre strike of the replication options
        if use_forward is True:
            sstar = fwd
        else:
            sstar = stock_price

        # Replication argument from Demeterfi, Derman, Kamal and Zhou from
        # Goldman Sachs Research notes March 1999. See Appendix A. This aim is
        # to use calls and puts to approximate the payoff of a log contract

        min_strike = sstar - (num_put_options + 1) * strike_spacing

        self.put_wts = []
        self.put_strikes = []
        self.call_wts = []
        self.call_strikes = []

        # if the lower strike is < 0 we go to as low as the strike spacing
        if min_strike < strike_spacing:
            k = sstar
            klist = [sstar]
            while k >= strike_spacing:
                k -= strike_spacing
                klist.append(k)
            put_k = np.array(klist)
            self.num_put_options = len(put_k) - 1
        else:
            put_k = np.linspace(sstar, min_strike, num_put_options + 2)

        self.put_strikes = put_k

        max_strike = sstar + (num_call_options + 1) * strike_spacing
        call_k = np.linspace(sstar, max_strike, num_call_options + 2)

        self.call_strikes = call_k

        option_total = (
            2.0 * (r * t_mat - (s0 * g / sstar - 1.0) - np.log(sstar / s0)) / t_mat
        )

        self.call_wts = np.zeros(num_call_options)
        self.put_wts = np.zeros(num_put_options)

        def f(x):
            return (2.0 / t_mat) * ((x - sstar) / sstar - np.log(x / sstar))

        sum_wts = 0.0
        for n in range(0, self.num_put_options):
            kp = put_k[n + 1]
            k = put_k[n]
            self.put_wts[n] = (f(kp) - f(k)) / (k - kp) - sum_wts
            sum_wts += self.put_wts[n]

        sum_wts = 0.0
        for n in range(0, self.num_call_options):
            kp = call_k[n + 1]
            k = call_k[n]
            self.call_wts[n] = (f(kp) - f(k)) / (kp - k) - sum_wts
            sum_wts += self.call_wts[n]

        pi_put = 0.0
        for n in range(0, num_put_options):
            k = put_k[n]
            vol = volatility_curve.volatility(k)
            opt = EquityVanillaOption(self.maturity_dt, k, put_type)
            model = BlackScholes(vol)
            v = opt.value(value_dt, s0, discount_curve, dividend_curve, model)
            pi_put += v * self.put_wts[n]

        pi_call = 0.0
        for n in range(0, num_call_options):
            k = call_k[n]
            vol = volatility_curve.volatility(k)
            opt = EquityVanillaOption(self.maturity_dt, k, call_type)
            model = BlackScholes(vol)
            v = opt.value(value_dt, s0, discount_curve, dividend_curve, model)
            pi_call += v * self.call_wts[n]

        pi = pi_call + pi_put
        option_total += g * pi
        var = option_total

        return var

    ###########################################################################

    def realised_variance(self, close_prices, use_logs=True):
        """Calculate the realised variance according to market standard
        calculations which can either use log or percentage returns."""

        num_observations = len(close_prices)

        for i in range(0, num_observations):
            if close_prices[i] <= 0.0:
                raise FinError("Stock prices must be greater than zero")

        cum_x2 = 0.0

        if use_logs is True:
            for i in range(1, num_observations):
                x = np.log(close_prices[i] / close_prices[i - 1])
                cum_x2 += x * x
        else:
            for i in range(1, num_observations):
                x = (close_prices[i] - close_prices[i - 1]) / close_prices[i - 1]
                cum_x2 += x * x

        var = cum_x2 * 252.0 / num_observations
        return var

    ###########################################################################

    def print_weights(self):
        """Print the list of puts and calls used to replicate the static
        replication component of the variance swap hedge."""

        if self.num_put_options == 0 and self.num_call_options == 0:
            print("No call or put options generated.")
            return

        print("TYPE", "STRIKE", "WEIGHT")
        for n in range(self.num_put_options - 1, -1, -1):
            k = self.put_strikes[n]
            wt = self.put_wts[n] * self.notional
            print("PUT %7.2f %10.3f" % (k, wt))

        for n in range(0, self.num_call_options):
            k = self.call_strikes[n]
            wt = self.call_wts[n] * self.notional
            print("CALL %7.2f %10.3f" % (k, wt))

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.start_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("STRIKE VARIANCE", self.strike_variance)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("PAY STRIKE FLAG", self.pay_strike_flag, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
