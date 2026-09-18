##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

import numpy as np

# TODO: Add perturbatory risk using the analytical methods !!
# TODO: Add Sobol to Monte Carlo

from ...utils.error import FinError

from ...utils.global_types import OptionTypes
from ...utils.global_types import AsianOptionValuationTypes

from ...utils.frequency import FrequencyTypes
from ...utils.helpers import check_argument_types, label_to_string
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve

from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.check_values import check_t_exp
from ...utils.helpers import option_years

from ...models.equity_asian_option_mc import equity_asian_value_mc_fast_cv_numba
from ...models.equity_asian_option_mc import equity_asian_value_mc_fast_numba
from ...models.equity_asian_option_mc import equity_asian_value_mc_numba

from ...models.equity_asian_option_bs import value_curran
from ...models.equity_asian_option_bs import value_turnbull_wakeman
from ...models.equity_asian_option_bs import value_geometric

########################################################################################


########################################################################################


########################################################################################
# An Asian option on an arithmetic average and strike K has a payoff
# Max(SA(T)-K,0) where SA is the arithmetic average
# We define three dates
# - Valuation date for which we want the price
# - Start Averaging Date for when the averaging starts
# - Expiry date for when the payoff is made and the option expires
#
# In the model we have
# tv = is the time now
# t0 = time to the start averaging date in years
# t = time to the expiry date in years
# tau = length of averaging period in years at the start of the option
#
# We can be before the start of the averaging period in which case t0 > 0
# We can be after the start of the averaging period in which case we set t0=0
# and we note that t <= tau
#
# If we are in the averaging period then we need to know the accrued average
# I call this AA and the new average is now given by the accrued average
# The option payoff is now Max( (AA x (tau-t) + SA(t0) x (t-t0))/tau - K,0)
# This simplifies to
#
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#  (1/tau) * Max( (AA x (tau-t) +  - K x tau + SA(t0) x (t-t0)),0)
#
########################################################################################


########################################################################################


class EquityAsianOption:
    """Class for an Equity Asian Option. This is an option with a final payoff
    linked to the averaging of the stock price over some specified period
    before the option expires. The valuation is done for both an arithmetic and
    a geometric average but the former can only be done either using an
    analytical approximation of the arithmetic average distribution or by using
    Monte-Carlo simulation."""

    def __init__(
        self,
        start_averaging_dt: Date,
        expiry_dt: Date,
        strike_price: float,
        opt_type: OptionTypes,
        num_obs: int = 100,
    ):
        """Create an EquityAsian option object which takes a start date for
        the averaging, an expiry date, a strike price, an OPTION_TYPE and a
        number of observations."""

        check_argument_types(self.__init__, locals())

        if start_averaging_dt > expiry_dt:
            raise FinError("Averaging starts after expiry date")

        self.start_averaging_date = start_averaging_dt
        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.opt_type = opt_type
        self.num_observations = num_obs

    ####################################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        method: AsianOptionValuationTypes,
        accrued_average: float = None,
    ):
        """Calculate the value of an Asian option using one of the specified
        analytical approximations for an average rate option. These are the
        three enumerated values in the enum AsianOptionValuationMethods. The
        choices of approximation are (i) GEOMETRIC - the average is a geometric
        one as in paper by Kenna and Worst (1990), (ii) TURNBULL_WAKEMAN -
        this is a value based on an edgeworth expansion of the moments of the
        arithmetic average, and (iii) CURRAN - another approximative approach
        by Curran based on conditioning on the geometric mean price. Just
        choose the corresponding enumerated value to switch between these
        different approaches.

        Note that the accrued average is only required if the value date is
        inside the averaging period for the option."""

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_avg = option_years(value_dt, self.start_averaging_date, fail=False)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        k = self.strike_price
        n = self.num_observations
        opt_type_value = self.opt_type.value

        if method == AsianOptionValuationTypes.GEOMETRIC:
            v = value_geometric(t_avg, t_exp, k, n, opt_type_value, stock_price, r, q, model, accrued_average)

        elif method == AsianOptionValuationTypes.TURNBULL_WAKEMAN:
            v = value_turnbull_wakeman(t_avg, t_exp, k, n, opt_type_value, stock_price, r, q, model, accrued_average)

        elif method == AsianOptionValuationTypes.CURRAN:
            v = value_curran(t_avg, t_exp, k, n, opt_type_value, stock_price, r, q, model, accrued_average)
        else:
            raise FinError("Unknown valuation model")

        return v

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_paths: int,
        seed: int,
        accrued_average: float,
    ):
        """Monte Carlo valuation of the Asian Average option using standard
        Monte Carlo code enhanced by Numba. I have discontinued the use of this
        as it is both slow and has limited variance reduction."""

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_avg = option_years(value_dt, self.start_averaging_date, fail=False)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        tau = t_exp - t_avg

        volatility = model.volatility

        k = self.strike_price
        n = self.num_observations

        v = equity_asian_value_mc_numba(
            t_avg,
            t_exp,
            tau,
            k,
            n,
            self.opt_type.value,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
        )

        return v

    ####################################################################################

    def value_mc_fast(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,  # Model
        num_paths,  # Numpaths integer
        seed,
        accrued_average,
    ):
        """Monte Carlo valuation of the Asian Average option. This method uses
        a lot of Numpy vectorisation. It is also helped by Numba."""

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)
        check_stock_price(stock_price)

        t_exp = option_years(value_dt, self.expiry_dt)
        t_avg = option_years(value_dt, self.start_averaging_date, fail=False)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        tau = t_exp - t_avg

        k = self.strike_price
        n = self.num_observations

        volatility = model.volatility

        v = equity_asian_value_mc_fast_numba(
            t_avg,
            t_exp,
            tau,
            k,
            n,
            self.opt_type.value,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
        )

        return v

    ####################################################################################

    def value_mc_fast_vc_numba(
        self,
        t_avg,
        t_exp,
        stock_price: float,
        r: float,
        q: float,
        model,
        num_paths: int,
        seed: int,
        accrued_average: float,
    ):
        """Monte Carlo valuation of the Asian Average option using a control
        variate method that improves accuracy and reduces the variance of the
        price. This uses Numpy and Numba. This is the standard MC pricer."""

        tau = t_exp - t_avg

        k = self.strike_price
        n = self.num_observations

        volatility = model.volatility

        # For control variate we price a Geometric average option exactly
        v_g_exact = self.value_geometric(
            t_avg,
            t_exp,
            stock_price,
            r,
            q,
            model,
            accrued_average,
        )

        v = equity_asian_value_mc_fast_cv_numba(
            t_avg,
            t_exp,
            tau,
            k,
            n,
            self.opt_type.value,
            stock_price,
            r,
            q,
            volatility,
            num_paths,
            seed,
            accrued_average,
            v_g_exact,
        )

        return v

    ####################################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START AVERAGING DATE", self.start_averaging_date)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("NUM OBSERVATIONS", self.num_observations, "")
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
