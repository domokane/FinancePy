##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO - MUST ADD ACCRUED INTEREST TO MODEL!!!!

from math import exp, sqrt
from typing import List

from numba import njit
import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.math import test_monotonicity
from ...utils.global_vars import g_days_in_year
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types

from ...utils.schedule import Schedule
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes

from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.interpolator import InterpTypes, _uinterpolate

###############################################################################


@njit(fastmath=True, cache=True)
def _value_convertible(
    t_mat,
    face_amount,
    cpn_times,
    cpn_flows,
    call_times,
    call_prices,
    put_times,
    put_prices,
    conv_ratio,
    start_convert_time,
    # Market inputs
    stock_price,
    df_times,
    df_values,
    dividend_times,
    dividend_yields,
    stock_volatility,
    credit_spread,
    recovery_rate,
    # Tree details
    num_steps_per_year,
):

    interp = InterpTypes.FLAT_FWD_RATES.value

    if len(cpn_times) > 0:
        if cpn_times[-1] > t_mat:
            raise FinError("Coupon after maturity")

    if len(call_times) > 0:
        if call_times[-1] > t_mat:
            raise FinError("Call times after maturity")

    if len(put_times) > 0:
        if put_times[-1] > t_mat:
            raise FinError("Put times after maturity")

    if len(df_times) > 0:
        if df_times[-1] > t_mat:
            raise FinError("Discount times after maturity")

    if len(dividend_times) > 0:
        if dividend_times[-1] > t_mat:
            raise FinError("Dividend times after maturity")

    if credit_spread < 0.0:
        raise FinError("Credit spread negative.")

    if recovery_rate < 0.0 or recovery_rate > 1.0:
        raise FinError("Recovery rate should be between 0 and 1.")

    if stock_volatility < 0.0:
        raise FinError("Stock volatility cannot be negative.")

    if num_steps_per_year < 1:
        raise FinError("Num Steps per year must more than 1.")

    if len(dividend_times) > 0.0:
        if dividend_times[-1] > t_mat:
            raise FinError("Last dividend is after bond maturity.")

    if recovery_rate > 0.999 or recovery_rate < 0.0:
        raise FinError("Recovery rate must be between 0 and 0.999.")

    num_times = int(num_steps_per_year * t_mat) + 1  # add one for today time 0
    num_times = num_steps_per_year  # XXXXXXXX!!!!!!!!!!!!!!!!!!!!!

    if num_times < 5:
        raise FinError("Numsteps must be greater than 5.")

    num_levels = num_times

    # this is the size of the step
    dt = t_mat / (num_times - 1)

    tree_times = np.linspace(0.0, t_mat, num_times)
    tree_dfs = np.zeros(num_times)
    for i in range(0, num_times):
        df = _uinterpolate(tree_times[i], df_times, df_values, interp)
        tree_dfs[i] = df

    h = credit_spread / (1.0 - recovery_rate)
    survival_prob = exp(-h * dt)

    # map coupons onto tree but preserve their present value using risky dfs
    tree_flows = np.zeros(num_times)
    num_cpns = len(cpn_times)
    for i in range(0, num_cpns):
        flow_time = cpn_times[i]
        n = int(round(flow_time / dt, 0))
        tree_time = tree_times[n]
        df_flow = _uinterpolate(flow_time, df_times, df_values, interp)
        df_flow *= exp(-h * flow_time)
        df_tree = _uinterpolate(tree_time, df_times, df_values, interp)
        df_tree *= exp(-h * tree_time)
        tree_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree

    # map call onto tree - must have no calls at high value
    tree_call_value = np.ones(num_times) * face_amount * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(round(call_time / dt, 0))
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    tree_put_value = np.zeros(num_times)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(round(put_time / dt, 0))
        tree_put_value[n] = put_prices[i]

    # map discrete dividend yields onto tree dates when they are made
    tree_dividend_yld = np.zeros(num_times)
    numDividends = len(dividend_times)
    for i in range(0, numDividends):
        dividend_time = dividend_times[i]
        n = int(round(dividend_time / dt, 0))
        tree_dividend_yld[n] = dividend_yields[i]

    # Set up the tree of stock prices using a 2D matrix - half the matrix is
    # unused but this may be a cost worth bearing for simpler code. Review.
    tree_stock_value = np.zeros(shape=(num_times, num_levels))
    e = stock_volatility**2 - h
    if e < 0.0:
        raise FinError("Volatility squared minus the hazard rate is negative.")

    u = exp(sqrt(e * dt))
    d = 1.0 / u
    u2 = u * u
    tree_stock_value[0, 0] = stock_price
    for i_time in range(1, num_times):
        s = tree_stock_value[i_time - 1, 0] * d
        tree_stock_value[i_time, 0] = s

        for i_node in range(1, i_time + 1):
            s = s * u2
            tree_stock_value[i_time, i_node] = s

        # we now reduce all stocks by the same yield amount at the same date
        y = tree_dividend_yld[i_time]
        for i_node in range(0, i_time + 1):
            tree_stock_value[i_time, i_node] *= 1.0 - y

    # set up the tree of conversion values. Before allowed to convert the
    # conversion value must be set equal to zero

    tree_convert_value = np.zeros(shape=(num_times, num_levels))
    for i_time in range(0, num_times):
        if tree_times[i_time] >= start_convert_time:
            for i_node in range(0, i_time + 1):
                s = tree_stock_value[i_time, i_node]
                tree_convert_value[i_time, i_node] = s * conv_ratio * 1.0

    #    print_tree(tree_convert_value)

    tree_convert_bond_value = np.zeros(shape=(num_times, num_levels))

    # store probability of up move as a function of time on the tree
    tree_probs_up = np.zeros(num_times)
    tree_probs_dn = np.zeros(num_times)
    q = 0.0  # we have discrete dividends paid as dividend yields only
    for i_time in range(1, num_times):
        a = tree_dfs[i_time - 1] / tree_dfs[i_time] * exp(-q * dt)
        tree_probs_up[i_time] = (a - d * survival_prob) / (u - d)
        tree_probs_dn[i_time] = (u * survival_prob - a) / (u - d)
    #        r = log(a)/dt
    #        n_min = r*r / stock_volatility / stock_volatility

    if np.any(tree_probs_up > 1.0):
        raise FinError("p_up > 1.0. Increase time steps.")

    ###########################################################################
    # work backwards by first setting values at bond maturity date
    ###########################################################################

    flow = tree_flows[num_times - 1]
    bullet_pv = (1.0 + flow) * face_amount
    for i_node in range(0, num_levels):
        convValue = tree_convert_value[num_times - 1, i_node]
        tree_convert_bond_value[num_times - 1, i_node] = max(
            bullet_pv, convValue
        )

    #  begin backward steps from expiry
    for i_time in range(num_times - 2, -1, -1):

        p_up = tree_probs_up[i_time + 1]
        p_dn = tree_probs_dn[i_time + 1]
        pDef = 1.0 - survival_prob
        df = tree_dfs[i_time + 1] / tree_dfs[i_time]
        call = tree_call_value[i_time]
        put = tree_put_value[i_time]
        flow = tree_flows[i_time]

        for i_node in range(0, i_time + 1):
            fut_value_up = tree_convert_bond_value[i_time + 1, i_node + 1]
            fut_value_dn = tree_convert_bond_value[i_time + 1, i_node]
            hold = (
                p_up * fut_value_up + p_dn * fut_value_dn
            )  # p_up already embeds Q
            hold_pv = (
                df * hold
                + pDef * df * recovery_rate * face_amount
                + flow * face_amount
            )
            conv = tree_convert_value[i_time, i_node]
            value = min(max(hold_pv, conv, put), call)
            tree_convert_bond_value[i_time, i_node] = value

        bullet_pv = df * bullet_pv * survival_prob
        bullet_pv += pDef * df * recovery_rate * face_amount
        bullet_pv += flow * face_amount

    price = tree_convert_bond_value[0, 0]
    delta = (tree_convert_bond_value[1, 1] - tree_convert_bond_value[1, 0]) / (
        tree_stock_value[1, 1] - tree_stock_value[1, 0]
    )
    delta_up = (
        tree_convert_bond_value[2, 3] - tree_convert_bond_value[2, 2]
    ) / (tree_stock_value[2, 3] - tree_stock_value[2, 2])
    delta_dn = (
        tree_convert_bond_value[2, 2] - tree_convert_bond_value[2, 1]
    ) / (tree_stock_value[2, 2] - tree_stock_value[2, 1])
    gamma = (delta_up - delta_dn) / (
        tree_stock_value[1, 1] - tree_stock_value[1, 0]
    )
    theta = (tree_convert_bond_value[2, 2] - tree_convert_bond_value[0, 0]) / (
        2.0 * dt
    )
    results = np.array([price, bullet_pv, delta, gamma, theta])
    return results


###############################################################################


class BondConvertible:
    """Class for convertible bonds. These bonds embed rights to call and put
    the bond in return for equity. Until then, they are bullet bonds which
    means they have regular coupon payments of a known size that are paid on
    known dates plus a payment of par at maturity. As the options are price
    based, the decision to convert to equity depends on the stock price,
    the credit quality of the issuer and the level of interest rates."""

    def __init__(
        self,
        maturity_dt: Date,  # bond maturity date
        coupon: float,  # annual coupon
        freq_type: FrequencyTypes,  # coupon frequency type
        start_convert_dt: Date,  # conversion starts on this date
        conversion_ratio: float,  # num shares per face of notional
        call_dts: List[Date],  # list of call dates
        call_prices: List[float],  # list of call prices
        put_dts: List[Date],  # list of put dates
        put_prices: List[float],  # list of put prices
        dc_type: DayCountTypes,  # day count type for accrued
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
    ):
        """Create BondConvertible object by providing the bond Maturity
        date, coupon, frequency type, accrual convention type and then all
        the details regarding the conversion option including the list of the
        call and put dates and the corresponding list of call and put prices.
        """

        check_argument_types(self.__init__, locals())

        if start_convert_dt > maturity_dt:
            raise FinError("Start convert date is after bond maturity.")

        self.maturity_dt = maturity_dt
        self.cpn = coupon
        self.dc_type = dc_type
        self.freq = annual_frequency(freq_type)
        self.freq_type = freq_type
        self.cal_type = cal_type
        self.call_dts = call_dts
        self.call_prices = call_prices

        if len(self.call_dts) != len(self.call_prices):
            raise FinError("Call dates and prices not same length.")

        self.put_dts = put_dts
        self.put_prices = put_prices

        if len(self.put_dts) != len(self.put_prices):
            raise FinError("Put dates and prices not same length.")

        if len(put_dts) > 0:
            if put_dts[-1] > maturity_dt:
                raise FinError("Last put is after bond maturity.")

        if len(call_dts) > 0:
            if call_dts[-1] > maturity_dt:
                raise FinError("Last call is after bond maturity.")

        self.start_convert_dt = start_convert_dt

        if conversion_ratio < 0.0:
            raise FinError("Conversion ratio is negative.")

        self.conversion_ratio = conversion_ratio
        self.par = 100.0

        self.settle_dt = Date(1, 1, 1900)
        """ I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. """

        self.accrued_int = None
        self.accrued_days = 0.0
        self.alpha = 0.0

    ###########################################################################

    def _calculate_cpn_dts(self, settle_dt: Date):
        """Determine the convertible bond cash flow payment dates."""

        # No need to generate flows if settlement date has not changed
        if settle_dt == self.settle_dt:
            return

        self.settle_dt = settle_dt

        bd_type = BusDayAdjustTypes.NONE
        dg_type = DateGenRuleTypes.BACKWARD

        self.cpn_dts = Schedule(
            settle_dt,
            self.maturity_dt,
            self.freq_type,
            self.cal_type,
            bd_type,
            dg_type,
        ).generate()

        self.pcd = self.cpn_dts[0]
        self.ncd = self.cpn_dts[1]

        self.accrued_int = None
        self.accrued_interest(settle_dt, 1.0)

    ###########################################################################

    def value(
        self,
        settle_dt: Date,
        stock_price: float,
        stock_volatility: float,
        dividend_dts: List[Date],
        dividend_yields: List[float],
        discount_curve: DiscountCurve,
        credit_spread: float,
        recovery_rate: float = 0.40,
        num_steps_per_year: int = 100,
    ):
        """
        A binomial tree valuation model for a convertible bond that captures
        the embedded equity option due to the existence of a conversion option
        which can be invoked after a specific date.

        The model allows the user to enter a schedule of dividend payment
        dates but the size of the payments must be in yield terms i.e. a known
        percentage of currently unknown future stock price is paid. Not a
        fixed amount. A fixed yield. Following this payment the stock is
        assumed to drop by the size of the dividend payment.

        The model also captures the stock dependent credit risk of the cash
        flows in which the bond price can default at any time with a hazard
        rate implied by the credit spread and an associated recovery rate.
        This is the model proposed by Hull (OFODS 6th edition,.page 522).

        The model captures both the issuer's call schedule which is assumed
        to apply on a list of dates provided by the user, along with a call
        price. It also captures the embedded owner's put schedule of prices.
        """

        if stock_price <= 0.0:
            stock_price = 1e-10  # Avoid overflows in delta calc

        if stock_volatility <= 0.0:
            stock_volatility = 1e-10  # Avoid overflows in delta calc

        self._calculate_cpn_dts(settle_dt)

        t_mat = (self.maturity_dt - settle_dt) / g_days_in_year

        if t_mat <= 0.0:
            raise FinError("Maturity must not be on or before the value date.")

        # We include time zero in the coupon times and flows
        cpn_times = [0.0]
        cpn_flows = [0.0]

        cpn = self.cpn / self.freq

        for dt in self.cpn_dts[1:]:
            flow_time = (dt - settle_dt) / g_days_in_year
            cpn_times.append(flow_time)
            cpn_flows.append(cpn)

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        if np.any(cpn_times < 0.0):
            raise FinError("No coupon times can be before the value date.")

        if np.any(cpn_times > t_mat):
            raise FinError("No coupon times can be after the maturity date.")

        call_times = []

        for dt in self.call_dts:
            call_time = (dt - settle_dt) / g_days_in_year
            call_times.append(call_time)

        call_times = np.array(call_times)
        call_prices = np.array(self.call_prices)

        if np.any(call_times < 0.0):
            raise FinError("No call times can be before the value date.")

        if np.any(call_times > t_mat):
            raise FinError("No call times can be after the maturity date.")

        put_times = []

        for dt in self.put_dts:
            put_time = (dt - settle_dt) / g_days_in_year
            put_times.append(put_time)

        put_times = np.array(put_times)
        put_prices = np.array(self.put_prices)

        if np.any(put_times > t_mat):
            raise FinError("No put times can be after the maturity date.")

        if np.any(put_times <= 0.0):
            raise FinError("No put times can be on or before value date.")

        if len(dividend_yields) != len(dividend_dts):
            raise FinError("Number of dividend yields and dates not same.")

        dividend_times = []
        for dt in dividend_dts:
            dividend_time = (dt - settle_dt) / g_days_in_year
            dividend_times.append(dividend_time)
        dividend_times = np.array(dividend_times)
        dividend_yields = np.array(dividend_yields)

        # If it's before today it starts today
        tconv = (self.start_convert_dt - settle_dt) / g_days_in_year
        tconv = max(tconv, 0.0)

        discount_factors = []
        for t in cpn_times:
            df = discount_curve.df_t(t)
            discount_factors.append(df)

        discount_times = np.array(cpn_times)
        discount_factors = np.array(discount_factors)

        if test_monotonicity(cpn_times) is False:
            raise FinError("Coupon times not monotonic")

        if test_monotonicity(call_times) is False:
            raise FinError("Coupon times not monotonic")

        if test_monotonicity(put_times) is False:
            raise FinError("Coupon times not monotonic")

        if test_monotonicity(discount_times) is False:
            raise FinError("Coupon times not monotonic")

        if test_monotonicity(dividend_times) is False:
            raise FinError("Coupon times not monotonic")

        v1 = _value_convertible(
            t_mat,
            self.par,
            cpn_times,
            cpn_flows,
            call_times,
            call_prices,
            put_times,
            put_prices,
            self.conversion_ratio,
            tconv,
            # Market inputs
            stock_price,
            discount_times,
            discount_factors,
            dividend_times,
            dividend_yields,
            stock_volatility,
            credit_spread,
            recovery_rate,
            # Tree details
            num_steps_per_year,
        )

        v2 = _value_convertible(
            t_mat,
            self.par,
            cpn_times,
            cpn_flows,
            call_times,
            call_prices,
            put_times,
            put_prices,
            self.conversion_ratio,
            tconv,
            # Market inputs
            stock_price,
            discount_times,
            discount_factors,
            dividend_times,
            dividend_yields,
            stock_volatility,
            credit_spread,
            recovery_rate,
            # Tree details
            num_steps_per_year + 1,
        )

        cbprice = (v1[0] + v2[0]) / 2.0
        bond = (v1[1] + v2[1]) / 2.0
        delta = (v1[2] + v2[2]) / 2.0
        gamma = (v1[3] + v2[3]) / 2.0
        theta = (v1[4] + v2[4]) / 2.0

        results = {
            "cbprice": cbprice,
            "bond": bond,
            "delta": delta,
            "gamma": gamma,
            "theta": theta,
        }

        return results

    ###########################################################################

    def accrued_days(self, settle_dt: Date):
        """Calculate number days from previous coupon date to settlement."""
        self._calculate_cpn_dts(settle_dt)

        if len(self.cpn_dts) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        return settle_dt - self.pcd

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float):
        """Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date."""

        if settle_dt != self.settle_dt:
            self._calculate_cpn_dts(settle_dt)

        if len(self.cpn_dts) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self.dc_type)

        (acc_factor, num, _) = dc.year_frac(
            self.pcd, settle_dt, self.ncd, self.freq
        )

        self.alpha = 1.0 - acc_factor * self.freq

        self.accrued_int = acc_factor * face * self.cpn
        self.accrued_days = num
        return self.accrued_int

    ###########################################################################

    def current_yield(self, clean_price: float):
        """Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)"""

        y = self.cpn * self.par / clean_price
        return y

    ###########################################################################

    def __repr__(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("COUPON", self.cpn)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        s += label_to_string("CONVERSION RATIO", self.conversion_ratio)
        s += label_to_string("START CONVERT DATE", self.start_convert_dt)
        s += label_to_string("CALL", "DATES")

        for i in range(0, len(self.call_dts)):
            s += label_to_string(self.call_dts[i], self.call_prices[i])

        s += label_to_string("PUT", "DATES")

        for i in range(0, len(self.put_dts)):
            s += label_to_string(self.put_dts[i], self.put_prices[i])

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
###############################################################################
# TEST PV OF CASHFLOW MAPPING
#    if 1==0:
#        pv = 0.0
#        for i in range(0, num_cpns):
#            t = cpn_times[i]
#            df = uinterpolate(t, discount_times, discount_factors, interp)
#            pv += df * cpn_amounts[i]
#            print(i, t, cpn_amounts[i], df, pv)
#        pv += df
#
#        print("ACTUAL PV",pv)
#
#        pv = 0.0
#        for i in range(0, num_times):
#            t = tree_times[i]
#            df = uinterpolate(t, discount_times, discount_factors, interp)
#            pv += df * tree_flows[i]
#            print(i, t, tree_flows[i], df, pv)
#        pv += df
#
#        print("ACTUAL PV",pv)
###############################################################################
###############################################################################
###############################################################################


def print_tree(array):
    n1, n2 = array.shape
    for i in range(0, n1):
        for j in range(0, n2):
            x = array[j, n1 - 1 - i]
            if x != 0.0:
                print("%10.2f" % array[j, n1 - i - 1], end="")
            else:
                print("%10s" % "-", end="")
        print("")


###############################################################################
