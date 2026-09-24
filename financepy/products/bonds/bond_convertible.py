##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO - MUST ADD ACCRUED INTEREST TO MODEL!!!!

from typing import List

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.math import test_monotonicity
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types

from ...utils.schedule import Schedule
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.calendar import Calendar

from ...market.curves.discount_curve import DiscountCurve
from ...models.bond_convertible_bs_tree import value_convertible
from ...utils.check_values import check_curve_dt

########################################################################################


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
        call_prices: np.ndarray,  # list of call prices
        put_dts: List[Date],  # list of put dates
        put_prices: np.ndarray,  # list of put prices
        accrual_dc_type: DayCountTypes,  # day count type for accrued
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
    ) -> None:
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
        self.accrual_dc_type = accrual_dc_type
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

        self._pcd = None
        self._ncd = None

        self.cal_type = cal_type

        self.cpn_dts = []
        self.payment_dts = []

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

        self._pcd = self.cpn_dts[0]
        self._ncd = self.cpn_dts[1]

        calendar = Calendar(self.cal_type)

        self.payment_dts = []

        # I do not adjust the first date as it is the issue date
        self.payment_dts.append(self.cpn_dts[0])

        for cpn_dt in self.cpn_dts[1:]:
            pmt_dt = calendar.adjust(cpn_dt, bd_type)
            self.payment_dts.append(pmt_dt)

        self.accrued_int = None
        self.accrued_interest(settle_dt, 1.0)

    ###########################################################################

    def value(
        self,
        settle_dt: Date,
        stock_price: float,
        stock_volatility: float,
        dividend_dts: List[Date],
        dividend_yields: np.ndarray,
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

        check_curve_dt(settle_dt, discount_curve)

        if stock_price <= 0.0:
            stock_price = 1e-10  # Avoid overflows in delta calc

        if stock_volatility <= 0.0:
            stock_volatility = 1e-10  # Avoid overflows in delta calc

        self._calculate_cpn_dts(settle_dt)

        t_mat = (self.maturity_dt - settle_dt) / G_DAYS_IN_YEAR

        if t_mat <= 0.0:
            raise FinError("Maturity must not be on or before the value date.")

        # We include time zero in the coupon times and flows
        cpn_times = [0.0]
        cpn_flows = [0.0]

        cpn = self.cpn / self.freq

        for dt in self.payment_dts[1:]:
            flow_time = (dt - settle_dt) / G_DAYS_IN_YEAR
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
            call_time = (dt - settle_dt) / G_DAYS_IN_YEAR
            call_times.append(call_time)

        call_times = np.array(call_times)
        call_prices = np.array(self.call_prices)

        if np.any(call_times < 0.0):
            raise FinError("No call times can be before the value date.")

        if np.any(call_times > t_mat):
            raise FinError("No call times can be after the maturity date.")

        put_times = []

        for dt in self.put_dts:
            put_time = (dt - settle_dt) / G_DAYS_IN_YEAR
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
            dividend_time = (dt - settle_dt) / G_DAYS_IN_YEAR
            dividend_times.append(dividend_time)
        dividend_times = np.array(dividend_times)
        dividend_yields = np.array(dividend_yields)

        # If it's before today it starts today
        tconv = (self.start_convert_dt - settle_dt) / G_DAYS_IN_YEAR
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

        v1 = value_convertible(
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

        v2 = value_convertible(
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

    # def accrued_days(self, settle_dt: Date):
    #     """Calculate number days from previous coupon date to settlement."""
    #     self._calculate_cpn_dts(settle_dt)

    #     if len(self.cpn_dts) <= 2:
    #         raise FinError("Accrued interest - not enough flow dates.")

    #     return settle_dt - self._pcd

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float):
        """Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date."""

        if settle_dt != self.settle_dt:
            self._calculate_cpn_dts(settle_dt)

        if len(self.cpn_dts) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self.accrual_dc_type)

        acc_factor, num, _ = dc.year_frac(self._pcd, settle_dt, self._ncd, self.freq)

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
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("MATURITY_DATE", self.maturity_dt)
        s += label_to_string("COUPON", self.cpn)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DC_TYPE", self.accrual_dc_type)
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


########################################################################################
########################################################################################
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
########################################################################################
########################################################################################
########################################################################################


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


########################################################################################
