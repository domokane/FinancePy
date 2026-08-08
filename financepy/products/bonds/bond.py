########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

########################################################################################
# TODO: - ROUNDING CONVENTIONS FOR ACCRUED
# TODO: - CHECK OAS CALCULATION
# TODO:  - Check how first coupon on floating leg is sized on asset swaps. """
########################################################################################

# https://www.dmo.gov.uk/media/15004/convention_changes.pdf

########################################################################################
# Conventions:
#  GILTS - SEMI ANNUAL ACT/ACT
#  US TREASURIES
########################################################################################

########################################################################################
# NOTE THAT I ASSUME THAT IF YOU SETTLE A SWAP ON A COUPON PAYMENT DATE YOU
# GET THE COUPON AND THE ACCRUED INTEREST EQUALS THE COUPON.
########################################################################################

from typing import List
from enum import Enum
import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.global_vars import G_SMALL
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.schedule import Schedule
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types
from ...utils.helpers import times_from_dates
from ...utils.math import npv
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.interpolator import InterpTypes
from ...market.curves.pwf_onf_discount_curve import PWFONFDiscountCurve
from ...market.curves.composite_discount_curve import CompositeDiscountCurve
from ...market.curves.bond_bootstrap_discount_curve import BondBootstrapDiscountCurve

# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

########################################################################################


class CouponType(Enum):
    FIXED = 0  # USES C/F
    ACCRUED = 1  # USES Accrual_Period * C


########################################################################################


class YTMCalcType(Enum):
    ZERO = 0
    UK_DMO = 1
    US_STREET = 2
    US_TREASURY = 3
    CFETS = 4  # China Foreign Exchange Trade System
    CALCULUS = 5  # Using Calculus for duration


########################################################################################


def _f(ytm, *args):
    """Function used to do root search in price to yield calculation."""
    bond = args[0]
    settle_dt = args[1]
    price = args[2]
    convention = args[3]
    px = bond.dirty_price_from_ytm(settle_dt, ytm, convention)
    obj_fn = px - price
    return obj_fn


###############################################################################


def _g(oas, *args):
    """Function used to do root search in price to OAS calculation."""
    bond = args[0]
    settle_dt = args[1]
    price = args[2]
    discount_curve = args[3]
    px = bond.dirty_price_from_oas(settle_dt, discount_curve, oas)
    obj_fn = px - price
    return obj_fn


###############################################################################


def validate_yield(ytm):

    ytms = []

    if isinstance(ytm, float) or isinstance(ytm, np.float64):
        ytms = np.array([ytm])
    elif isinstance(ytm, list) or isinstance(ytm, np.ndarray):
        ytms = np.array(ytm)

    for ytm in ytms:
        if ytm < -0.20 or ytm > 2.0:
            raise FinError("YTM of " + str(ytm) + " must be between -20% and +200%")


###############################################################################


def validate_price(px):
    if px < 0.0 or px > 300.0:
        raise FinError("Price of " + str(px) + " must be between 0 and 300")


###############################################################################


def vectorise_price(price):

    if isinstance(price, int):
        price = float(price)

    if isinstance(price, float) or isinstance(price, np.float64):
        prices = np.array([price])
    elif isinstance(price, list) or isinstance(price, np.ndarray):
        prices = np.array(price)
    else:
        raise FinError("Unknown type for price " + str(type(price)))

    return prices


########################################################################################


class Bond:
    """Class for fixed coupon bonds and performing related analytics. These
    are bullet bonds which means they have regular coupon payments of a known
    size that are paid on known dts plus a payment of par at maturity."""

    def __init__(
        self,
        issue_dt: Date,
        maturity_dt: Date,
        coupon: float,  # Annualised bond coupon
        freq_type: FrequencyTypes,
        accrual_dc_type: DayCountTypes,
        ex_div_days: int = 0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type=BusDayAdjustTypes.FOLLOWING,
        dg_type=DateGenRuleTypes.BACKWARD,
        cpn_type=CouponType.FIXED,
    ):
        """Create Bond object by providing the issue date, maturity Date,
        coupon frequency, annualised coupon, the accrual convention type, face
        amount and the number of ex-dividend days. A calendar type is used
        to determine holidays from which coupon dts might be shifted."""

        check_argument_types(self.__init__, locals())

        if issue_dt >= maturity_dt:
            raise FinError("Issue Date must preceded maturity date.")

        # If the maturity date falls on the last day of the month we assume
        # that earlier flows also fall on month ends
        self.end_of_month = False
        if maturity_dt.is_eom():
            self.end_of_month = True

        self.issue_dt = issue_dt
        self.maturity_dt = maturity_dt

        if coupon == 0.0:
            raise FinError("Zero coupon bonds must use BondZero class.")

        if freq_type == FrequencyTypes.ZERO:
            raise FinError("Zero coupon bonds must use BondZero class.")

        self.cpn = coupon
        self.cpn_type = cpn_type
        self.freq_type = freq_type

        self.accrual_dc_type = accrual_dc_type
        self.freq = annual_frequency(freq_type)

        if ex_div_days > 90:
            raise FinError("Ex dividend days cannot be more than 90" + str(ex_div_days))

        self.ex_div_dt = None
        self.ex_div_days = ex_div_days

        self.par = 100.0  # This is how price is quoted and amount at maturity
        self.cal_type = cal_type

        self.cpn_dts = []  # can be holidays or weekend
        self.payment_dts = []  # Actual pay dts are adjusted to bus days
        self.flow_amounts = []

        self.bd_type = bd_type
        self.dg_type = dg_type

        self._calculate_cpn_dts()
        self._calculate_payment_dts()
        self._calculate_flow_amounts()

        self._pcd = None
        self._ncd = None

        # Private
        self.accrued_int = None
        self.accrued_days = 0.0
        self.alpha = 0.0

    ####################################################################################

    def _calculate_cpn_dts(self):
        """Determine the unadjusted bond coupon dts from the issue date to the
        maturity date even if they have passed. Note that for analytical
        calculations these are not usually adjusted and so may fall on a
        weekend or holiday.
        """

        # This should only be called once from init
        # bd_type = BusDayAdjustTypes.FOLLOWING
        # dg_type = DateGenRuleTypes.BACKWARD

        self.cpn_dts = Schedule(
            self.issue_dt,
            self.maturity_dt,
            self.freq_type,
            CalendarTypes.NONE,
            self.bd_type,
            self.dg_type,
            end_of_month=self.end_of_month,
        ).generate()

    ####################################################################################

    def _calculate_payment_dts(self):
        """For the actual payment dts, they are adjusted,
        and so we then use the calendar payment dts. Although payments are
        calculated as though coupon periods are the same length, payments that
        fall on a Saturday or Sunday can only be made on the next business day
        """

        # dts are adjusted forward to the next business day
        bus_day_adj_type = BusDayAdjustTypes.FOLLOWING
        calendar = Calendar(self.cal_type)

        # Expect at least an issue date and a maturity date - if not - problem
        if len(self.cpn_dts) < 2:
            raise FinError("Cannot calculate payment dts with one payment")

        self.payment_dts = []

        # I do not adjust the first date as it is the issue date
        self.payment_dts.append(self.cpn_dts[0])

        for cpn_dt in self.cpn_dts[1:]:
            pmt_dt = calendar.adjust(cpn_dt, bus_day_adj_type)
            self.payment_dts.append(pmt_dt)

    ############################################################################

    def times(self, settle_dt, time_dc_type:DayCountTypes = DayCountTypes.ACT_365F):
        """Years from settlement to payments using day count convention"""
        times = times_from_dates(settle_dt, self.payment_dts, time_dc_type)

        print("TIMES")
        print(times)
        print(len(times))
        return times

    ############################################################################

    def flows(self, settle_dt):
        """Times from settlement to payments"""
        n_flows = len(self.flow_amounts)
        flows = []
        for i in range(0, n_flows):
            flow = self.flow_amounts[i]
            if self.payment_dts[i] >= settle_dt:
                flows.append(flow)

        print("Flows")
        print(flows)
        print(len(flows))
        return flows

    ###########################################################################

    def _calculate_flow_amounts(self):
        """Determine the bond cash flow payment amounts without principal.
        There is no adjustment based on the adjusted payment dts."""

        self.flow_amounts = [0.0]

        if self.cpn_type == CouponType.FIXED:

            for _ in self.cpn_dts[1:]:

                flow = self.cpn / self.freq
                self.flow_amounts.append(flow)

        elif self.cpn_type == CouponType.ACCRUED:

            self_dc = DayCount(self.accrual_dc_type)
            prev_dt = self.cpn_dts[0]

            for next_dt in self.cpn_dts[1:]:

                yf, _, _ = self_dc.year_frac(prev_dt, next_dt)
                prev_dt = next_dt
                flow = yf * self.cpn
                self.flow_amounts.append(flow)

        else:
            raise FinError("Error: Unknown Coupon Type")

        return

    ###########################################################################

    def reset_flows(
        self,
        cpn_dts: List[Date],
        payment_dts: List[Date],
        flow_amounts: np.ndarray,
    ):
        """Set the flows of the bond externally. Coupon dates are for accrued
        while payment dates are calendar adjusted. Flows are on payment dates
        """

        n_cpn_dts = len(cpn_dts)
        n_payment_dts = len(payment_dts)
        n_flows = len(flow_amounts)

        if n_cpn_dts != n_payment_dts:
            raise FinError("Number of coupon dates not equal to number payments")

        if n_cpn_dts != n_flows:
            raise FinError("Number of coupon dates not equal to number flows")

        for i in range(1, n_cpn_dts):

            if cpn_dts[i] < cpn_dts[i - 1]:
                raise FinError("Coupon dates not in order")

            if payment_dts[i] < payment_dts[i - 1]:
                raise FinError("Payment dates not in order")

        self.cpn_dts = cpn_dts
        self.payment_dts = payment_dts
        self.flow_amounts = flow_amounts

    ###########################################################################

    def dirty_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the dirty price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM."""

        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settlement is after maturity date")

        if convention not in YTMCalcType:
            raise FinError("Yield convention unknown." + str(convention))

        # We MUST call this to update value of self.alpha and ex-div date
        self.accrued_interest(settle_dt, 1.0)

        #######################################################################
        # HANDLE EX_DIVIDEND DATES
        #######################################################################

        pay_first_cpn = 1.0
        if settle_dt > self.ex_div_dt:
            pay_first_cpn = 0.0

        #######################################################################

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        f = annual_frequency(self.freq_type)
        c = self.cpn

        if convention == YTMCalcType.ZERO:
            raise FinError("Zero coupon bonds must use BondZero class.")

        v = 1.0 / (1.0 + ytm / f)

        # n is the number of flows after the next coupon date
        # Do not use payment dates here as convention is to use coupon dates
        n = 0
        for dt in self.cpn_dts:
            if dt > settle_dt:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No coupons left")

        if convention == YTMCalcType.UK_DMO:
            if n == 0:
                dp = (v ** (self.alpha)) * (1.0 + pay_first_cpn * c / f)
            else:
                term1 = (c / f) * pay_first_cpn
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = v**n
                dp = (v ** (self.alpha)) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.US_TREASURY:
            if n == 0:
                dp = (v ** (self.alpha)) * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_cpn
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = v**n
                vw = 1.0 / (1.0 + self.alpha * ytm / f)
                dp = (vw) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.US_STREET:
            if n == 0:
                vw = 1.0 / (1.0 + self.alpha * ytm / f)
                dp = vw * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_cpn
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = v**n
                dp = (v ** (self.alpha)) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.CFETS:
            if n == 0:
                last_year = self.maturity_dt.add_tenor("-12M")

                dc = DayCount(DayCountTypes.ACT_365L)

                alpha = (
                    1
                    - dc.year_frac(
                        last_year,
                        settle_dt,
                        self.maturity_dt,
                        freq_type=FrequencyTypes.ANNUAL,
                    )[0]
                )

                vw = 1.0 / (1.0 + alpha * ytm)
                dp = vw * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_cpn
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = v**n
                dp = (v ** (self.alpha)) * (term1 + term2 + term3 + term4)

        elif convention == YTMCalcType.CALCULUS:

            dp = 0.0
            d = 1.0 + ytm / f
            # g starts at the discount factor for the NEXT coupon date
            g = 1.0 / np.pow(d, self.alpha)
            flow = self.cpn / self.freq

            n_next = 0
            for dt in self.cpn_dts:
                if dt > settle_dt:
                    break
                n_next += 1

            n_start = n_next

            if settle_dt > self.ex_div_dt:
                n_start = n_next + 1
                g = g / d

            last_g = 0.0
            # n represents 'periods from today'
            for dt in self.cpn_dts[n_start:]:
                dp += flow * g
                last_g = g
                g = g / d

            dp += 1.0 * last_g

        else:
            raise FinError("Unknown yield convention")

        return dp * self.par

    ###########################################################################

    def dirty_price_from_ytm_vector(
        self,
        settle_dts: (Date, list),
        ytms: (float, list, np.ndarray),
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """
        Calculate dirty prices for vectors of settlement dates or yields.
        """
        # Standardize inputs to lists
        if isinstance(settle_dts, list) or isinstance(settle_dts, np.ndarray):
            settle_dt_list = settle_dts
        else:
            settle_dt_list = [settle_dts]

        if isinstance(ytms, list) or isinstance(ytms, np.ndarray):
            ytm_list = ytms
        else:
            ytm_list = [ytms]

        validate_yield(ytms)

        num_dates = len(settle_dt_list)
        num_ytms = len(ytm_list)

        # Create vectors of equal length based on broadcasting rules
        settle_dt_vector = []
        ytm_vector = []

        if num_dates == 1 and num_ytms > 1:
            # One date matched with many yields
            for i in range(num_ytms):
                settle_dt_vector.append(settle_dt_list[0])
                ytm_vector.append(ytm_list[i])

        elif num_ytms == 1 and num_dates > 1:
            # One yield matched with many dates
            for i in range(num_dates):
                settle_dt_vector.append(settle_dt_list[i])
                ytm_vector.append(ytm_list[0])

        elif num_dates == num_ytms:
            # Matched pairs of dates and yields
            for i in range(num_dates):
                settle_dt_vector.append(settle_dt_list[i])
                ytm_vector.append(ytm_list[i])

        else:
            raise FinError(
                "Number of dates and yields must match, or one must be a scalar."
            )

        # Execute the loop using the standardized vectors
        dps = []
        num_to_calculate = len(settle_dt_vector)

        for i in range(num_to_calculate):
            s_dt = settle_dt_vector[i]
            y = ytm_vector[i]

            # Call the core pricing function
            price = self.dirty_price_from_ytm(s_dt, y, convention)

            dps.append(price)

        if num_to_calculate == 1:
            return dps[0]
        else:
            return np.array(dps)

    ###########################################################################

    def forward_price(
        self,
        settle_dt: Date,
        forward_dt: Date,
        clean_price: float,
        repo_rate: float,
        repo_dc_type: DayCountTypes = DayCountTypes.ACT_360,
    ):

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        if clean_price < 0 or repo_rate < 0:
            raise ValueError("Prices and repo rate must be non-negative")

        if settle_dt >= forward_dt:
            raise ValueError("Settlement date must be before forward date")

        repo_dc = DayCount(repo_dc_type)
        t_fwd, _, _ = repo_dc.year_frac(settle_dt, forward_dt)

        accrued_settle = self.accrued_interest(settle_dt)
        accrued_forward = self.accrued_interest(forward_dt)

        fv_cpns = 0.0
        for dt, amt in zip(self.payment_dts[1:], self.flow_amounts[1:]):
            if settle_dt < dt <= forward_dt:
                t_cpn_to_forward, _, _ = repo_dc.year_frac(dt, forward_dt)
                g_cpn = 1.0 + repo_rate * t_cpn_to_forward
                fv_cpns += amt * self.par * g_cpn

        dirty_price = clean_price + accrued_settle
        g = 1.0 + repo_rate * t_fwd
        fwd_dirty_price = dirty_price * g - fv_cpns
        fwd_clean_price = fwd_dirty_price - accrued_forward
        return fwd_clean_price

    ###########################################################################

    def principal(
        self, settle_dt: Date, ytm: float, face: float, convention: YTMCalcType
    ):
        """Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates."""

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        dirty_price = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        principal = dirty_price * face / self.par
        principal = principal - self.accrued_int
        return principal

    ###########################################################################

    def dollar_duration(
        self,
        settle_dt: (Date, List),
        ytm: (float, List),
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate dollar duration from a 1bp bumped-yield DV01.

        Dollar duration is the dirty-price change for a 100bp change in yield,
        approximated as:

            dollar_duration = 10000 * DV01

        where DV01 is the positive dirty-price change for a +1bp yield shift.
        """

        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        dd = 10000 * self.dv01(settle_dt, ytm, convention)
        return dd

    ###########################################################################

    def dv01(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """
        Change in dirty price for a +1 basis point increase in yield.
        Bloomberg-style DV01 / PV01 is positive for a long fixed-rate bond.
        """

        # Change in price for a
        validate_yield(ytm)

        dy = 0.0001  # 1 bp
        p_0 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p_1 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        dP = p_0 - p_1
        return dP

    ###########################################################################

    def macaulay_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the Macaulay duration of the bond on a settlement date
        given its yield to maturity."""

        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        d = 1.0 + ytm / self.freq

        if convention is not YTMCalcType.CALCULUS:
            dd = self.dollar_duration(settle_dt, ytm, convention)
            dp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
            md = (dd / dp) * d

        else:

            # We call this to update value of self.alpha and ex-dividend date
            self.accrued_interest(settle_dt, 1.0)

            md = 0.0
            dp = 0.0
            # g starts at the discount factor for the NEXT coupon date
            g = 1.0 / np.pow(d, self.alpha)
            flow = self.cpn / self.freq

            # 1. Find the index of the next coupon date
            n_next = 0
            for dt in self.cpn_dts:
                if dt > settle_dt:
                    break
                n_next += 1

            n_start = n_next
            k = 0

            if settle_dt > self.ex_div_dt:
                n_start = n_next + 1
                # Move the discount factor and the period counter forward by one
                g = g / d
                k = 1

            # n represents 'periods from today'
            for dt in self.cpn_dts[n_start:]:
                t = (self.alpha + k) / self.freq

                dp += flow * g
                md += flow * t * g

                last_g = g
                last_t = t

                g = g / d
                k += 1

            dp += 1.0 * last_g
            md += 1.0 * last_t * last_g

            dp = dp * self.par
            md = md * self.par

            md = md / dp

        return md

    ###########################################################################

    def modified_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the modified duration of the bond on a settlement
        date given its yield to maturity."""
        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        dd = self.dollar_duration(settle_dt, ytm, convention)
        dp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd / dp
        return md

    ###########################################################################

    def theta(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
        day_step: int = 1,
    ):
        """Calculate the change in the dirty price due to moving calendar time
        i.e. the settlement date, forward by day_step days."""
        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        fp1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        next_dt = settle_dt.add_days(day_step)
        fp2 = self.dirty_price_from_ytm(next_dt, ytm, convention)

        theta = fp2 - fp1
        return theta

    ###########################################################################

    def key_rate_durations_zero_tent(
        self,
        settle_dt: Date,
        zero_curve,
        key_rate_tenors: list = None,
        shift: float = None,
    ):
        """
        Calculate key-rate durations using the textbook zero-rate tent formulation.

        This method perturbs the zero curve directly with piecewise-linear
        tent functions:
            - interior key rates have triangular weights
            - the first key rate is flat to the left
            - the last key rate is flat to the right

        Parameters
        ----------
        settle_dt : Date
            Settlement date.
        zero_curve : DiscountCurve
            Base zero/discount curve used for valuation. It must provide
            discount factors through df(date).
        key_rate_tenors : list of float, optional
            Key-rate maturities in years. If None, defaults to
            [0.5, 1, 2, 3, 5, 7, 10, 20, 30].
        shift : float, optional
            Size of the zero-rate bump. If None, defaults to 1 bp = 0.0001.

        Returns
        -------
        tuple
            (key_rate_tenors, key_rate_durations)

        Notes
        -----
        The bumped zero curve is defined by

            z_i^{up/down}(t) = z(t) +/- shift * w_i(t)

        where w_i(t) is the tent weight for key node i.

        This implementation assumes continuously-compounded zero rates:
            z(t) = -log(df(t))/t

        For t = 0, the bump is taken to be zero.
        """

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        if key_rate_tenors is None:
            key_rate_tenors = np.array([0.5, 1, 2, 3, 5, 7, 10, 20, 30], dtype=float)
        else:
            key_rate_tenors = np.array(key_rate_tenors, dtype=float)

        if np.any(np.diff(key_rate_tenors) <= 0.0):
            raise FinError("key_rate_tenors must be strictly increasing")

        if shift is None:
            shift = 1.0e-4

        if shift <= 0.0:
            raise FinError("shift must be positive")

        # ------------------------------------------------------------------
        # Helper: textbook key-rate tent weights
        # ------------------------------------------------------------------
        def tent_weight(t, i, tenors):
            n = len(tenors)

            # first node: flat to the left, linearly down to next node
            if i == 0:
                t1 = tenors[0]
                t2 = tenors[1]
                if t <= t1:
                    return 1.0
                elif t1 < t < t2:
                    return (t2 - t) / (t2 - t1)
                else:
                    return 0.0

            # last node: linearly up from previous node, flat to the right
            elif i == n - 1:
                tnm1 = tenors[n - 2]
                tn = tenors[n - 1]
                if t <= tnm1:
                    return 0.0
                elif tnm1 < t < tn:
                    return (t - tnm1) / (tn - tnm1)
                else:
                    return 1.0

            # interior nodes: standard triangular tent
            else:
                tL = tenors[i - 1]
                tC = tenors[i]
                tR = tenors[i + 1]

                if t <= tL or t >= tR:
                    return 0.0
                elif tL < t <= tC:
                    return (t - tL) / (tC - tL)
                else:
                    return (tR - t) / (tR - tC)

        # ------------------------------------------------------------------
        # Helper: convert discount factor to continuous zero rate
        # ------------------------------------------------------------------
        def zero_rate_from_df(df, t):
            if t <= 0.0:
                return 0.0
            if df <= 0.0:
                raise FinError("Discount factor must be positive")
            return -np.log(df) / t

        # ------------------------------------------------------------------
        # Helper: build a bumped discount curve wrapper
        # ------------------------------------------------------------------
        class BumpedZeroTentCurve:
            def __init__(self, base_curve, settle_dt, tenors, node_index, bump):
                self.base_curve = base_curve
                self.settle_dt = settle_dt
                self.tenors = tenors
                self.node_index = node_index
                self.bump = bump

            def df(self, dt):
                t = (dt - self.settle_dt) / G_DAYS_IN_YEAR
                if t <= 0.0:
                    return 1.0

                base_df = self.base_curve.df(dt)
                z = zero_rate_from_df(base_df, t)
                w = tent_weight(t, self.node_index, self.tenors)
                z_bumped = z + self.bump * w
                return np.exp(-z_bumped * t)

        # ------------------------------------------------------------------
        # Base price
        # ------------------------------------------------------------------
        p0 = self.dirty_price_from_discount_curve(settle_dt, zero_curve)

        krds = []

        for i in range(len(key_rate_tenors)):
            curve_up = BumpedZeroTentCurve(
                zero_curve, settle_dt, key_rate_tenors, i, +shift
            )
            curve_dn = BumpedZeroTentCurve(
                zero_curve, settle_dt, key_rate_tenors, i, -shift
            )

            p_up = self.dirty_price_from_discount_curve(settle_dt, curve_up)
            p_dn = self.dirty_price_from_discount_curve(settle_dt, curve_dn)

            krd = (p_dn - p_up) / (2.0 * shift * p0)
            krds.append(krd)

        return key_rate_tenors, np.array(krds, dtype=float)

    ###########################################################################

    def key_rate_durations(
        self,
        settle_dt: Date,
        ytm: float,
        key_rate_tenors: list = None,
        shift: float = None,
        rates: list = None,
    ):
        """
        Calculates the key rate durations for a bond.

        Parameters
        ----------
        settle_dt : FinancePy Date object
            The settlement date.
        ytm : float
            The yield to maturity.
        key_rate_tenors : list of float, optional
            The tenors of the key rates. If None, defaults to
            [0.5, 1, 2, 3, 5, 7, 10, 20, 30].
        shift : float, optional
            The shift used to calculate the key rate durations.
            If None, defaults to 0.0001.
        rates : list of float, optional
            Corresponding yield curve data aligned with key_rate_tenors.
            If None, a flat yield curve at `ytm` is used.

        Returns
        -------
        tuple of (numpy array of float, numpy array of float)
            A tuple containing the key rate tenors and the key rate durations.
        """
        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        # check if key_rate_tenors is None
        # if it is None, create an array of key rates from 0.5 to 30 years

        if key_rate_tenors is None:
            key_rate_tenors = np.array([0.5, 1, 2, 3, 5, 7, 10, 20, 30], dtype=float)
        else:
            key_rate_tenors = np.array(key_rate_tenors, dtype=float)

        # set the shift to a small value if not give
        if shift is None:
            shift = 0.0001

        if shift <= 0.0:
            raise FinError("Shift must be positive")

        # Base par-rate vector
        if rates is None:
            base_rates = np.ones(len(key_rate_tenors), dtype=float) * ytm
        else:
            base_rates = np.array(rates, dtype=float)
            if len(base_rates) != len(key_rate_tenors):
                raise FinError("rates and key_rate_tenors must have the same length")

        lin_zero_interp = InterpTypes.LINEAR_ZERO_RATES
        us_street = YTMCalcType.US_STREET

        def build_par_curve(rate_vec):
            """Construct BondZeroCurve from par bonds implied by rate_vec."""
            par_bonds = []
            clean_prices = []

            for tenor, cpn in zip(key_rate_tenors, rate_vec):
                mat_dt = settle_dt.add_years(float(tenor))
                par_bond = Bond(
                    settle_dt,
                    mat_dt,
                    cpn,
                    self.freq_type,
                    self.accrual_dc_type,
                )
                par_bonds.append(par_bond)

            for par_bond, par_rate in zip(par_bonds, rate_vec):
                clean_price = par_bond.clean_price_from_ytm(
                    settle_dt, par_rate, us_street
                )
                clean_prices.append(clean_price)

            return BondBootstrapDiscountCurve(
                settle_dt, par_bonds, clean_prices, lin_zero_interp
            )

        # Base curve and base price: compute once
        par_crv = build_par_curve(base_rates)
        p_zero = self.dirty_price_from_discount_curve(settle_dt, par_crv)

        key_rate_durations = []

        for ind in range(len(key_rate_tenors)):

            # Create independent shocked copies
            rates_up = base_rates.copy()
            rates_down = base_rates.copy()

            rates_up[ind] += shift
            rates_down[ind] -= shift

            # Rebuild shocked curves
            par_crv_up = build_par_curve(rates_up)
            par_crv_down = build_par_curve(rates_down)

            # Reprice target bond
            p_up = self.dirty_price_from_discount_curve(settle_dt, par_crv_up)
            p_down = self.dirty_price_from_discount_curve(settle_dt, par_crv_down)

            # Central-difference KRD
            krd = (p_down - p_up) / (2.0 * shift * p_zero)
            key_rate_durations.append(krd)

        return key_rate_tenors, np.array(key_rate_durations, dtype=float)

    ###########################################################################

    def convexity_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. This is defined
        as 1/P_D partial 2 P_D / partial dy 2"""
        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1  # / self.par
        return conv

    ###########################################################################

    def clean_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input."""
        validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        dp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        accrued_int = self.accrued_interest(settle_dt, self.par)
        cp = dp - accrued_int
        return cp

    ###########################################################################

    def clean_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        """Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date."""

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        self.accrued_interest(settle_dt, self.par)

        dirty_price = self.dirty_price_from_discount_curve(settle_dt, discount_curve)

        clean_price = dirty_price - self.accrued_int
        return clean_price

    ###########################################################################

    def dirty_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        # 1. Validation checks
        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")
        if settle_dt < discount_curve.value_dt:
            raise FinError("Bond settles before Discount curve date")
        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles after it matures.")

        # 2. Ex-dividend logic
        self._calc_pcd_ncd(settle_dt)
        cal = Calendar(self.cal_type)
        self.ex_div_dt = cal.add_business_days(self._ncd, -self.ex_div_days)

        px_pv = 0.0
        n_flows = len(self.cpn_dts)
        cpn_amount = self.cpn / self.freq

        # 3. Iterate through all potential flows
        for i in range(n_flows):
            cpn_dt = self.cpn_dts[i]
            pmt_dt = self.payment_dts[i]

            # Only consider flows where the coupon date is in the future
            if cpn_dt > settle_dt:

                df = discount_curve.df(pmt_dt)

                # Check for ex-dividend status on the next immediate coupon
                # This is only relevant for the first future coupon encountered
                is_first_future_cpn = cpn_dt == self._ncd

                if is_first_future_cpn and settle_dt > self.ex_div_dt:
                    # Buyer does not receive this coupon
                    flow = 0.0
                else:
                    flow = cpn_amount

                # Add principal if it's the maturity date
                if i == n_flows - 1:
                    flow += 1.0  # Normalized to par=1.0 for now

                px_pv += flow * df

        # 4. Re-base to settlement date
        df_settle = discount_curve.df(settle_dt)
        dirty_price = (px_pv / df_settle) * self.par
        return dirty_price

    ###########################################################################

    def yield_to_maturity(
        self,
        settle_dt: Date,
        clean_price: int | float | list | np.ndarray,
        convention: YTMCalcType = YTMCalcType.US_TREASURY,
    ):
        """Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver."""

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        clean_prices = vectorise_price(clean_price)

        self.accrued_interest(settle_dt, 1.0)

        accrued_amount = self.accrued_int * self.par
        dirty_prices = clean_prices + accrued_amount
        ytms = []

        for dirty_price in dirty_prices:
            argtuple = (self, settle_dt, dirty_price, convention)

            try:
                ytm = optimize.brentq(
                    _f,
                    a=-0.2,  # REDUCED THESE TO -20%
                    b=2.0,  # To +200%
                    args=argtuple,
                    xtol=1e-8,
                    maxiter=50,
                )
            except RuntimeError:
                print(f"Warning: YTM calculation did not converge for price {dirty_price}")
                ytm = np.nan

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        else:
            return np.array(ytms)

    ###########################################################################

    def current_yield(
        self,
        settle_dt: Date,
        clean_price: float | list | np.ndarray
    ):
        """Calculate the bond's simple yield."""

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        self.accrued_interest(settle_dt, 1.0)
        accrued_amount = self.accrued_int * self.par

        clean_prices = vectorise_price(clean_price)
        dirty_prices = clean_prices + accrued_amount
        simple_ys = []

        for dirty_price in dirty_prices:
            simple_y = self.cpn * self.par / dirty_price
            simple_ys.append(simple_y)

        if len(simple_ys) == 1:
            return simple_ys[0]
        else:
            return np.array(simple_ys)

    ###########################################################################

    def _calc_pcd_ncd(self, settle_dt: Date):

        num_flows = len(self.cpn_dts)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dts.")

        for i_flow in range(1, num_flows):
            # coupons paid on a settlement date are paid to the seller
            if self.cpn_dts[i_flow] > settle_dt:
                self._pcd = self.cpn_dts[i_flow - 1]
                self._ncd = self.cpn_dts[i_flow]
                break

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float = 100.0):
        """Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. If the
        bond trades with ex-coupon dates then you need to use the number of
        days before the coupon date the ex-coupon date is. You can specify the
        calendar to be used in the bond constructor - NONE means only calendar
        days, WEEKEND is only weekends, or you can specify a country calendar
        for business days."""

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        self._calc_pcd_ncd(settle_dt)

        dc = DayCount(self.accrual_dc_type)
        cal = Calendar(self.cal_type)

        # Calculation of the ex-dividend date
        self.ex_div_dt = cal.add_business_days(self._ncd, -1 * self.ex_div_days)

        acc_factor, num, _ = dc.year_frac(
            self._pcd, settle_dt, self._ncd, self.freq_type
        )

        self.alpha = 1.0 - acc_factor * self.freq

        if settle_dt > self.ex_div_dt:
            self.accrued_int = acc_factor - 1.0 / self.freq
        else:
            self.accrued_int = acc_factor

        self.accrued_int *= self.cpn * face
        self.accrued_days = num

        return self.accrued_int

    ###########################################################################

    def asset_swap_spread(
        self,
        settle_dt: Date,
        clean_price: float,
        discount_curve: DiscountCurve,
        swap_float_day_count_convention_type=DayCountTypes.ACT_360,
        swap_float_frequency_type=FrequencyTypes.SEMI_ANNUAL,
        swap_float_calendar_type=CalendarTypes.WEEKEND,
        swap_float_bus_day_adjust_rule_type=BusDayAdjustTypes.FOLLOWING,
        swap_float_date_gen_rule_type=DateGenRuleTypes.BACKWARD,
    ):
        """Calculate the par asset swap spread of the bond. The discount curve
        is an Ibor curve that is passed in. This function is vectorised with
        respect to the clean price."""

        clean_price = np.array(clean_price)
        self.accrued_interest(settle_dt, 1.0)
        accrued_amount = self.accrued_int * self.par
        bond_price = clean_price + accrued_amount

        # Calculate the price of the bond discounted on the Ibor curve
        pv_ibor = 0.0
        prev_dt = self._pcd

        for dt in self.payment_dts[1:]:
            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                pv_ibor += df * self.cpn / self.freq

        pv_ibor += discount_curve.df(self.payment_dts[-1])

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prev_dt = self._pcd
        schedule = Schedule(
            settle_dt,
            self.maturity_dt,
            swap_float_frequency_type,
            swap_float_calendar_type,
            swap_float_bus_day_adjust_rule_type,
            swap_float_date_gen_rule_type,
        )

        day_count = DayCount(swap_float_day_count_convention_type)

        prev_dt = self._pcd
        pv01 = 0.0
        for dt in schedule.adjusted_dts[1:]:
            df = discount_curve.df(dt)
            year_frac = day_count.year_frac(prev_dt, dt)[0]
            pv01 = pv01 + year_frac * df
            prev_dt = dt

        asw = (pv_ibor - bond_price / self.par) / pv01
        return asw

    ###########################################################################

    def z_spread(
        self,
        settle_dt: Date,
        clean_price: float,
        discount_curve: DiscountCurve,
    ):
        """Calculate the z-spread of the bond. The discount curve
        is a Ibor curve that is passed in."""

        self.accrued_int = self.accrued_interest(settle_dt, 1.0)
        accrued_amount = self.accrued_int * self.par
        bond_price = clean_price + accrued_amount

        def _bond_price_diff_from_z_spread(z_spr_try):
            flat_curve = PWFONFDiscountCurve.flat_curve(settle_dt, z_spr_try)
            bumped_curve = CompositeDiscountCurve([discount_curve, flat_curve])
            curve_bond_price = self.dirty_price_from_discount_curve(
                settle_dt, bumped_curve
            )
            return curve_bond_price - bond_price

        z_spread = optimize.newton(
            _bond_price_diff_from_z_spread,
            x0=0.0,  # guess initial value of 0%
            fprime=None,
            tol=1e-8,
            maxiter=50,
            fprime2=None,
        )

        return z_spread

    # ###########################################################################

    # def dirty_price_from_oas(
    #     self, settle_dt: Date, discount_curve: DiscountCurve, oas: float
    # ):
    #     """Calculate the dirty price of the bond from its OAS given the bond
    #     settlement date, a discount curve and the oas as a number."""

    #     self.accrued_interest(settle_dt, 1.0)
    #     f = self.freq
    #     c = self.cpn
    #     df_adjusted = 1.0

    #     pv = 0.0
    #     for dt in self.payment_dts[1:]:

    #         # coupons paid on a settlement date are paid to the seller
    #         if dt > settle_dt:
    #             t = (dt - settle_dt) / G_DAYS_IN_YEAR
    #             t = np.maximum(t, G_SMALL)
    #             df = discount_curve.df(dt)

    #             # determine the Ibor implied zero rate
    #             r = f * (np.power(df, -1.0 / t / f) - 1.0)

    #             # determine the OAS adjusted zero rate
    #             df_adjusted = np.power(1.0 + (r + oas) / f, -t * f)
    #             pv = pv + (c / f) * df_adjusted

    #     pv = pv + df_adjusted
    #     pv *= self.par
    #     return pv

    ###########################################################################

    def dirty_price_from_oas(
        self, settle_dt: Date, discount_curve: DiscountCurve, oas: float
    ):
        """Calculate the price of the bond by adding a spread (OAS)
        to the curve."""

        self.accrued_interest(settle_dt, 1.0)
        f = self.freq
        cpn_flow = self.cpn / f
        pv = 0.0
        time_dc_type = discount_curve.time_dc_type

        # We need the settle DF to re-base the price to the settlement date
        df_settle = discount_curve.df(settle_dt)

        n_flows = len(self.payment_dts)
        for i in range(n_flows):

            pmt_dt = self.payment_dts[i]
            cpn_dt = self.cpn_dts[i]  # Use coupon date for ownership logic

            if cpn_dt > settle_dt:

                t = times_from_dates(settle_dt, pmt_dt, time_dc_type)
                t = np.maximum(t, G_SMALL)

                # Get base discount factor from curve
                df_base = discount_curve.df(pmt_dt)

                # Apply OAS (Discrete compounding matching bond frequency)
                df_oas = np.power(1.0 + oas / f, -t * f)
                df_adj = (df_base / df_settle) * df_oas

                # Ownership logic
                if cpn_dt == self._ncd and settle_dt > self.ex_div_dt:
                    flow = 0.0
                else:
                    flow = cpn_flow

                # Add principal at maturity
                if i == n_flows - 1:
                    flow += 1.0

                pv += flow * df_adj

        pv = pv * self.par
        return pv

    ###########################################################################

    def option_adjusted_spread(
        self,
        settle_dt: Date,
        clean_price: float,
        discount_curve: DiscountCurve,
    ):
        """Return OAS for bullet bond given settlement date, clean bond price
        and the discount relative to which the spread is to be computed."""

        if isinstance(clean_price, float) or isinstance(clean_price, np.float64):
            clean_prices = np.array([clean_price])
        elif isinstance(clean_price, list) or isinstance(clean_price, np.ndarray):
            clean_prices = np.array(clean_price)
        else:
            raise FinError("Unknown type for clean_price " + str(type(clean_price)))

        self.accrued_interest(settle_dt, 1.0)

        accrued_amount = self.accrued_int * self.par
        dirty_prices = clean_prices + accrued_amount

        oass = []

        for dirty_price in dirty_prices:
            argtuple = (self, settle_dt, dirty_price, discount_curve)
            try:
                oas = optimize.newton(
                    _g,
                    x0=0.01,  # initial value of 1%
                    fprime=None,
                    args=argtuple,
                    tol=1e-8,
                    maxiter=50,
                    fprime2=None,
                )

                oass.append(oas)
            except RuntimeError:
                print(f"Warning: OAS calculation did not converge for price {dirty_price}")
                oass.append(np.nan)

        if len(oass) == 1:
            return oass[0]
        else:
            return np.array(oass)

    ###########################################################################

    def dirty_price_from_survival_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        survival_curve: DiscountCurve,
        recovery_rate: float,
    ):
        """Calculate discounted present value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the
        defaulting principal we discretize the time steps using the coupon
        payment times. A finer discretization may handle the time value with
        more accuracy. I reduce any error by averaging period start and period
        end payment present values."""

        f = self.freq
        cpn_flow = self.cpn / f

        self.accrued_interest(settle_dt, 1.0)

        df_settle = discount_curve.df(settle_dt)
        q_settle = survival_curve.survival_prob(settle_dt)

        prev_df = df_settle
        prev_q = q_settle

        pv_coupons = 0.0
        pv_recovery = 0.0

        n_flows = len(self.payment_dts)
        for i in range(n_flows):
            pmt_dt = self.payment_dts[i]
            cpn_dt = self.cpn_dts[i]

            if cpn_dt > settle_dt:
                df = discount_curve.df(pmt_dt)
                q = survival_curve.survival_prob(pmt_dt)

                # --- Coupon Leg ---
                # Buyer receives coupon only if bond survives AND it's not ex-div
                if not (cpn_dt == self._ncd and settle_dt > self.ex_div_dt):
                    pv_coupons += cpn_flow * df * q

                # --- Recovery Leg (Protection) ---
                # Probability of default in this specific interval [prev_dt, pmt_dt]
                dq = prev_q - q

                # Average the PV of the recovery payment (Start vs End of period)
                pv_rec_start = dq * recovery_rate * prev_df
                pv_rec_end = dq * recovery_rate * df
                pv_recovery += 0.5 * (pv_rec_start + pv_rec_end)

                # --- Principal (Maturity) ---
                if i == n_flows - 1:
                    pv_coupons += 1.0 * df * q  # Pay par if survives to maturity

                prev_df = df
                prev_q = q

        # Price = (PV of all future contingent flows) / DF at settlement
        # We also divide by q_settle because the price is contingent on
        # the bond having already survived to the settlement date.
        total_pv = (pv_coupons + pv_recovery) / (df_settle * q_settle)

        return total_pv * self.par

    ###########################################################################

    def clean_price_from_survival_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        survival_curve: DiscountCurve,
        recovery_rate: float,
    ):
        """Calculate clean price value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default."""

        self.accrued_interest(settle_dt, 1.0)

        dirty_price = self.dirty_price_from_survival_curve(
            settle_dt, discount_curve, survival_curve, recovery_rate
        )

        clean_price = dirty_price - self.accrued_int
        return clean_price

    ###########################################################################

    def calc_ror(
        self,
        begin_dt: Date,
        end_dt: Date,
        begin_ytm: float,
        end_ytm: float,
        convention: YTMCalcType = YTMCalcType.US_STREET,
    ):
        """
        Calculate the rate of total return(capital return and interest) given a
        BUY YTM and a SELL YTM of this bond.
        This function computes the dirty prices at buying and selling, plus the
        coupon payments during the period.
        It returns a tuple which includes a simple rate of return, a compounded
        IRR and the PnL.
        """
        buy_price = self.dirty_price_from_ytm(begin_dt, begin_ytm, convention)
        sell_price = self.dirty_price_from_ytm(end_dt, end_ytm, convention)
        dts_cfs = zip(self.payment_dts, self.flow_amounts)

        # The coupon or par payments on buying date belong to the buyer. The
        # coupon or par payments on selling date are given to the new buyer.
        dts_cfs = [
            (d, c * self.par) for (d, c) in dts_cfs if (d >= begin_dt) and (d < end_dt)
        ]

        dts_cfs.append((begin_dt, -buy_price))
        dts_cfs.append((end_dt, sell_price))
        times_cfs = [((d - begin_dt) / 365, c) for (d, c) in dts_cfs]
        pnl = sum(c for (t, c) in times_cfs)
        simple_return = (pnl / buy_price) * 365 / (end_dt - begin_dt)
        brentq_up_lim = 5
        brentq_dn_lim = -0.9999

        # in case brentq cannot find the irr root
        if simple_return > brentq_up_lim or simple_return < brentq_dn_lim:
            irr = simple_return
        else:
            irr = optimize.brentq(
                npv,
                # f(a) and f(b) must have opposite signs
                a=brentq_dn_lim,
                b=brentq_up_lim,
                xtol=1e-8,
                args=(np.array(times_cfs),),
            )

        return simple_return, irr, pnl

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("MATURITY_DATE", self.maturity_dt)
        s += label_to_string("COUPON (%)", self.cpn * 100.0)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY_COUNT", self.accrual_dc_type)
        s += label_to_string("EX-DIVIDEND DAYS", self.ex_div_days)
        s += label_to_string("CALENDAR TYPE", self.cal_type)
        s += label_to_string("BUS DAYS ADJUST", self.bd_type)
        s += label_to_string("DATE GEN RULE", self.dg_type)
        s += label_to_string("COUPON TYPE", self.cpn_type, "")

        return s

    ###########################################################################

    def print_payments(self, settle_dt: Date, face: float = 100):
        """Print a list of the unadjusted coupon dts and the adjusted
        payment dates used in the analytic calculations for the bond."""

        # We MUST call this to update value of self.alpha and ex-div date
        self.accrued_interest(settle_dt, 1.0)

        n_flows = len(self.cpn_dts)
        dw = len(str(settle_dt))  # date width
        sw = 14  # status
        aw = 14  # amount

        header = "%*s \t %*s \t %*s  %*s\n" % (
            dw,
            "Coupon Date",
            dw,
            "Payment Date",
            sw,
            "Status",
            aw,
            "Amount",
        )
        flow_str = header
        flow = 0.0

        #        flow_str += "%*s\n" % (dw, str(settle_dt))

        status = "SETTLEMENT"

        flow_str += "%*s \t %*s \t %*s \n" % (
            dw,
            str(settle_dt),
            dw,
            "",
            sw,
            status,
        )

        n_start = -1
        for i in range(0, n_flows):

            cpn_dt = self.cpn_dts[i]
            pmt_dt = self.payment_dts[i]

            if cpn_dt > settle_dt:
                n_start = i
                break

        # If no future coupons found, bond is matured
        if n_start == -1:
            print(flow_str + "Bond has matured.")
            return

        if settle_dt >= self.ex_div_dt:

            cpn_dt = self.cpn_dts[n_start]
            pmt_dt = self.payment_dts[n_start]
            flow_str += "%*s \t  %*s \t %*s \t %*.2f\n" % (
                dw,
                str(cpn_dt),
                dw,
                str(pmt_dt),
                sw,
                "EX-DIV",
                aw,
                0.0,
            )
            n_start += 1

        for i in range(n_start, n_flows):
            cpn_dt = self.cpn_dts[i]
            pmt_dt = self.payment_dts[i]
            flow = self.flow_amounts[i] * face

            status = "UNCHANGED"

            if cpn_dt != pmt_dt:
                status = "HOLIDAY ROLL"

            if cpn_dt == self.maturity_dt:
                status = "MATURITY"

            # Check if this is the final payment (Maturity)
            if i == n_flows - 1:
                flow += face  # Add principal to the last payment

            flow_str += "%*s \t %*s \t %*s \t  %*.2f \n" % (
                dw,
                str(cpn_dt),
                dw,
                str(pmt_dt),
                sw,
                status,
                aw,
                flow,
            )

        print(flow_str)
        return

    ###########################################################################

    def npv(self, settle_dt: Date, discount_curve: DiscountCurve, face: float = 100):
        """Calculate full NPV of bond."""

        print("Under construction")

        # We MUST call this to update value of self.alpha and ex-div date
        self.accrued_interest(settle_dt, 1.0)

        n_flows = len(self.cpn_dts)
        dw = len(str(settle_dt))  # date width
        sw = 14  # status
        aw = 14  # amount

        header = "%*s \t %*s \t %*s  %*s\n" % (
            dw,
            "Coupon Date",
            dw,
            "Payment Date",
            sw,
            "Status",
            aw,
            "Amount",
        )
        flow_str = header
        flow = 0.0

        #        flow_str += "%*s\n" % (dw, str(settle_dt))

        status = "SETTLEMENT"

        flow_str += "%*s \t %*s \t %*s \n" % (
            dw,
            str(settle_dt),
            dw,
            "",
            sw,
            status,
        )

        n_start = -1
        for i in range(0, n_flows):

            cpn_dt = self.cpn_dts[i]
            pmt_dt = self.payment_dts[i]

            if cpn_dt > settle_dt:
                n_start = i
                break

        # If no future coupons found, bond is matured
        if n_start == -1:
            print(flow_str + "Bond has matured.")
            return

        if settle_dt >= self.ex_div_dt:

            cpn_dt = self.cpn_dts[n_start]
            pmt_dt = self.payment_dts[n_start]
            flow_str += "%*s \t  %*s \t %*s \t %*.2f\n" % (
                dw,
                str(cpn_dt),
                dw,
                str(pmt_dt),
                sw,
                "EX-DIV",
                aw,
                0.0,
            )
            n_start += 1

        for i in range(n_start, n_flows):
            cpn_dt = self.cpn_dts[i]
            pmt_dt = self.payment_dts[i]
            flow = self.flow_amounts[i] * face

            status = "UNCHANGED"

            if cpn_dt != pmt_dt:
                status = "HOLIDAY ROLL"

            if cpn_dt == self.maturity_dt:
                status = "MATURITY"

            # Check if this is the final payment (Maturity)
            if i == n_flows - 1:
                flow += face  # Add principal to the last payment

            flow_str += "%*s \t %*s \t %*s \t  %*.2f \n" % (
                dw,
                str(cpn_dt),
                dw,
                str(pmt_dt),
                sw,
                status,
                aw,
                flow,
            )

        print(flow_str)
        return

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dts used in
        analytic calculations for the bond."""
        print(self)


########################################################################################
