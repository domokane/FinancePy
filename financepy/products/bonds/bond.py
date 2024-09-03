###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# TODO: - ROUNDING CONVENTIONS FOR ACCRUED
# TODO: - CHECK OAS CALCULATION
# TODO:  - Check how first coupon on floating leg is sized on asset swaps. """
###############################################################################

# https://www.dmo.gov.uk/media/15004/convention_changes.pdf

###############################################################################
# Conventions:
#  GILTS - SEMI ANNUAL ACT/ACT
#  US TREASURIES
###############################################################################

###############################################################################
# NOTE THAT I ASSUME THAT IF YOU SETTLE A SWAP ON A COUPON PAYMENT DATE YOU
# GET THE COUPON AND THE ACCRUED INTEREST EQUALS THE COUPON.
###############################################################################

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.global_vars import g_days_in_year, g_small
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.schedule import Schedule
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.math import npv
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.interpolator import InterpTypes
from .bond_zero_curve import BondZeroCurve
from ...market.curves.discount_curve_pwf_onf import DiscountCurvePWFONF
from ...market.curves.composite_discount_curve import CompositeDiscountCurve


# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

###############################################################################


from enum import Enum


class YTMCalcType(Enum):
    ZERO = 0
    UK_DMO = 1
    US_STREET = 2
    US_TREASURY = 3
    CFETS = 4  # China Foreign Exchange Trade System


###############################################################################


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
        dc_type: DayCountTypes,
        ex_div_days: int = 0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type=BusDayAdjustTypes.FOLLOWING,
        dg_type=DateGenRuleTypes.BACKWARD,
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

        if dc_type == dc_type.ZERO:
            raise FinError("Zero coupon bonds must use BondZero class.")

        self.cpn = coupon
        self.freq_type = freq_type

        self.dc_type = dc_type
        self.freq = annual_frequency(freq_type)

        if ex_div_days > 90:
            raise FinError(
                "Ex dividend days cannot be more than 90" + str(ex_div_days)
            )

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
        self._calculate_flows()

        self.pcd = None
        self.ncd = None

        # Private
        self.accrued_int = None
        self.accrued_days = 0.0
        self.alpha = 0.0

    ###########################################################################

    def _calculate_cpn_dts(self):
        """Determine the bond coupon dts. Note that for analytical
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

    ###########################################################################

    def _calculate_payment_dts(self):
        """For the actual payment dts, they are adjusted,
        and so we then use the calendar payment dts. Although payments are
        calculated as though coupon periods are the same length, payments that
        fall on a Saturday or Sunday can only be made on the next business day
        """

        # dts are adjusted forward to the next business day
        bus_day_adj_type = BusDayAdjustTypes.FOLLOWING
        calendar = Calendar(self.cal_type)

        self._calculate_cpn_dts()

        self.payment_dts = []

        # Expect at least an issue date and a maturity date - if not - problem
        if len(self.cpn_dts) < 2:
            raise FinError("Cannot calculate payment dts with one payment")

        # I do not adjust the first date as it is the issue date
        self.payment_dts.append(self.cpn_dts[0])

        for cpn_dt in self.cpn_dts[1:]:
            pmt_dt = calendar.adjust(cpn_dt, bus_day_adj_type)

            self.payment_dts.append(pmt_dt)

    ###########################################################################

    def _calculate_flows(self):
        """Determine the bond cash flow payment amounts without principal.
        There is no adjustment based on the adjusted payment dts."""

        self.flow_amounts = [0.0]

        for _ in self.cpn_dts[1:]:
            cpn = self.cpn / self.freq
            self.flow_amounts.append(cpn)

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

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settlement is after maturity date")

        if convention not in YTMCalcType:
            raise FinError("Yield convention unknown." + str(convention))

        # TODO check that no unnecessary calculations are being done
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

        # n is the number of flows after the next coupon
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
        #                print(term1, term2, term3, term4, v, self.alpha, dp)
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
        else:
            raise FinError("Unknown yield convention")

        return dp * self.par

    ###########################################################################

    def principal(
        self, settle_dt: Date, ytm: float, face: float, convention: YTMCalcType
    ):
        """Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates."""

        dirty_price = self.dirty_price_from_ytm(settle_dt, ytm, convention)

        principal = dirty_price * face / self.par
        principal = principal - self.accrued_int
        return principal

    ###########################################################################

    def dollar_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg."""

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

    ###########################################################################

    def macauley_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity."""

        dd = self.dollar_duration(settle_dt, ytm, convention)
        fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd * (1.0 + ytm / self.freq) / fp
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

        dd = self.dollar_duration(settle_dt, ytm, convention)
        fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd / fp
        return md

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
        bond : FinancePy Bond object

        settle_dt : FinancePy Date object
            The settlement date.
        ytm : float
            The yield to maturity.
        key_rate_tenors : list of float, optional
            The tenors of the key rates, default is None which will generate
            the tenors from 0.25 to 30 years.
        shift : float, optional
            The shift used to calculate the key rate durations, default is None
            which will set the shift to 0.0001.
        rates: list of float, optional
            Corresponding yield curve data in line with key_rate_tenors
            If None, flat yield curve is used

        Returns
        -------
        tuple of (numpy array of float, numpy array of float)
            A tuple containing the key rate tenors and the key rate durations.
        """

        # check if key_rate_tenors is None
        # if it is None, create an array of key rates from 0.5 to 30 years

        if key_rate_tenors is None:
            key_rate_tenors = np.array([0.5, 1, 2, 3, 5, 7, 10, 20, 30])

        # set the shift to a small value if not give
        if not shift:
            shift = 0.0001

        lin_zero_interp = InterpTypes.LINEAR_ZERO_RATES
        us_street = YTMCalcType.US_STREET

        # initialize an empty list for the key rate durations
        key_rate_durations = []

        # iterate over each key rate (tenor) and calculate key rate duration
        for ind, _ in enumerate(key_rate_tenors):

            # if rates are not given create an array of rates where each
            # rate is equal to the ytm value
            if rates is None:
                rates = np.ones(len(key_rate_tenors)) * ytm

            # Create set of par bonds to be used in BondZeroCurve
            # ytm and coupons are equal
            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):
                mat_dt = settle_dt.add_years(tenor)

                par_bond = Bond(
                    settle_dt, mat_dt, cpn, self.freq_type, self.dc_type
                )

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):
                clean_price = par_bond.clean_price_from_ytm(
                    settle_dt, ytm, us_street
                )
                clean_prices.append(clean_price)

            par_crv = BondZeroCurve(
                settle_dt, par_bonds, clean_prices, lin_zero_interp
            )

            # calculate the dirty price of the bond using the discount curve
            p_zero = self.dirty_price_from_discount_curve(settle_dt, par_crv)

            # shift up by the yield of corresponding par bond
            rates[ind] += shift

            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):
                mat = settle_dt.add_years(tenor)

                par_bond = Bond(
                    settle_dt, mat, cpn, self.freq_type, self.dc_type
                )

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):
                clean_price = par_bond.clean_price_from_ytm(
                    settle_dt, ytm, us_street
                )
                clean_prices.append(clean_price)

            par_crv_up = BondZeroCurve(
                settle_dt, par_bonds, clean_prices, lin_zero_interp
            )

            # calculate the full price of the bond
            # using the discount curve with the key rate shifted up
            p_up = self.dirty_price_from_discount_curve(settle_dt, par_crv_up)

            # create a curve again with the key rate shifted down
            # by twice the shift value.
            rates[ind] -= shift * 2

            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):
                mat = settle_dt.add_years(tenor)

                par_bond = Bond(
                    settle_dt, mat, cpn, self.freq_type, self.dc_type
                )

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):
                clean_price = par_bond.clean_price_from_ytm(
                    settle_dt, ytm, us_street
                )
                clean_prices.append(clean_price)

            par_crv_down = BondZeroCurve(
                settle_dt, par_bonds, clean_prices, lin_zero_interp
            )

            # calculate the full price of the bond using
            p_down = self.dirty_price_from_discount_curve(
                settle_dt, par_crv_down
            )

            # calculate the key rate duration
            # using the formula (P_down - P_up) / (2 * shift * P_zero)
            key_rate_duration = (p_down - p_up) / (2 * shift * p_zero)

            # append the key rate duration to the key_rate_durations list
            key_rate_durations.append(key_rate_duration)

        return key_rate_tenors, np.array(key_rate_durations)

    ###########################################################################

    def convexity_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ):
        """Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input."""

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self.par
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

        self.accrued_interest(settle_dt, self.par)

        dirty_price = self.dirty_price_from_discount_curve(
            settle_dt, discount_curve
        )

        clean_price = dirty_price - self.accrued_int
        return clean_price

    ###########################################################################

    def dirty_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        """Calculate the bond price using a provided discount curve to PV the
        bond's cash flows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        """

        if settle_dt < discount_curve.value_dt:
            raise FinError("Bond settles before Discount curve date")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles after it matures.")

        self._calc_pcd_ncd(settle_dt)

        cal = Calendar(self.cal_type)

        self.ex_div_dt = cal.add_business_days(self.ncd, -self.ex_div_days)

        pay_first_cpn = 1.0
        if settle_dt > self.ex_div_dt:
            pay_first_cpn = 0.0

        px = 0.0
        df = 1.0
        df_settle_dt = discount_curve.df(settle_dt)

        dt = self.cpn_dts[1]
        if dt > settle_dt:
            df = discount_curve.df(dt)
            flow = self.cpn / self.freq
            pv = flow * df
            px += pv * pay_first_cpn

        for dt in self.cpn_dts[2:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                flow = self.cpn / self.freq
                pv = flow * df
                px += pv

        px += df
        px = px / df_settle_dt

        return px * self.par

    ###########################################################################

    def current_yield(self, clean_price):
        """Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)"""

        y = self.cpn * self.par / clean_price
        return y

    ###########################################################################

    def yield_to_maturity(
        self,
        settle_dt: Date,
        clean_price: float,
        convention: YTMCalcType = YTMCalcType.US_TREASURY,
    ):
        """Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver."""

        if isinstance(clean_price, float) or isinstance(
            clean_price, np.float64
        ):
            clean_prices = np.array([clean_price])
        elif isinstance(clean_price, list) or isinstance(
            clean_price, np.ndarray
        ):
            clean_prices = np.array(clean_price)
        else:
            raise FinError(
                "Unknown type for clean_price " + str(type(clean_price))
            )

        self.accrued_interest(settle_dt, 1.0)

        accrued_amount = self.accrued_int * self.par
        dirty_prices = clean_prices + accrued_amount
        ytms = []

        for dirty_price in dirty_prices:
            argtuple = (self, settle_dt, dirty_price, convention)

            ytm = optimize.newton(
                _f,
                x0=0.05,  # guess initial value of 5%
                fprime=None,
                args=argtuple,
                tol=1e-8,
                maxiter=50,
                fprime2=None,
            )

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        else:
            return np.array(ytms)

    ###########################################################################

    def _calc_pcd_ncd(self, settle_dt: Date):

        num_flows = len(self.cpn_dts)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dts.")

        for i_flow in range(1, num_flows):
            # coupons paid on a settlement date are paid to the seller
            if self.cpn_dts[i_flow] > settle_dt:
                self.pcd = self.cpn_dts[i_flow - 1]
                self.ncd = self.cpn_dts[i_flow]
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

        self._calc_pcd_ncd(settle_dt)

        dc = DayCount(self.dc_type)
        cal = Calendar(self.cal_type)

        self.ex_div_dt = cal.add_business_days(self.ncd, -1 * self.ex_div_days)

        (acc_factor, num, _) = dc.year_frac(
            self.pcd, settle_dt, self.ncd, self.freq_type
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
        prev_dt = self.pcd

        for dt in self.cpn_dts[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                pv_ibor += df * self.cpn / self.freq

        pv_ibor += df

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prev_dt = self.pcd
        schedule = Schedule(
            settle_dt,
            self.maturity_dt,
            swap_float_frequency_type,
            swap_float_calendar_type,
            swap_float_bus_day_adjust_rule_type,
            swap_float_date_gen_rule_type,
        )

        day_count = DayCount(swap_float_day_count_convention_type)

        prev_dt = self.pcd
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
        settlement_date: Date,
        clean_price: float,
        discount_curve: DiscountCurve,
    ):
        """Calculate the z-spread of the bond. The discount curve
        is a Ibor curve that is passed in."""

        self.accrued_int = self.accrued_interest(settlement_date, 1.0)
        accrued_amount = self.accrued_int * self.par
        bondPrice = clean_price + accrued_amount

        def _bond_price_diff_from_z_spread(z_spr_try):
            flat_curve = DiscountCurvePWFONF.flat_curve(
                settlement_date, z_spr_try
            )
            bumped_curve = CompositeDiscountCurve([discount_curve, flat_curve])
            curve_bond_price = self.dirty_price_from_discount_curve(
                settlement_date, bumped_curve
            )
            return curve_bond_price - bondPrice

        z_spread = optimize.newton(
            _bond_price_diff_from_z_spread,
            x0=0.0,  # guess initial value of 0%
            fprime=None,
            tol=1e-8,
            maxiter=50,
            fprime2=None,
        )[0]
        return z_spread

    ###########################################################################

    def dirty_price_from_oas(
        self, settle_dt: Date, discount_curve: DiscountCurve, oas: float
    ):
        """Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number."""

        self.accrued_interest(settle_dt, 1.0)
        f = self.freq
        c = self.cpn
        df_adjusted = 1.0

        pv = 0.0
        for dt in self.cpn_dts[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                t = (dt - settle_dt) / g_days_in_year

                t = np.maximum(t, g_small)

                df = discount_curve.df(dt)

                # determine the Ibor implied zero rate
                r = f * (np.power(df, -1.0 / t / f) - 1.0)

                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas) / f, -t * f)
                pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted
        pv *= self.par
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

        if isinstance(clean_price, float) or isinstance(
            clean_price, np.float64
        ):
            clean_prices = np.array([clean_price])
        elif isinstance(clean_price, list) or isinstance(
            clean_price, np.ndarray
        ):
            clean_prices = np.array(clean_price)
        else:
            raise FinError(
                "Unknown type for clean_price " + str(type(clean_price))
            )

        self.accrued_interest(settle_dt, 1.0)

        accrued_amount = self.accrued_int * self.par
        dirty_prices = clean_prices + accrued_amount

        oass = []

        for dirty_price in dirty_prices:
            argtuple = (self, settle_dt, dirty_price, discount_curve)

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

        if len(oass) == 1:
            return oass[0]
        else:
            return np.array(oass)

    ###########################################################################

    def bond_payments(self, settle_dt: Date, face: float):
        """Print a list of the unadjusted coupon payment dts used in
        analytic calculations for the bond."""

        flow = face * self.cpn / self.freq

        flow_str = ""

        for dt in self.cpn_dts[1:-1]:
            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                flow_str += "%12s %12.5f \n" % (dt, flow)

        redemption_amount = face + flow
        flow_str += "%12s %12.5f \n" % (self.cpn_dts[-1], redemption_amount)

        return flow_str

    ###########################################################################

    def print_payments(self, settle_dt: Date, face: float = 100.0):
        """Print a list of the unadjusted coupon payment dts used in
        analytic calculations for the bond."""

        print(self.bond_payments(settle_dt, face))

    def print_coupon_dates(self, settle_dt: Date):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""

        print(self.coupon_dates(settle_dt))

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
        c = self.cpn

        pv = 0.0
        prev_q = 1.0
        prev_df = 1.0

        defaulting_pv_pay_start = 0.0
        defaulting_pv_pay_end = 0.0

        for dt in self.cpn_dts[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                q = survival_curve.survival_prob(dt)

                # Add PV of coupon conditional on surviving to payment date
                # Any default results in all subsequent coupons being lost
                # with zero recovery

                pv = pv + (c / f) * df * q
                dq = q - prev_q

                defaulting_pv_pay_start += -dq * recovery_rate * prev_df
                defaulting_pv_pay_end += -dq * recovery_rate * df

                # Add on PV of principal if default occurs in coupon period
                prev_q = q
                prev_df = df

        pv = pv + 0.50 * defaulting_pv_pay_start
        pv = pv + 0.50 * defaulting_pv_pay_end
        pv = pv + df * q
        pv *= self.par
        return pv

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
        This function computes the full prices at buying and selling, plus the
        coupon payments during the period.
        It returns a tuple which includes a simple rate of return, a compounded
        IRR and the PnL.
        """
        buy_price = self.dirty_price_from_ytm(begin_dt, begin_ytm, convention)
        sell_price = self.dirty_price_from_ytm(end_dt, end_ytm, convention)
        dts_cfs = zip(self.cpn_dts, self.flow_amounts)

        # The coupon or par payments on buying date belong to the buyer. The
        # coupon or par payments on selling date are given to the new buyer.
        dts_cfs = [
            (d, c * self.par)
            for (d, c) in dts_cfs
            if (d >= begin_dt) and (d < end_dt)
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

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("COUPON (%)", self.cpn * 100.0)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        s += label_to_string("EX_DIV DAYS", self.ex_div_days, "")
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dts used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
