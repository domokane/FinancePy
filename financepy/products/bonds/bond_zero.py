import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
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
from ...utils.frequency import FrequencyTypes
from ...products.bonds.bond import YTMCalcType

###############################################################################
# TO DO - THIS CLASS NEEDS TO INHERIT FROM BOND CLASS
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


class BondZero:
    """A zero cpn bond is a bond which doesn't pay any periodic payments.
    Instead, it is issued at a discount. The entire face value of the bond is
    paid out at maturity. It is issued as a deep discount bond.

    There is a special convention for accrued interest in which

        Accrued_interest = (par - issue price) * D

    where D = (settle_dt - issue_dt)/(maturity_dt - issue_dt).
    """

    def __init__(
        self,
        issue_dt: Date,
        maturity_dt: Date,
        issue_price: float,  # Issue price usually discounted
    ):
        """Create BondZero object by providing the issue date, maturity Date,
        face amount and issue price."""

        check_argument_types(self.__init__, locals())

        if issue_dt >= maturity_dt:
            raise FinError("Issue Date must preceded maturity date.")

        self.issue_dt = issue_dt
        self.maturity_dt = maturity_dt
        self.dc_type = DayCountTypes.ZERO
        self.issue_price = issue_price  # Price of issue, usually discounted
        self.par = 100.0  # This is how price is quoted and amount at maturity
        self.freq_type = FrequencyTypes.ZERO
        self.cpn_dts = [issue_dt, maturity_dt]
        self.payment_dts = [issue_dt, maturity_dt]
        self.flow_amounts = [0.0, 0.0]  # cpn payments are zero
        self.cal_type = CalendarTypes.WEEKEND
        self.ex_div_days = 0

        self.accrued_int = None
        self.accrued_days = None
        self.alpha = None
        self.pcd = None
        self.ncd = None

    ###########################################################################

    def dirty_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the full price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM."""

        if convention != YTMCalcType.ZERO:
            raise FinError("Need to use YTMCalcType.ZERO for zero cpn bond")

        self.accrued_interest(settle_dt, 1.0)

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        # n is the number of flows after the next cpn
        n = 0
        for dt in self.cpn_dts:
            if dt > settle_dt:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No cpns left")
        # A zero cpn bond has a price equal to the discounted principal
        # assuming an annualised rate raised to the power of years

        dc = DayCount(self.dc_type)
        (acc_factor, _, _) = dc.year_frac(
            settle_dt, self.maturity_dt, self.maturity_dt, FrequencyTypes.ZERO
        )
        if acc_factor <= 1:
            pv = self.par / (1.0 + ytm * acc_factor)
        else:
            pv = self.par / (1.0 + ytm) ** acc_factor

        return pv

    ###########################################################################

    def principal(
        self,
        settle_dt: Date,
        ytm: float,
        face: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
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
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg."""

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        d = -(p2 - p0) / dy / 2.0
        return d

    ###########################################################################

    def macauley_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity."""

        dd = self.dollar_duration(settle_dt, ytm, convention)
        fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd * (1.0 + ytm) / fp
        return md

    ###########################################################################

    def modified_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the modified duration of the bondon a settlement date
        given its yield to maturity."""

        dd = self.dollar_duration(settle_dt, ytm, convention)
        fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd / fp
        return md

    ###########################################################################

    def convexity_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input."""

        dy = 0.0001
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
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input."""

        dirty_price = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        accrued = self.accrued_interest(settle_dt, self.par)
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def clean_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        """Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date."""

        dirty_price = self.dirty_price_from_discount_curve(
            settle_dt, discount_curve
        )

        accrued = self.accrued_interest(settle_dt, self.par)
        clean_price = dirty_price - accrued
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

        px = 0.0
        df = 1.0
        df_settle = discount_curve.df(
            settle_dt,
        )

        for dt in self.cpn_dts[1:]:

            # cpns paid on the settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                flow = 0
                pv = flow * df
                px += pv

        px += df * self.par
        px = px / df_settle

        return px * self.par

    ###########################################################################

    def current_yield(self, clean_price):
        """
        Calculate the current yield of the bond which is the
        cpn divided by the clean price (not the full price).
        The cpn of a zero cpn bond is defined as:
        (par - issue_price) / tenor
        """
        dc = DayCount(self.dc_type)
        tenor, _, _ = dc.year_frac(
            self.issue_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )
        virtual_cpn = (self.par - self.issue_price) / tenor
        y = virtual_cpn / clean_price
        return y

    ###########################################################################

    def yield_to_maturity(
        self,
        settle_dt: Date,
        clean_price: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
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

        accrued_amount = self.accrued_interest(settle_dt, self.par)
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

    def accrued_interest(self, settle_dt: Date, face: float):
        """Calculate the amount of cpn that has accrued between the
        previous cpn date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous cpn payment date and settlement date. If the
        bond trades with ex-cpn dates then you need to supply the number of
        days before the cpn date the ex-cpn date is. You can specify the
        calendar to be used - NONE means only calendar days, WEEKEND is only
        weekends or you can specify a country calendar for business days."""

        num_flows = len(self.cpn_dts)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond Zero settlement after maturity date")

        for i_flow in range(1, num_flows):

            # cpns paid on the settlement date are paid to the seller
            if self.cpn_dts[i_flow] > settle_dt:
                self.pcd = self.cpn_dts[i_flow - 1]
                self.ncd = self.cpn_dts[i_flow]
                break

        dc = DayCount(self.dc_type)
        cal = Calendar(self.cal_type)

        ex_dividend_dt = cal.add_business_days(self.ncd, -self.ex_div_days)

        (acc_factor, num, _) = dc.year_frac(
            self.pcd, settle_dt, self.ncd, FrequencyTypes.ZERO
        )

        if settle_dt > ex_dividend_dt:
            acc_factor = acc_factor - 1.0

        self.alpha = 1.0 - acc_factor

        num = settle_dt - self.issue_dt
        den = self.maturity_dt - self.issue_dt

        f = num / den
        g = ((self.par - self.issue_price)) / self.par

        self.accrued_int = f * g * face
        self.accrued_days = num

        return self.accrued_int

    ###########################################################################

    def asset_swap_spread(
        self,
        settle_dt: Date,
        clean_price: float,
        discount_curve: DiscountCurve,
        swapFloatDayCountConventionType=DayCountTypes.ACT_360,
        swap_float_freq_type=FrequencyTypes.SEMI_ANNUAL,
        swap_float_cal_type=CalendarTypes.WEEKEND,
        swap_float_bus_day_adjust_rule_type=BusDayAdjustTypes.FOLLOWING,
        swapFloatDateGenRuleType=DateGenRuleTypes.BACKWARD,
    ):
        """Calculate the par asset swap spread of the bond. The discount curve
        is a Ibor curve that is passed in. This function is vectorised with
        respect to the clean price."""

        clean_price = np.array(clean_price)
        self.accrued_interest(settle_dt, 1.0)
        accrued_amount = self.accrued_int * self.par
        bond_price = clean_price + accrued_amount
        # Calculate the price of the bond discounted on the Ibor curve
        pv_ibor = 0.0
        prev_dt = self.pcd

        for dt in self.cpn_dts[1:]:

            # cpns paid on the settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                # pv_ibor += df * self.cpn / self.freq

        pv_ibor += df * self.par

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the cpn starts accruing on the settlement date
        prev_dt = self.pcd
        schedule = Schedule(
            settle_dt,
            self.maturity_dt,
            swap_float_freq_type,
            swap_float_cal_type,
            swap_float_bus_day_adjust_rule_type,
            swapFloatDateGenRuleType,
        )

        day_count = DayCount(swapFloatDayCountConventionType)

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

    def dirty_price_from_oas(
        self, settle_dt: Date, discount_curve: DiscountCurve, oas: float
    ):
        """Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number."""

        self.accrued_interest(settle_dt, 1.0)

        pv = 0.0
        for dt in self.cpn_dts[1:]:

            # cpns paid on the settlement date are paid to the seller
            if dt > settle_dt:
                t = (dt - settle_dt) / g_days_in_year

                t = np.maximum(t, g_small)

                df = discount_curve.df(dt)
                # determine the Ibor implied zero rate
                r = np.power(df, -1.0 / t) - 1.0
                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas), -t)

        pv = pv + df_adjusted * self.par
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
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""
        flow_str = ""
        flow_str += "%12s %12.2f \n" % (self.cpn_dts[-1], face)

        return flow_str

    ###########################################################################

    def print_bond_payments(self, settle_dt: Date, face: float = 100.0):
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""

        print(self.bond_payments(settle_dt, face))

    ###########################################################################

    def dirty_price_from_survival_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        survival_curve: DiscountCurve,
        recovery_rate: float,
    ):
        """Calculate discounted present value of flows assuming default model.
        The survival curve treats the cpns as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the
        defaulting principal we discretize the time steps using the cpn
        payment times. A finer discretization may handle the time value with
        more accuracy. I reduce any error by averaging period start and period
        end payment present values."""

        pv = 0.0
        prev_q = 1.0
        prev_df = 1.0

        defaulting_principal_pv_pay_start = 0.0
        defaulting_principal_pv_pay_end = 0.0

        for dt in self.cpn_dts[1:]:

            # cpns paid on the settlement date are paid to the seller
            if dt > settle_dt:
                df = discount_curve.df(dt)
                q = survival_curve.survival_prob(dt)

                # Add PV of cpn conditional on surviving to payment date
                # Any default results in all subsequent cpns being lost
                # with zero recovery

                dq = q - prev_q

                defaulting_principal_pv_pay_start += (
                    -dq * recovery_rate * prev_df
                )
                defaulting_principal_pv_pay_start += -dq * recovery_rate * df

                # Add on PV of principal if default occurs in cpn period
                prev_q = q
                prev_df = df

        pv = pv + 0.50 * defaulting_principal_pv_pay_start
        pv = pv + 0.50 * defaulting_principal_pv_pay_end
        pv = pv + df * q * self.par
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
        The survival curve treats the cpns as zero recovery payments while
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
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """TODO: TIDY THIS UP !!!!!!!!!!!!!!!!!

        Calculates the rate of total return (capital return and interest) given
        a BUY YTM and a SELL YTM of this bond.

        This function computes the full prices at buying and selling, plus the
        cpn payments during the period.

        It returns a tuple which includes a simple rate of return, a compounded
        IRR and the PnL.
        """

        buy_price = self.dirty_price_from_ytm(begin_dt, begin_ytm, convention)
        sell_price = self.dirty_price_from_ytm(end_dt, end_ytm, convention)

        dates_cfs = zip(self.cpn_dts, self.flow_amounts)

        # The cpn or par payments on buying date belong to the buyer.
        # The cpn or par payments on selling date are given to the new buyer
        dates_cfs = [
            (d, c * self.par)
            for (d, c) in dates_cfs
            if (d >= begin_dt) and (d < end_dt)
        ]

        dates_cfs.append((begin_dt, -buy_price))

        dates_cfs.append((end_dt, sell_price))
        times_cfs = [((d - begin_dt) / 365, c) for (d, c) in dates_cfs]

        pnl = sum(c for (t, c) in times_cfs)
        simple_rtn = (pnl / buy_price) * 365 / (end_dt - begin_dt)
        brentq_up_bound = 5
        brentq_down_bound = -0.9999

        # in case brentq cannot find the irr root

        if simple_rtn > brentq_up_bound or simple_rtn < brentq_down_bound:

            irr = simple_rtn

        else:

            irr = optimize.brentq(
                npv,
                # f(a) and f(b) must have opposite signs
                a=brentq_down_bound,
                b=brentq_up_bound,
                xtol=1e-8,
                args=(np.array(times_cfs),),
            )

        return simple_rtn, irr, pnl

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("COUPON (%)", 0)
        s += label_to_string("ISSUE PRICE", self.issue_price)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
