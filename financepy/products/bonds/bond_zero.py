from typing import List

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...products.bonds.bond import YTMCalcType


def validate_yield(ytm):
    if ytm < -0.20 or ytm > 2.0:
        raise FinError("YTM must be between -20% and +200%")


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
        issue_price: float,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,  # ONLY NEEDED FOR YIELD CALCULATIONS
    ):
        """Create BondZero object by providing the issue date, maturity Date,
        face amount and issue price."""

        check_argument_types(self.__init__, locals())

        if issue_dt >= maturity_dt:
            raise FinError("Issue Date must preceded maturity date.")

        self.issue_dt = issue_dt
        self.maturity_dt = maturity_dt
        self.issue_price = issue_price
        self.par = 100.0  # This is how price is quoted and amount at maturity
        self.freq_type = FrequencyTypes.ZERO
        self.time_dc_type = time_dc_type

        self.accrued_int = None
        self.accrued_days = None
        self.alpha = None
        self._pcd = None
        self._ncd = None

    ###########################################################################

    def dirty_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        if convention != YTMCalcType.ZERO:
            raise FinError("Need YTMCalcType.ZERO for zero coupon bond")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles on or after maturity")

        validate_yield(ytm)

        # year fraction from settlement to maturity
        dc = DayCount(self.time_dc_type)
        T, _, _ = dc.year_frac(
            settle_dt, self.maturity_dt, self.maturity_dt, FrequencyTypes.ZERO
        )

        ytm = np.asarray(ytm)  # keep vectorisation

        if T <= 1.0:
            # money market convention
            return self.par / (1.0 + ytm * T)
        else:
            # annual compounding
            return self.par / (1.0 + ytm) ** T

    ###########################################################################

    def principal(
        self,
        settle_dt: Date,
        ytm: float,
        face: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Return BBG principal value for a given face amount."""

        validate_yield(ytm)

        dirty_price = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        accrued = self.accrued_interest(settle_dt, face)

        return dirty_price * face / self.par - accrued

    ###########################################################################

    def accretion_yield(self, settle_dt: Date, clean_price):
        """Return average annual accretion yield from settlement."""

        if settle_dt > self.maturity_dt:
            raise FinError("Settlement on or after maturity")

        clean_price = np.asarray(clean_price)

        dc = DayCount(self.time_dc_type)
        T, _, _ = dc.year_frac(
            settle_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )

        if T <= 0.0:
            raise FinError("Invalid time to maturity")

        accrued = self.accrued_interest(settle_dt, self.par)
        dirty_price = clean_price + accrued

        annual_accretion = (self.par - dirty_price) / T
        y = annual_accretion / clean_price

        if y.ndim:
            return y
        else:
            return float(y)

    ###########################################################################

    def current_yield(self, clean_price):
        """
        Calculate the current yield of the bond which is the
        cpn divided by the clean price (not the full price).
        The cpn of a zero cpn bond is defined as:
        (par - issue_price) / tenor
        """
        dc = DayCount(self.dc_type)
        T, _, _ = dc.year_frac(
            self.issue_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )
        virtual_cpn = (self.par - self.issue_price) / T
        y = virtual_cpn / clean_price
        return y

    ###########################################################################

    def yield_to_maturity(
        self,
        settle_dt: Date,
        clean_price: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Return YTM from clean price for a zero coupon bond."""

        if convention != YTMCalcType.ZERO:
            raise FinError("Need YTMCalcType.ZERO for zero coupon bond")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles on or after maturity")

        clean_price = np.asarray(clean_price)

        # convert to dirty price
        accrued = self.accrued_interest(settle_dt, self.par)
        dirty_price = clean_price + accrued

        # time to maturity
        dc = DayCount(self.time_dc_type)
        T, _, _ = dc.year_frac(
            settle_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )

        if T <= 0:
            raise FinError("Invalid time to maturity")

        # analytic inversion
        if T <= 1.0:
            # simple interest convention
            ytm = (self.par / dirty_price - 1.0) / T
        else:
            # annual compounding
            ytm = (self.par / dirty_price) ** (1.0 / T) - 1.0

        return ytm if ytm.ndim else float(ytm)

    ###########################################################################

    def discount_rate(
        self,
        settle_dt: Date,
        clean_price: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Return discount rate for a zero coupon bond."""

        if convention != YTMCalcType.ZERO:
            raise FinError("Need YTMCalcType.ZERO for zero coupon bond")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles on or after maturity")

        clean_price = np.asarray(clean_price)

        # convert to dirty price
        accrued = self.accrued_interest(settle_dt, self.par)
        dirty_price = clean_price + accrued

        # time to maturity
        dc = DayCount(self.time_dc_type)
        T, _, _ = dc.year_frac(
            settle_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )

        if T <= 0:
            raise FinError("Invalid time to maturity")

        # analytic inversion
        if T <= 1.0:
            # simple interest convention
            dr = (1.0 - dirty_price / self.par) / T
        else:
            # annual compounding
            dr = (1.0 - dirty_price / self.par) ** (1.0 / T) - 1.0

        return dr if dr.ndim else float(dr)

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float = 100.0):
        """Return accrued discount on a zero coupon bond.

        Convention:
            accrued = face * (par - issue_price) / par
                      * (settle_dt - issue_dt) / (maturity_dt - issue_dt)
        """

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date before issue date")

        if settle_dt > self.maturity_dt:
            raise FinError("Settlement date on or after maturity date")

        elapsed_days = settle_dt - self.issue_dt
        total_days = self.maturity_dt - self.issue_dt

        accrual_fraction = elapsed_days / total_days
        total_discount = (self.par - self.issue_price) / self.par

        self.accrued_int = face * total_discount * accrual_fraction
        self.accrued_days = elapsed_days

        self._pcd = self.issue_dt
        self._ncd = self.maturity_dt
        self.alpha = 1.0 - accrual_fraction

        return self.accrued_int

    ###########################################################################

    def dollar_duration(
        self,
        settle_dt: (Date, List),
        ytm: (float, List),
        convention: YTMCalcType = YTMCalcType.ZERO,
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
        convention: YTMCalcType = YTMCalcType.ZERO,
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
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Return Macaulay duration of a zero coupon bond."""
        validate_yield(ytm)

        dc = DayCount(self.time_dc_type)
        T, _, _ = dc.year_frac(
            settle_dt,
            self.maturity_dt,
            self.maturity_dt,
            FrequencyTypes.ZERO,
        )

        return T

    ###########################################################################

    def modified_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the modified duration of the bondon a settlement date
        given its yield to maturity."""
        validate_yield(ytm)

        dd = self.dollar_duration(settle_dt, ytm, convention)
        fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        md = dd / fp * 10000
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

        validate_yield(ytm)

        dy = 0.0001
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / (dy * dy * p1)
        return conv

    ###########################################################################

    def clean_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Calculate the clean price of a zero coupon bond from its YTM."""
        validate_yield(ytm)

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
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
    ):
        """Return dirty price per 100 nominal using a discount curve."""

        if settle_dt < discount_curve.value_dt:
            raise FinError("Bond settles before discount curve date")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles on or after maturity")

        df_settle = discount_curve.df(settle_dt)
        df_mat = discount_curve.df(self.maturity_dt)

        return self.par * df_mat / df_settle

    def bond_payments(self, settle_dt: Date, face: float):
        """Return the zero bond redemption payment."""

        if settle_dt > self.maturity_dt:
            return ""

        return "%12s %12.2f \n" % (self.maturity_dt, face)

    ###########################################################################

    def print_payments(self, settle_dt: Date, face: float = 100.0):
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
        """Return risky dirty price per 100 nominal for a zero coupon bond."""

        if settle_dt < discount_curve.value_dt:
            raise FinError("Bond settles before discount curve date")

        if settle_dt > self.maturity_dt:
            raise FinError("Bond settles on or after maturity")

        df_settle = discount_curve.df(settle_dt)
        df_mat = discount_curve.df(self.maturity_dt)

        q_settle = survival_curve.survival_prob(settle_dt)

        if q_settle <= 0.0:
            raise FinError(
                "Survival probability at settlement must be positive"
            )

        q_mat = survival_curve.survival_prob(self.maturity_dt)

        df_forward = df_mat / df_settle
        q_forward = q_mat / q_settle

        if q_forward < 0.0 or q_forward > 1.0:
            raise FinError("Invalid conditional survival probability")

        survival_pv = self.par * df_forward * q_forward

        default_pv = self.par * recovery_rate * df_forward * (1.0 - q_forward)

        return survival_pv + default_pv

    ###########################################################################

    def clean_price_from_survival_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        survival_curve: DiscountCurve,
        recovery_rate: float,
    ):
        """Return risky clean price per 100 nominal for a zero coupon bond."""

        dirty_price = self.dirty_price_from_survival_curve(
            settle_dt,
            discount_curve,
            survival_curve,
            recovery_rate,
        )

        accrued = self.accrued_interest(settle_dt, self.par)

        return dirty_price - accrued

    ###########################################################################

    def calc_ror(
        self,
        begin_dt: Date,
        end_dt: Date,
        begin_ytm: float,
        end_ytm: float,
        convention: YTMCalcType = YTMCalcType.ZERO,
    ):
        """Return annualised simple return, IRR, and PnL for a zero coupon bond."""

        if begin_dt >= end_dt:
            raise FinError("Begin date must be before end date")

        if begin_dt > self.maturity_dt:
            raise FinError("Begin date on or after maturity")

        buy_price = self.dirty_price_from_ytm(begin_dt, begin_ytm, convention)

        if end_dt >= self.maturity_dt:
            sell_price = self.par
            end_dt = self.maturity_dt
        else:
            sell_price = self.dirty_price_from_ytm(end_dt, end_ytm, convention)

        holding_period = (end_dt - begin_dt) / G_DAYS_IN_YEAR

        pnl = sell_price - buy_price
        simple_rtn = pnl / buy_price / holding_period
        irr = (sell_price / buy_price) ** (1.0 / holding_period) - 1.0

        return simple_rtn, irr, pnl

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("COUPON (%)", 0)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("TIME DAY COUNT TYPE", self.time_dc_type)
        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""
        print(self)


########################################################################################
