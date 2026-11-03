##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################


import numpy as np
from scipy import optimize

from ...market.curves.discount_curve import DiscountCurve
from ...products.bonds.bond import (
    YTMCalcType,
    _validate_price,
    _validate_yield,
    _vectorise_price,
)
from ...products.product import Product
from ...utils.calendar import BusDayAdjustTypes, Calendar, CalendarTypes
from ...utils.check_values import check_curve_dt
from ...utils.date import Date
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.error import FinError
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...utils.helpers import check_argument_types, label_to_string


def _f(ytm, *args):
    """Objective function used for root-solving in yield to maturity."""
    bond = args[0]
    settle_dt = args[1]
    price = args[2]
    convention = args[3]
    px = bond.dirty_price_from_ytm(settle_dt, ytm, convention)
    return px - price


class BondPerpetual(Product):
    """Class for perpetual bonds (Consols / undated debt). These instruments
    have regular periodic coupon payments of a known rate and have no maturity
    date (T = infinity). The principal is never repaid."""

    def __init__(
        self,
        issue_dt: Date,
        coupon: float,
        freq_type: FrequencyTypes,
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_ACT_ICMA,
        first_coupon_dt: Date | None = None,
        ex_div_days: int = 0,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
    ) -> None:
        """Create a BondPerpetual object by providing issue date, annualized
        coupon, payment frequency, accrual day count convention, and optional
        first coupon date, ex-dividend days, calendar, and business day rule."""

        check_argument_types(self.__init__, locals())

        if coupon <= 0.0:
            raise FinError("Coupon must be strictly positive for a perpetual bond.")

        if freq_type == FrequencyTypes.ZERO:
            raise FinError("Zero coupon perpetual bond is not valid.")

        if ex_div_days > 90:
            raise FinError("Ex-dividend days cannot exceed 90, got " + str(ex_div_days))

        self.freq_type = freq_type
        self.freq = annual_frequency(freq_type)

        if first_coupon_dt is not None:
            if first_coupon_dt <= issue_dt:
                raise FinError("First coupon date must be strictly after issue date.")
            self.first_coupon_dt = first_coupon_dt
        else:
            months_step = 12 // self.freq
            self.first_coupon_dt = issue_dt.add_months(months_step)

        self.issue_dt = issue_dt
        self.coupon = float(coupon)
        self.cpn = float(coupon)
        self.accrual_dc_type = accrual_dc_type
        self.ex_div_days = ex_div_days
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.par = 100.0

        self.ex_div_dt = None
        self._pcd = None
        self._ncd = None
        self.alpha = 0.0
        self.accrued_int = 0.0
        self.accrued_days = 0

    ###########################################################################

    def _calc_pcd_ncd(self, settle_dt: Date) -> None:
        """Calculate Previous Coupon Date (PCD) and Next Coupon Date (NCD) for
        a given settlement date."""
        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        if settle_dt < self.first_coupon_dt:
            self._pcd = self.issue_dt
            self._ncd = self.first_coupon_dt
            return

        months_step = 12 // self.freq
        months_diff = (settle_dt.y - self.first_coupon_dt.y) * 12 + (
            settle_dt.m - self.first_coupon_dt.m
        )
        k = months_diff // months_step

        cand_pcd = self.first_coupon_dt.add_months(k * months_step)
        if cand_pcd > settle_dt:
            k -= 1
            cand_pcd = self.first_coupon_dt.add_months(k * months_step)

        while True:
            next_cand = self.first_coupon_dt.add_months((k + 1) * months_step)
            if next_cand <= settle_dt:
                k += 1
                cand_pcd = next_cand
            else:
                cand_ncd = next_cand
                break

        self._pcd = cand_pcd
        self._ncd = cand_ncd

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float = 100.0) -> float:
        """Calculate the amount of coupon accrued between the previous coupon
        date and the settlement date."""
        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        self._calc_pcd_ncd(settle_dt)

        dc = DayCount(self.accrual_dc_type)
        cal = Calendar(self.cal_type)

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

    def dirty_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float | list[float] | np.ndarray,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float | np.ndarray:
        """Calculate the dirty price of a perpetual bond from its yield to
        maturity. Supports vectorization with respect to yield inputs."""
        _validate_yield(ytm)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        if convention not in YTMCalcType:
            raise FinError("Yield convention unknown: " + str(convention))

        self.accrued_interest(settle_dt, 1.0)

        pay_first_cpn = 1.0
        if settle_dt > self.ex_div_dt:
            pay_first_cpn = 0.0

        is_scalar = np.isscalar(ytm)
        ytm_arr = np.atleast_1d(np.array(ytm, dtype=float))

        if np.any(ytm_arr <= 0.0):
            raise FinError("Yield to maturity must be strictly positive for a perpetual bond.")

        f = self.freq
        c = self.cpn
        v = 1.0 / (1.0 + ytm_arr / f)

        # In a perpetual bond, the remaining coupons form an infinite annuity:
        # sum_{k=1}^infinity (c/f) * v^k = (c/f) * v / (1 - v) = c / ytm
        cash_at_ncd = (c / f) * pay_first_cpn + (c / ytm_arr)

        if convention in (
            YTMCalcType.UK_DMO,
            YTMCalcType.US_STREET,
            YTMCalcType.CALCULUS,
            YTMCalcType.CFETS,
        ):
            dp = (v ** self.alpha) * cash_at_ncd
        elif convention == YTMCalcType.US_TREASURY:
            # 31 CFR 356 Appendix B fractional discounting
            vw = 1.0 / (1.0 + self.alpha * ytm_arr / f)
            dp = vw * cash_at_ncd
        else:
            dp = (v ** self.alpha) * cash_at_ncd

        res = dp * self.par
        if is_scalar:
            return float(res[0])
        return res

    ###########################################################################

    def clean_price_from_ytm(
        self,
        settle_dt: Date,
        ytm: float | list[float] | np.ndarray,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float | np.ndarray:
        """Calculate the clean price of a perpetual bond from its yield to
        maturity."""
        _validate_yield(ytm)
        dirty_price = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        accrued = self.accrued_interest(settle_dt, self.par)
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def yield_to_maturity(
        self,
        settle_dt: Date,
        clean_price: float | list[float] | np.ndarray,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float | np.ndarray:
        """Calculate the yield to maturity (yield to perpetuity) from clean price."""
        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        clean_prices = _vectorise_price(clean_price)
        _validate_price(clean_prices[0])

        self.accrued_interest(settle_dt, 1.0)
        accrued_amount = self.accrued_int * self.par
        dirty_prices = clean_prices + accrued_amount

        ytms = []
        for dirty_price in dirty_prices:
            if dirty_price <= 0.0:
                raise FinError("Dirty price must be strictly positive for yield calculation.")

            argtuple = (self, settle_dt, dirty_price, convention)

            try:
                ytm = optimize.brentq(
                    _f,
                    a=1e-6,
                    b=2.0,
                    args=argtuple,
                    xtol=1e-8,
                    maxiter=100,
                )
            except RuntimeError:
                print(f"Warning: YTM calculation did not converge for price {dirty_price}")
                ytm = np.nan

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        return np.array(ytms)

    ###########################################################################

    def current_yield(
        self, settle_dt: Date, clean_price: float | list[float] | np.ndarray
    ) -> float | np.ndarray:
        """Calculate the bond's current yield: coupon / clean_price."""
        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        clean_prices = _vectorise_price(clean_price)
        _validate_price(clean_prices[0])

        current_ys = (self.cpn * self.par) / clean_prices
        if len(current_ys) == 1:
            return float(current_ys[0])
        return current_ys

    ###########################################################################

    def dv01(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float:
        """Calculate the DV01 (dollar value of 1 basis point increase in yield)."""
        _validate_yield(ytm)
        dy = 0.0001  # 1 bp
        p0 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p1 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        return float(p0 - p1)

    ###########################################################################

    def dollar_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float:
        """Calculate dollar duration: 10,000 * DV01."""
        return 10000.0 * self.dv01(settle_dt, ytm, convention)

    ###########################################################################

    def modified_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float:
        """Calculate modified duration: dollar_duration / dirty_price."""
        dd = self.dollar_duration(settle_dt, ytm, convention)
        dp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        return float(dd / dp)

    ###########################################################################

    def macaulay_duration(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float:
        """Calculate Macaulay duration: modified_duration * (1 + ytm / freq)."""
        d = 1.0 + ytm / self.freq
        dd = self.dollar_duration(settle_dt, ytm, convention)
        dp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        return float((dd / dp) * d)

    ###########################################################################

    def convexity_from_ytm(
        self,
        settle_dt: Date,
        ytm: float,
        convention: YTMCalcType = YTMCalcType.UK_DMO,
    ) -> float:
        """Calculate convexity from yield to maturity: (1 / P) * (d^2 P / dy^2)."""
        _validate_yield(ytm)
        dy = 0.0001
        p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
        p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / (dy * dy * p1)
        return float(conv)

    ###########################################################################

    def dirty_price_from_discount_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        num_coupons: int = 200,
    ) -> float:
        """Calculate the dirty price of the perpetual bond from an arbitrary
        discount curve using coupon stream summation with a terminal perpetuity
        tail adjustment."""
        check_curve_dt(settle_dt, discount_curve)

        if settle_dt < self.issue_dt:
            raise FinError("Settlement date falls before issue date")

        self._calc_pcd_ncd(settle_dt)
        cal = Calendar(self.cal_type)
        self.ex_div_dt = cal.add_business_days(self._ncd, -self.ex_div_days)

        months_step = 12 // self.freq
        cpn_amount = (self.cpn / self.freq)

        px_pv = 0.0
        last_dt = self._ncd

        for i in range(num_coupons):
            flow_dt = self._ncd.add_months(i * months_step)
            last_dt = flow_dt
            if flow_dt > settle_dt:
                df = discount_curve.df(flow_dt)
                if i == 0 and settle_dt > self.ex_div_dt:
                    flow = 0.0
                else:
                    flow = cpn_amount
                px_pv += flow * df

        # Terminal perpetuity tail beyond num_coupons
        df_terminal = discount_curve.df(last_dt)
        r_terminal = discount_curve.zero_rate_cc(last_dt)
        if r_terminal > 1e-6:
            y_terminal = self.freq * (np.exp(r_terminal / self.freq) - 1.0)
            tail_pv = (cpn_amount / (y_terminal / self.freq)) * df_terminal
            px_pv += tail_pv

        df_settle = discount_curve.df(settle_dt)
        dirty_price = (px_pv / df_settle) * self.par
        return float(dirty_price)

    ###########################################################################

    def clean_price_from_discount_curve(
        self,
        settle_dt: Date,
        discount_curve: DiscountCurve,
        num_coupons: int = 200,
    ) -> float:
        """Calculate the clean price of the perpetual bond from a discount curve."""
        check_curve_dt(settle_dt, discount_curve)
        dirty_price = self.dirty_price_from_discount_curve(
            settle_dt, discount_curve, num_coupons
        )
        accrued = self.accrued_interest(settle_dt, self.par)
        return float(dirty_price - accrued)

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        discount_curve: DiscountCurve,
        *args,
        **kwargs,
    ) -> float:
        """Evaluate the clean price of the perpetual bond on value_dt using
        the provided discount curve, adhering to the Product interface."""
        return self.clean_price_from_discount_curve(value_dt, discount_curve)

    ###########################################################################

    def __repr__(self) -> str:
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("FIRST COUPON DATE", self.first_coupon_dt)
        s += label_to_string("COUPON", self.cpn)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("ACCRUAL DAY COUNT", self.accrual_dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST", self.bd_type)
        s += label_to_string("EX DIVIDEND DAYS", self.ex_div_days)
        return s

    def _print(self) -> None:
        """Print the perpetual bond details."""
        print(self)
