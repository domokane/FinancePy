##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.schedule import Schedule
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.helpers import label_to_string, check_argument_types


###############################################################################
# TODO: Need to complete and verify the risk sensitivity calculations.
###############################################################################


def _f(dm, *args):
    """Function used to do solve root search in DM calculation"""

    self = args[0]
    settle_dt = args[1]
    next_cpn = args[2]
    current_ibor = args[3]
    future_ibor = args[4]
    dirty_price = args[5]

    px = self.dirty_price_from_dm(
        settle_dt, next_cpn, current_ibor, future_ibor, dm
    )

    obj_fn = px - dirty_price
    return obj_fn


###############################################################################


class BondFRN:
    """Class for managing floating rate notes that pay a floating index plus a
    quoted margin."""

    def __init__(
        self,
        issue_dt: Date,
        maturity_dt: Date,
        quoted_margin: float,  # Fixed spread paid on top of index
        freq_type: FrequencyTypes,
        dc_type: DayCountTypes,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
    ):
        """Create FinFloatingRateNote object given its maturity date, its
        quoted margin, coupon frequency, DAY COUNT TYPE. Face is the size of
        the position and par is the notional on which price is quoted."""

        check_argument_types(self.__init__, locals())

        self.issue_dt = issue_dt
        self.maturity_dt = maturity_dt
        self.quoted_margin = quoted_margin
        self.freq_type = freq_type
        self.dc_type = dc_type
        self.freq = annual_frequency(freq_type)
        self.cal_type = cal_type

        self.settle_dt = Date(1, 1, 1900)
        self.accrued = None
        self.accrued_days = 0.0
        self.accrual_factor = 0.0
        self.accrued_int = 0.0

        self.cpn_dts = []
        self.flow_amounts = []
        self.par = 100.0  # This is how price is quoted

        self._calculate_cpn_dts()

    ###########################################################################

    def _calculate_cpn_dts(self):
        """Determine the bond cashflow payment dates."""

        # This should only be called once from init

        bd_type = BusDayAdjustTypes.NONE
        dg_type = DateGenRuleTypes.BACKWARD

        self.cpn_dts = Schedule(
            self.issue_dt,
            self.maturity_dt,
            self.freq_type,
            self.cal_type,
            bd_type,
            dg_type,
        ).generate()

    ###########################################################################

    def dirty_price_from_dm(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the full price of the bond from its discount margin (DM)
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally, there
        is the level of subsequent future Ibor payments and the discount
        margin."""

        day_counter = DayCount(self.dc_type)

        q = self.quoted_margin
        num_flows = len(self.cpn_dts)

        # We discount using Libor over the period from settlement to the ncd
        (alpha, _, _) = day_counter.year_frac(settle_dt, self.ncd)

        df = 1.0 / (1.0 + alpha * (current_ibor + dm))

        # A full coupon is paid
        (alpha, _, _) = day_counter.year_frac(self.pcd, self.ncd)
        pv = next_cpn * alpha * df

        # Now do all subsequent coupons that fall after the ncd
        for i_flow in range(1, num_flows):

            if self.cpn_dts[i_flow] > self.ncd:

                pcd = self.cpn_dts[i_flow - 1]
                ncd = self.cpn_dts[i_flow]
                (alpha, _, _) = day_counter.year_frac(pcd, ncd)

                df = df / (1.0 + alpha * (future_ibor + dm))
                c = future_ibor + q
                pv = pv + c * alpha * df

        pv += df
        pv = pv * self.par
        return pv

    ###########################################################################

    def principal(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
        face: float = 100.0,
    ):
        """Calculate the clean trade price of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates."""

        dirty_price = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        self.accrued = self.accrual_factor * next_cpn * 1.0
        principal = dirty_price * face / self.par - self.accrued
        return principal

    ###########################################################################

    def dollar_duration(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg."""

        dy = 0.0001  # 1 basis point

        p0 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor + dy, future_ibor, dm
        )

        p2 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor - dy, future_ibor, dm
        )

        durn = (p2 - p0) / dy / 2.0
        return durn

    ###########################################################################

    def dollar_credit_duration(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the risk or dP/dy of the bond by bumping."""

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dy = 0.0001

        p0 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm + dy
        )

        p2 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        durn = (p2 - p0) / dy
        return durn

    ###########################################################################

    def macauley_duration(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the Macauley duration of the FRN on a settlement date
        given its yield to maturity."""

        dd = self.dollar_duration(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        fp = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        md = dd * (1.0 + (next_cpn + dm) / self.freq) / fp
        return md

    ###########################################################################

    def modified_duration(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally, there
        is the level of subsequent future Ibor payments and the discount
        margin."""

        dd = self.dollar_duration(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        fp = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )
        md = dd / fp
        return md

    ###########################################################################

    def modified_credit_duration(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally, there
        is the level of subsequent future Ibor payments and the discount
        margin."""

        dd = self.dollar_credit_duration(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        fp = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )
        md = dd / fp
        return md

    ###########################################################################

    def convexity_from_dm(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the bond convexity from the discount margin (DM) using a
        standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally, there
        is the level of subsequent future Ibor payments and the discount
        margin."""

        dy = 0.0001

        p0 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor - dy, future_ibor, dm
        )

        p1 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        p2 = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor + dy, future_ibor, dm
        )

        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self.par
        return conv

    ###########################################################################

    def clean_price_from_dm(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        dm: float,
    ):
        """Calculate the bond clean price from the discount margin
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally, there
        is the level of subsequent future Ibor payments and the discount
        margin."""

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dirty_price = self.dirty_price_from_dm(
            settle_dt, next_cpn, current_ibor, future_ibor, dm
        )

        self.accrued_interest(settle_dt, next_cpn)

        self.accrued = self.accrual_factor * next_cpn * self.par

        clean_price = dirty_price - self.accrued
        return clean_price

    ###########################################################################

    def discount_margin(
        self,
        settle_dt: Date,
        next_cpn: float,
        current_ibor: float,
        future_ibor: float,
        clean_price: float,
    ):
        """Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver."""

        self.accrued_interest(settle_dt, next_cpn)

        # Accrued needs to be adjusted to par notional
        self.accrued = self.accrual_factor * next_cpn * self.par

        dirty_price = clean_price + self.accrued

        argtuple = (
            self,
            settle_dt,
            next_cpn,
            current_ibor,
            future_ibor,
            dirty_price,
        )

        dm = optimize.newton(
            _f,
            x0=0.01,  # initial value of 10%
            fprime=None,
            args=argtuple,
            tol=1e-12,
            maxiter=50,
            fprime2=None,
        )

        return dm

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, next_cpn: float):
        """Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Ex-dividend dates are
        not handled. Contact me if you need this functionality."""

        num_flows = len(self.cpn_dts)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self.dc_type)

        for i in range(1, num_flows):
            if self.cpn_dts[i] > settle_dt:
                self.pcd = self.cpn_dts[i - 1]
                self.ncd = self.cpn_dts[i]
                break

        (acc_factor, num, _) = dc.year_frac(
            self.pcd, settle_dt, self.ncd, self.freq_type
        )

        self.alpha = 1.0 - acc_factor * self.freq
        self.accrual_factor = acc_factor
        self.accrued_days = num
        self.accrued_int = acc_factor * next_cpn * 1.0

        return self.accrued_int

    ###########################################################################

    def print_payments(self, settle_dt: Date):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        self._calculate_cpn_dts()
        for dt in self.cpn_dts[1:-1]:
            print(dt)

        print(self.cpn_dts[-1])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self.issue_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string(
            "QUOTED MARGIN (bp)", self.quoted_margin * 10000.0
        )
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
