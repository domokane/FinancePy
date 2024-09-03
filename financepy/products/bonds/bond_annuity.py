##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.date import Date
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.schedule import Schedule
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.error import FinError
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve


###############################################################################


class BondAnnuity:
    """An annuity is a vector of dates and flows generated according to ISDA
    standard rules which starts on the next date after the start date
    (effective date) and runs up to an end date with no principal repayment.
    Dates are then adjusted according to a specified calendar."""

    def __init__(
        self,
        maturity_dt: Date,
        cpn: float,
        freq_type: FrequencyTypes,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
    ):

        check_argument_types(self.__init__, locals())

        self.maturity_dt = maturity_dt
        self.cpn = cpn
        self.freq_type = freq_type
        self.freq = annual_frequency(freq_type)

        # ISDA Style conventions
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type
        self.dc_type = dc_type

        self.par = 100.0

        self.cpn_dts = []
        self.settle_dt = Date(1, 1, 1900)
        self.accrued_int = None
        self.accrued_days = 0.0
        self.alpha = 0.0

        self.pcd = None
        self.ncd = None
        self.flow_amounts = None

    ###########################################################################

    def clean_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        """Calculate the bond price using some discount curve to present-value
        the bond's cash flows."""

        dirty_price = self.dirty_price_from_discount_curve(
            settle_dt, discount_curve
        )
        accrued = self.accrued_int * self.par
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def dirty_price_from_discount_curve(
        self, settle_dt: Date, discount_curve: DiscountCurve
    ):
        """Calculate the bond price using some discount curve to present-value
        the bond's cash flows."""

        self.calculate_payments(settle_dt, 1.0)
        pv = 0.0

        num_flows = len(self.cpn_dts)

        for i in range(1, num_flows):
            dt = self.cpn_dts[i]
            df = discount_curve.df(dt)
            flow = self.flow_amounts[i]
            pv = pv + flow * df

        return pv * self.par

    ###########################################################################

    def calculate_payments(self, settle_dt: Date, face: float):
        """Calculate bond payments"""
        # No need to generate flows if settlement date has not changed
        if settle_dt == self.settle_dt:
            return

        if settle_dt == self.maturity_dt:
            raise FinError("Settlement date is maturity date.")

        self.settle_dt = settle_dt
        bd_type = BusDayAdjustTypes.FOLLOWING
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
        self.accrued_interest(settle_dt, 1.0)

        self.flow_amounts = [0.0]
        basis = DayCount(self.dc_type)

        prev_dt = self.pcd

        for next_dt in self.cpn_dts[1:]:
            alpha = basis.year_frac(prev_dt, next_dt)[0]
            flow = self.cpn * alpha * face
            self.flow_amounts.append(flow)
            prev_dt = next_dt

    ###########################################################################

    def accrued_interest(self, settle_dt: Date, face: float):
        """Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date."""

        if settle_dt != self.settle_dt:
            self.calculate_payments(settle_dt, 1.0)

        if len(self.cpn_dts) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self.dc_type)

        (acc_factor, num, _) = dc.year_frac(
            self.pcd, settle_dt, self.ncd, self.freq
        )

        self.alpha = 1.0 - acc_factor * self.freq

        self.accrual_factor = acc_factor
        self.accrued_int = acc_factor * face * self.cpn
        self.accrued_days = num

        return self.accrued_int

    ###########################################################################

    def print_payments(self, settle_dt: Date, face: float):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""

        self.calculate_payments(settle_dt, face)

        num_flows = len(self.cpn_dts)
        for i in range(1, num_flows):
            dt = self.cpn_dts[i]
            flow = self.flow_amounts[i]
            print(dt, ",", flow)

    ###########################################################################

    def __repr__(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS_DAY_RULE", self.bd_type)
        s += label_to_string("DATE_GEN_RULE", self.dg_type)

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
