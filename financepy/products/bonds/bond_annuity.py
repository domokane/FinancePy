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
    """ An annuity is a vector of dates and flows generated according to ISDA
    standard rules which starts on the next date after the start date
    (effective date) and runs up to an end date with no principal repayment.
    Dates are then adjusted according to a specified calendar. """

    def __init__(self,
                 maturity_date: Date,
                 cpn: float,
                 freq_type: FrequencyTypes,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 day_count_convention_type: DayCountTypes = DayCountTypes.ACT_360):

        check_argument_types(self.__init__, locals())

        self._maturity_date = maturity_date
        self._cpn = cpn
        self._freq_type = freq_type
        self._freq = annual_frequency(freq_type)

        # ISDA Style conventions
        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type
        self._day_count_convention_type = day_count_convention_type

        self._par = 100.0

        self._cpn_dates = []
        self._settle_date = Date(1, 1, 1900)
        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

    ###########################################################################

    def clean_price_from_discount_curve(self,
                                        settle_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cash flows. """

        dirty_price = self.dirty_price_from_discount_curve(settle_date,
                                                           discount_curve)
        accrued = self._accrued_interest * self._par
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def dirty_price_from_discount_curve(self,
                                        settle_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cash flows. """

        self.calculate_payments(settle_date, 1.0)
        pv = 0.0

        num_flows = len(self._cpn_dates)

        for i in range(1, num_flows):
            dt = self._cpn_dates[i]
            df = discount_curve.df(dt)
            flow = self._flow_amounts[i]
            pv = pv + flow * df

        return pv * self._par

    ###########################################################################

    def calculate_payments(self,
                           settle_date: Date,
                           face: (float)):

        # No need to generate flows if settlement date has not changed
        if settle_date == self._settle_date:
            return

        if settle_date == self._maturity_date:
            raise FinError("Settlement date is maturity date.")

        self._settle_date = settle_date
        bd_type = BusDayAdjustTypes.FOLLOWING
        dg_type = DateGenRuleTypes.BACKWARD

        self._cpn_dates = Schedule(settle_date,
                                   self._maturity_date,
                                   self._freq_type,
                                   self._cal_type,
                                   bd_type,
                                   dg_type)._generate()

        self._pcd = self._cpn_dates[0]
        self._ncd = self._cpn_dates[1]
        self.accrued_interest(settle_date, 1.0)

        self._flow_amounts = [0.0]
        basis = DayCount(self._day_count_convention_type)

        prev_dt = self._pcd

        for next_dt in self._cpn_dates[1:]:
            alpha = basis.year_frac(prev_dt, next_dt)[0]
            flow = self._cpn * alpha * face
            self._flow_amounts.append(flow)
            prev_dt = next_dt

    ###########################################################################

    def accrued_interest(self,
                         settle_date: Date,
                         face: (float)):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. """

        if settle_date != self._settle_date:
            self.calculate_payments(settle_date, 1.0)

        if len(self._cpn_dates) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self._day_count_convention_type)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settle_date,
                                            self._ncd,
                                            self._freq)

        self._alpha = 1.0 - acc_factor * self._freq

        self._accrued_interest = acc_factor * face * self._cpn
        self._accrued_days = num
        return self._accrued_interest

    ###########################################################################

    def print_payments(self,
                       settle_date: Date,
                       face: (float)):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        self.calculate_payments(settle_date, face)

        num_flows = len(self._cpn_dates)
        for i in range(1, num_flows):
            dt = self._cpn_dates[i]
            flow = self._flow_amounts[i]
            print(dt, ",", flow)

    ###########################################################################

    def __repr__(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUS_DAY_RULE", self._bd_type)
        s += label_to_string("DATE_GEN_RULE", self._dg_type)

        return s

    ###########################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
