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
                 coupon: float,
                 freq_type: FrequencyTypes,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 day_count_convention_type: DayCountTypes = DayCountTypes.ACT_360,
                 face: float = 100.0):

        check_argument_types(self.__init__, locals())

        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._frequency = annual_frequency(freq_type)

        # ISDA Style conventions
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type
        self._day_count_convention_type = day_count_convention_type

        self._face = face
        self._par = 100.0

        self._flow_dates = []
        self._settlement_date = Date(1, 1, 1900)
        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

    ###########################################################################

    def clean_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cash flows. """

        full_price = self.full_price_from_discount_curve(settlement_date,
                                                         discount_curve)
        accrued = self._accrued_interest * self._par / self._face
        clean_price = full_price - accrued
        return clean_price

    ###########################################################################

    def full_price_from_discount_curve(self,
                                       settlement_date: Date,
                                       discount_curve: DiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cash flows. """

        self.calculate_payments(settlement_date)
        pv = 0.0

        num_flows = len(self._flow_dates)

        for i in range(1, num_flows):
            dt = self._flow_dates[i]
            df = discount_curve.df(dt)
            flow = self._flow_amounts[i]
            pv = pv + flow * df

        return pv * self._par / self._face

    ###########################################################################

    def calculate_payments(self,
                           settlement_date: Date):

        # No need to generate flows if settlement date has not changed
        if settlement_date == self._settlement_date:
            return

        if settlement_date == self._maturity_date:
            raise FinError("Settlement date is maturity date.")

        self._settlement_date = settlement_date
        calendar_type = CalendarTypes.NONE
        bus_day_rule_type = BusDayAdjustTypes.NONE
        date_gen_rule_type = DateGenRuleTypes.BACKWARD

        self._flow_dates = Schedule(settlement_date,
                                    self._maturity_date,
                                    self._freq_type,
                                    calendar_type,
                                    bus_day_rule_type,
                                    date_gen_rule_type)._generate()

        self._pcd = self._flow_dates[0]
        self._ncd = self._flow_dates[1]
        self.calc_accrued_interest(settlement_date)

        self._flow_amounts = [0.0]
        basis = DayCount(self._day_count_convention_type)

        prev_dt = self._pcd

        for next_dt in self._flow_dates[1:]:
            alpha = basis.year_frac(prev_dt, next_dt)[0]
            flow = self._coupon * alpha * self._face
            self._flow_amounts.append(flow)
            prev_dt = next_dt

    ###########################################################################

    def calc_accrued_interest(self,
                              settlement_date: Date):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. """

        if settlement_date != self._settlement_date:
            self.calculate_payments(settlement_date)

        if len(self._flow_dates) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self._day_count_convention_type)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd,
                                            self._frequency)

        self._alpha = 1.0 - acc_factor * self._frequency

        self._accrued_interest = acc_factor * self._face * self._coupon
        self._accrued_days = num
        return self._accrued_interest

    ###########################################################################

    def print_flows(self,
                    settlement_date: Date):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        self.calculate_payments(settlement_date)

        num_flows = len(self._flow_dates)
        for i in range(1, num_flows):
            dt = self._flow_dates[i]
            flow = self._flow_amounts[i]
            print(dt, ",", flow)

    ###########################################################################

    def __repr__(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUS_DAY_RULE", self._bus_day_adjust_type)
        s += label_to_string("DATE_GEN_RULE", self._date_gen_rule_type)

        return s

    ###########################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
