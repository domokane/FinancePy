##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import ONE_MILLION
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes,  DateGenRuleTypes
from ...utils.calendar import Calendar, BusDayAdjustTypes
from ...utils.schedule import Schedule
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.global_types import SwapTypes
from ...market.curves.discount_curve import DiscountCurve

##########################################################################


class SwapFixedLeg:
    """ Class for managing the fixed leg of a swap. A fixed leg is a leg with
    a sequence of flows calculated according to an ISDA schedule and with a 
    coupon that is fixed over the life of the swap. """

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 end_date: (Date, str),  # Date contract ends
                 leg_type: SwapTypes,
                 coupon: (float),
                 freq_type: FrequencyTypes,
                 day_count_type: DayCountTypes,
                 notional: float = ONE_MILLION,
                 principal: float = 0.0,
                 payment_lag: int = 0,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 end_of_month: bool = False):
        """ Create the fixed leg of a swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional.  """

        check_argument_types(self.__init__, locals())

        if type(end_date) == Date:
            self._termination_date = end_date
        else:
            self._termination_date = effective_date.add_tenor(end_date)

        calendar = Calendar(calendar_type)

        self._maturity_date = calendar.adjust(self._termination_date,
                                              bus_day_adjust_type)

        if effective_date > self._maturity_date:
            raise FinError("Effective date after maturity date")

        self._effective_date = effective_date
        self._end_date = end_date
        self._leg_type = leg_type
        self._freq_type = freq_type
        self._payment_lag = payment_lag
        self._notional = notional
        self._principal = principal
        self._coupon = coupon

        self._day_count_type = day_count_type
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type
        self._end_of_month = end_of_month

        self._startAccruedDates = []
        self._endAccruedDates = []
        self._payment_dates = []
        self._payments = []
        self._year_fracs = []
        self._accrued_days = []
        self._rates = []

        self.generate_payments()

###############################################################################

    def generate_payments(self):
        ''' These are generated immediately as they are for the entire
        life of the swap. Given a valuation date we can determine
        which cash flows are in the future and value the swap
        The schedule allows for a specified lag in the payment date
        Nothing is paid on the swap effective date and so the first payment
        date is the first actual payment date. '''

        schedule = Schedule(self._effective_date,
                            self._termination_date,
                            self._freq_type,
                            self._calendar_type,
                            self._bus_day_adjust_type,
                            self._date_gen_rule_type,
                            end_of_month=self._end_of_month)

        scheduleDates = schedule._adjusted_dates

        if len(scheduleDates) < 2:
            raise FinError("Schedule has none or only one date")

        self._startAccruedDates = []
        self._endAccruedDates = []
        self._payment_dates = []
        self._payments = []
        self._year_fracs = []
        self._accrued_days = []
        self._rates = []

        prev_dt = scheduleDates[0]

        day_counter = DayCount(self._day_count_type)
        calendar = Calendar(self._calendar_type)

        for next_dt in scheduleDates[1:]:

            self._startAccruedDates.append(prev_dt)
            self._endAccruedDates.append(next_dt)

            if self._payment_lag == 0:
                payment_date = next_dt
            else:
                payment_date = calendar.add_business_days(next_dt,
                                                          self._payment_lag)

            self._payment_dates.append(payment_date)

            (year_frac, num, den) = day_counter.year_frac(prev_dt,
                                                          next_dt)

            self._rates.append(self._coupon)

            payment = year_frac * self._notional * self._coupon

            self._payments.append(payment)
            self._year_fracs.append(year_frac)
            self._accrued_days.append(num)

            prev_dt = next_dt

###############################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve):

        self._paymentDfs = []
        self._paymentPVs = []
        self._cumulativePVs = []

        notional = self._notional
        dfValue = discount_curve.df(valuation_date)
        legPV = 0.0
        numPayments = len(self._payment_dates)

        dfPmnt = 0.0

        for iPmnt in range(0, numPayments):

            pmntDate = self._payment_dates[iPmnt]
            pmntAmount = self._payments[iPmnt]

            if pmntDate > valuation_date:

                dfPmnt = discount_curve.df(pmntDate) / dfValue
                pmntPV = pmntAmount * dfPmnt
                legPV += pmntPV

                self._paymentDfs.append(dfPmnt)
                self._paymentPVs.append(pmntAmount*dfPmnt)
                self._cumulativePVs.append(legPV)

            else:

                self._paymentDfs.append(0.0)
                self._paymentPVs.append(0.0)
                self._cumulativePVs.append(0.0)

        if pmntDate > valuation_date:
            paymentPV = self._principal * dfPmnt * notional
            self._paymentPVs[-1] += paymentPV
            legPV += paymentPV
            self._cumulativePVs[-1] = legPV

        if self._leg_type == SwapTypes.PAY:
            legPV = legPV * (-1.0)

        return legPV

##########################################################################

    def print_payments(self):
        """ Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("COUPON (%):", self._coupon * 100)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))

        if len(self._payments) == 0:
            print("Payments not calculated.")
            return

        header = "PAY_DATE     ACCR_START   ACCR_END      DAYS  YEARFRAC"
        header += "    RATE      PAYMENT "
        print(header)

        num_flows = len(self._payment_dates)

        for iFlow in range(0, num_flows):
            print("%11s  %11s  %11s  %4d  %8.6f  %8.6f  %11.2f" %
                  (self._payment_dates[iFlow],
                   self._startAccruedDates[iFlow],
                   self._endAccruedDates[iFlow],
                   self._accrued_days[iFlow],
                   self._year_fracs[iFlow],
                   self._rates[iFlow] * 100.0,
                   self._payments[iFlow]))

###############################################################################

    def print_valuation(self):
        """ Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("COUPON (%):", self._coupon * 100)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))

        if len(self._payments) == 0:
            print("Payments not calculated.")
            return

        header = "PAY_DATE     ACCR_START   ACCR_END     DAYS  YEARFRAC"
        header += "    RATE      PAYMENT       DF          PV        CUM PV"
        print(header)

        num_flows = len(self._payment_dates)

        for iFlow in range(0, num_flows):
            print("%11s  %11s  %11s  %4d  %8.6f  %8.5f  % 11.2f  %10.8f  % 11.2f  % 11.2f" %
                  (self._payment_dates[iFlow],
                   self._startAccruedDates[iFlow],
                   self._endAccruedDates[iFlow],
                   self._accrued_days[iFlow],
                   self._year_fracs[iFlow],
                   self._rates[iFlow] * 100.0,
                   self._payments[iFlow],
                   self._paymentDfs[iFlow],
                   self._paymentPVs[iFlow],
                   self._cumulativePVs[iFlow]))

##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._effective_date)
        s += label_to_string("TERMINATION DATE", self._termination_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("PRINCIPAL", self._principal)
        s += label_to_string("LEG TYPE", self._leg_type)
        s += label_to_string("COUPON", self._coupon)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAY COUNT", self._day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUS DAY ADJUST", self._bus_day_adjust_type)
        s += label_to_string("DATE GEN TYPE", self._date_gen_rule_type)
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
