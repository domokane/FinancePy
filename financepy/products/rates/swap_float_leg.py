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


class SwapFloatLeg:
    """ Class for managing the floating leg of a swap. A float leg consists of
    a sequence of flows calculated according to an ISDA schedule and with a 
    coupon determined by an index curve which changes over life of the swap."""

    def __init__(self,
                 effective_date: Date,  # Date interest starts to accrue
                 end_date: (Date, str),  # Date contract ends
                 leg_type: SwapTypes,
                 spread: (float),
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
            raise FinError("Start date after maturity date")

        self._effective_date = effective_date
        self._end_date = end_date
        self._leg_type = leg_type
        self._freq_type = freq_type
        self._payment_lag = payment_lag
        self._notional = notional
        self._principal = 0.0
        self._spread = spread

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

        self.generate_payment_dates()

###############################################################################

    def generate_payment_dates(self):
        """ Generate the floating leg payment dates and accrual factors. The
        coupons cannot be generated yet as we do not have the index curve. """

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
        self._year_fracs = []
        self._accrued_days = []

        prev_dt = scheduleDates[0]

        day_counter = DayCount(self._day_count_type)
        calendar = Calendar(self._calendar_type)

        # All of the lists end up with the same length
        for next_dt in scheduleDates[1:]:

            self._startAccruedDates.append(prev_dt)
            self._endAccruedDates.append(next_dt)

            if self._payment_lag == 0:
                payment_date = next_dt
            else:
                payment_date = calendar.add_business_days(next_dt,
                                                          self._payment_lag)

            self._payment_dates.append(payment_date)

            (year_frac, num, _) = day_counter.year_frac(prev_dt,
                                                        next_dt)

            self._year_fracs.append(year_frac)
            self._accrued_days.append(num)

            prev_dt = next_dt

###############################################################################

    def value(self,
              valuation_date: Date,  # This should be the settlement date
              discount_curve: DiscountCurve,
              index_curve: DiscountCurve,
              firstFixingRate: float = None):
        """ Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve as of the valuation date
        supplied. For an existing swap, the user must enter the next fixing
        coupon. """

        if discount_curve is None:
            raise FinError("Discount curve is None")

        if index_curve is None:
            index_curve = discount_curve

        self._rates = []
        self._payments = []
        self._paymentDfs = []
        self._paymentPVs = []
        self._cumulativePVs = []

        notional = self._notional
        dfValue = discount_curve.df(valuation_date)
        legPV = 0.0
        numPayments = len(self._payment_dates)
        firstPayment = False

        index_basis = index_curve._day_count_type
        index_day_counter = DayCount(index_basis)

        for iPmnt in range(0, numPayments):

            pmntDate = self._payment_dates[iPmnt]

            if pmntDate > valuation_date:

                startAccruedDt = self._startAccruedDates[iPmnt]
                endAccruedDt = self._endAccruedDates[iPmnt]
                pay_alpha = self._year_fracs[iPmnt]

                (index_alpha, num, _) = index_day_counter.year_frac(startAccruedDt,
                                                                    endAccruedDt)

                if firstPayment is False and firstFixingRate is not None:

                    fwd_rate = firstFixingRate
                    firstPayment = True

                else:

                    dfStart = index_curve.df(startAccruedDt)
                    dfEnd = index_curve.df(endAccruedDt)
                    fwd_rate = (dfStart / dfEnd - 1.0) / index_alpha

                pmntAmount = (fwd_rate + self._spread) * pay_alpha * notional

                dfPmnt = discount_curve.df(pmntDate) / dfValue
                pmntPV = pmntAmount * dfPmnt
                legPV += pmntPV

                self._rates.append(fwd_rate)
                self._payments.append(pmntAmount)
                self._paymentDfs.append(dfPmnt)
                self._paymentPVs.append(pmntPV)
                self._cumulativePVs.append(legPV)

            else:

                self._rates.append(0.0)
                self._payments.append(0.0)
                self._paymentDfs.append(0.0)
                self._paymentPVs.append(0.0)
                self._cumulativePVs.append(legPV)

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
        print("SPREAD (bp):", self._spread * 10000)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))

        if len(self._payment_dates) == 0:
            print("Payments Dates not calculated.")
            return

        header = "PAY_DATE     ACCR_START   ACCR_END      DAYS  YEARFRAC"
        print(header)

        num_flows = len(self._payment_dates)

        for iFlow in range(0, num_flows):
            print("%11s  %11s  %11s  %4d  %8.6f  " %
                  (self._payment_dates[iFlow],
                   self._startAccruedDates[iFlow],
                   self._endAccruedDates[iFlow],
                   self._accrued_days[iFlow],
                   self._year_fracs[iFlow]))

###############################################################################

    def print_valuation(self):
        """ Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. """

        print("START DATE:", self._effective_date)
        print("MATURITY DATE:", self._maturity_date)
        print("SPREAD (BPS):", self._spread * 10000)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))

        if len(self._payments) == 0:
            print("Payments not calculated.")
            return

        header = "PAY_DATE     ACCR_START   ACCR_END     DAYS  YEARFRAC"
        header += "    IBOR      PAYMENT       DF          PV        CUM PV"
        print(header)

        num_flows = len(self._payment_dates)

        for iFlow in range(0, num_flows):
            print("%11s  %11s  %11s  %4d  %8.6f  %9.5f  % 11.2f  %10.8f  % 11.2f  % 11.2f" %
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

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._effective_date)
        s += label_to_string("TERMINATION DATE", self._termination_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("SWAP TYPE", self._leg_type)
        s += label_to_string("SPREAD (BPS)", self._spread*10000)
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
