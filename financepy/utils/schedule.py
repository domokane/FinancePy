##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from .error import FinError
from .date import Date
from .calendar import (Calendar, CalendarTypes)
from .calendar import (BusDayAdjustTypes, DateGenRuleTypes)
from .frequency import (annual_frequency, FrequencyTypes)
from .helpers import label_to_string
from .helpers import check_argument_types


###############################################################################
# TODO: Start and end date to allow for long stubs
###############################################################################


class Schedule:
    """ A schedule is a set of dates generated according to ISDA standard
    rules which starts on the next date after the effective date and runs up to
    a termination date. Dates are adjusted to a provided calendar. The zeroth
    element is the previous coupon date (PCD) and the first element is the
    Next Coupon Date (NCD). We reference ISDA 2006."""

    def __init__(self,
                 effective_date: Date,  # Also known as the start date
                 # This is UNADJUSTED (set flag to adjust it)
                 termination_date: Date,
                 freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 adjust_termination_date: bool = True,  # Default is to adjust
                 end_of_month: bool = False,  # All flow dates are EOM if True
                 first_date=None,  # First coupon date
                 next_to_last_date=None):  # Penultimate coupon date
        """ Create Schedule object which calculates a sequence of dates
        following the ISDA convention for fixed income products, mainly swaps. 

        If the date gen rule type is FORWARD we get the unadjusted dates by stepping 
        forward from the effective date in steps of months determined by the period
        tenor - i.e. the number of months between payments. We stop before we go past the
        termination date. 

        If the date gen rule type is BACKWARD we get the unadjusted dates by 
        stepping backward from the termination date in steps of months determined by 
        the period tenor - i.e. the number of months between payments. We stop 
        before we go past the effective date. 

        - If the EOM flag is false, and the start date is on the 31st then the
        the unadjusted dates will fall on the 30 if a 30 is a previous date. 
        - If the EOM flag is false and the start date is 28 Feb then all 
        unadjusted dates will fall on the 28th. 
        - If the EOM flag is false and the start date is 28 Feb then all 
        unadjusted dates will fall on their respective EOM.

        We then adjust all of the flow dates if they fall on a weekend or holiday 
        according to the calendar specified. These dates are adjusted in 
        accordance with the business date adjustment.

        The effective date is never adjusted as it is not a payment date.
        The termination date is not automatically business day adjusted in a 
        swap - assuming it is a holiday date. This must be explicitly stated in
        the trade confirm. However, it is adjusted in a CDS contract as standard. 

        Inputs first_date and next_to_last_date are for managing long payment stubs
        at the start and end of the swap but *have not yet been implemented*. All
        stubs are currently short, either at the start or end of swap. """

        check_argument_types(self.__init__, locals())

        if effective_date >= termination_date:
            raise FinError("Effective date must be before termination date.")

        self._effective_date = effective_date
        self._termination_date = termination_date

        if first_date is None:
            self._first_date = effective_date
        else:
            if first_date > effective_date and first_date < termination_date:
                self._first_date = first_date
                print("FIRST DATE NOT IMPLEMENTED")  # TODO
            else:
                raise FinError("First date must be after effective date and" +
                               " before termination date")

        if next_to_last_date is None:
            self._next_to_last_date = termination_date
        else:
            if next_to_last_date > effective_date and next_to_last_date < termination_date:
                self._next_to_last_date = next_to_last_date
                print("NEXT TO LAST DATE NOT IMPLEMENTED")  # TODO
            else:
                raise FinError("Next to last date must be after effective date and" +
                               " before termination date")

        self._freq_type = freq_type
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._adjust_termination_date = adjust_termination_date

        if end_of_month is True:
            self._end_of_month = True
        else:
            self._end_of_month = False

        self._adjusted_dates = None

        self._generate()

    ###############################################################################

    def schedule_dates(self):
        """ Returns a list of the schedule of Dates. """

        if self._adjusted_dates is None:
            self._generate()

        return self._adjusted_dates

    ###############################################################################

    def _generate(self):
        """ Generate schedule of dates according to specified date generation
        rules and also adjust these dates for holidays according to the
        specified business day convention and the specified calendar. """

        calendar = Calendar(self._calendar_type)
        frequency = annual_frequency(self._freq_type)
        num_months = int(12 / frequency)

        unadjusted_schedule_dates = []
        self._adjusted_dates = []

        if self._date_gen_rule_type == DateGenRuleTypes.BACKWARD:

            next_date = self._termination_date
            flow_num = 0

            while next_date > self._effective_date:

                unadjusted_schedule_dates.append(next_date)
                tot_num_months = num_months * (1 + flow_num)
                next_date = self._termination_date.add_months(-tot_num_months)

                if self._end_of_month is True:
                    next_date = next_date.eom()

                flow_num += 1

            # Add on the Previous Coupon Date
            unadjusted_schedule_dates.append(next_date)
            flow_num += 1

            # reverse order and holiday adjust dates
            # the first date is not adjusted as this was provided
            dt = unadjusted_schedule_dates[flow_num - 1]
            self._adjusted_dates.append(dt)

            # We adjust all flows after the effective date and before the
            # termination date to fall on business days according to their cal
            for i in range(1, flow_num - 1):
                dt = calendar.adjust(unadjusted_schedule_dates[flow_num - i - 1],
                                     self._bus_day_adjust_type)

                self._adjusted_dates.append(dt)

            self._adjusted_dates.append(self._termination_date)

        elif self._date_gen_rule_type == DateGenRuleTypes.FORWARD:

            # This needs checking
            next_date = self._effective_date
            flow_num = 0

            unadjusted_schedule_dates.append(next_date)
            flow_num = 1

            while next_date < self._termination_date:
                unadjusted_schedule_dates.append(next_date)
                tot_num_months = num_months * (flow_num)
                next_date = self._effective_date.add_months(tot_num_months)
                flow_num = flow_num + 1

            # The effective date is not adjusted as it is given
            for i in range(1, flow_num):
                dt = calendar.adjust(unadjusted_schedule_dates[i],
                                     self._bus_day_adjust_type)

                self._adjusted_dates.append(dt)

            self._adjusted_dates.append(self._termination_date)

        if self._adjusted_dates[0] < self._effective_date:
            self._adjusted_dates[0] = self._effective_date

        # The market standard for swaps is not to adjust the termination date
        # unless it is specified in the contract. It is standard for CDS.
        # We change it if the adjust_termination_date flag is True.
        if self._adjust_termination_date is True:
            self._termination_date = calendar.adjust(self._termination_date,
                                                     self._bus_day_adjust_type)

            self._adjusted_dates[-1] = self._termination_date

        #######################################################################
        # Check the resulting schedule to ensure that no two dates are the
        # same and that they are monotonic - this should never happen but ...
        #######################################################################

        if len(self._adjusted_dates) < 2:
            raise FinError("Schedule has two dates only.")

        prev_dt = self._adjusted_dates[0]
        for dt in self._adjusted_dates[1:]:

            if dt == prev_dt:
                raise FinError("Two matching dates in schedule")

            if dt < prev_dt:  # Dates must be ordered
                raise FinError("Dates are not monotonic")

            prev_dt = dt

        #######################################################################

        return self._adjusted_dates

    ##############################################################################

    def __repr__(self):
        """ Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self._effective_date)
        s += label_to_string("END DATE", self._termination_date)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUSDAYRULE", self._bus_day_adjust_type)
        s += label_to_string("DATEGENRULE", self._date_gen_rule_type)
        s += label_to_string("ADJUST TERM DATE", self._adjust_termination_date)
        s += label_to_string("END OF MONTH", self._end_of_month, "")

        if 1 == 0:
            if len(self._adjusted_dates) > 0:
                s += "\n\n"
                s += label_to_string("EFF", self._adjusted_dates[0], "")

            if len(self._adjusted_dates) > 1:
                s += "\n"
                s += label_to_string("FLW", self._adjusted_dates[1:], "",
                                     listFormat=True)

        return s

    ###############################################################################

    def _print(self):
        """ Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations. """
        print(self)

###############################################################################
