##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.schedule import Schedule
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.date import Date
from ...utils.helpers import label_to_string, check_argument_types

###############################################################################

from enum import Enum


class BondMortgageTypes(Enum):
    REPAYMENT = 1
    INTEREST_ONLY = 2

###############################################################################


class BondMortgage:
    """ A mortgage is a vector of dates and flows generated in order to repay
    a fixed amount given a known interest rate. Payments are all the same
    amount but with a varying mixture of interest and repayment of principal.
    """

    def __init__(self,
                 start_date: Date,
                 end_date: Date,
                 principal: float,
                 freq_type: FrequencyTypes = FrequencyTypes.MONTHLY,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 day_count_convention_type: DayCountTypes = DayCountTypes.ACT_360):
        """ Create the mortgage using start and end dates and principal. """

        check_argument_types(self.__init__, locals())

        if start_date > end_date:
            raise FinError("Start Date after End Date")

        self._start_date = start_date
        self._end_date = end_date
        self._principal = principal
        self._freq_type = freq_type
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type
        self._day_count_convention_type = day_count_convention_type

        self._schedule = Schedule(start_date,
                                  end_date,
                                  self._freq_type,
                                  self._calendar_type,
                                  self._bus_day_adjust_type,
                                  self._date_gen_rule_type)

###############################################################################

    def repayment_amount(self,
                         zero_rate: float):
        """ Determine monthly repayment amount based on current zero rate. """

        frequency = annual_frequency(self._freq_type)

        num_flows = len(self._schedule._adjusted_dates)
        p = (1.0 + zero_rate/frequency) ** (num_flows-1)
        m = zero_rate * p / (p - 1.0) / frequency
        m = m * self._principal
        return m

###############################################################################

    def generate_flows(self,
                       zero_rate: float,
                       mortgage_type: BondMortgageTypes):
        """ Generate the bond flow amounts. """

        self._mortgage_type = mortgage_type
        self._interest_flows = [0]
        self._principal_flows = [0]
        self._principal_remaining = [self._principal]
        self._total_flows = [0]

        num_flows = len(self._schedule._adjusted_dates)
        principal = self._principal
        frequency = annual_frequency(self._freq_type)

        if mortgage_type == BondMortgageTypes.REPAYMENT:
            monthly_flow = self.repayment_amount(zero_rate)
        elif mortgage_type == BondMortgageTypes.INTEREST_ONLY:
            monthly_flow = zero_rate * self._principal / frequency
        else:
            raise FinError("Unknown Mortgage type.")

        for i in range(1, num_flows):
            interestFlow = principal * zero_rate / frequency
            principalFlow = monthly_flow - interestFlow
            principal = principal - principalFlow
            self._interest_flows.append(interestFlow)
            self._principal_flows.append(principalFlow)
            self._principal_remaining.append(principal)
            self._total_flows.append(monthly_flow)

###############################################################################

    def print_leg(self):
        print("START DATE:", self._start_date)
        print("MATURITY DATE:", self._end_date)
        print("MORTGAGE TYPE:", self._mortgage_type)
        print("FREQUENCY:", self._freq_type)
        print("CALENDAR:", self._calendar_type)
        print("BUSDAYRULE:", self._bus_day_adjust_type)
        print("DATEGENRULE:", self._date_gen_rule_type)

        num_flows = len(self._schedule._adjusted_dates)

        print("%15s %12s %12s %12s %12s" %
              ("PAYMENT DATE", "INTEREST", "PRINCIPAL",
               "OUTSTANDING", "TOTAL"))

        print("")
        for i in range(0, num_flows):
            print("%15s %12.2f %12.2f %12.2f %12.2f" %
                  (self._schedule._adjusted_dates[i],
                   self._interest_flows[i],
                   self._principal_flows[i],
                   self._principal_remaining[i],
                   self._total_flows[i]))

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._start_date)
        s += label_to_string("MATURITY DATE", self._end_date)
        s += label_to_string("MORTGAGE TYPE", self._mortgage_type)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUSDAYRULE", self._bus_day_adjust_type)
        s += label_to_string("DATEGENRULE", self._date_gen_rule_type)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
