##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

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
                 start_dt: Date,
                 end_dt: Date,
                 principal: float,
                 freq_type: FrequencyTypes = FrequencyTypes.MONTHLY,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
                 dc_type: DayCountTypes = DayCountTypes.ACT_360):
        """ Create the mortgage using start and end dates and principal. """

        check_argument_types(self.__init__, locals())

        if start_dt > end_dt:
            raise FinError("Start Date after End Date")

        self._start_dt = start_dt
        self._end_dt = end_dt
        self._principal = principal
        self._freq_type = freq_type
        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type
        self._dc_type = dc_type

        self._schedule = Schedule(start_dt,
                                  end_dt,
                                  self._freq_type,
                                  self._cal_type,
                                  self._bd_type,
                                  self._dg_type)

###############################################################################

    def repayment_amount(self,
                         zero_rate: float):
        """ Determine monthly repayment amount based on current zero rate. """

        frequency = annual_frequency(self._freq_type)

        num_flows = len(self._schedule._adjusted_dts)
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

        num_flows = len(self._schedule._adjusted_dts)
        principal = self._principal
        frequency = annual_frequency(self._freq_type)

        if mortgage_type == BondMortgageTypes.REPAYMENT:
            monthly_flow = self.repayment_amount(zero_rate)
        elif mortgage_type == BondMortgageTypes.INTEREST_ONLY:
            monthly_flow = zero_rate * self._principal / frequency
        else:
            raise FinError("Unknown Mortgage type.")

        for _ in range(1, num_flows):
            interest_flow = principal * zero_rate / frequency
            principal_flow = monthly_flow - interest_flow
            principal = principal - principal_flow
            self._interest_flows.append(interest_flow)
            self._principal_flows.append(principal_flow)
            self._principal_remaining.append(principal)
            self._total_flows.append(monthly_flow)

###############################################################################

    def print_leg(self):
        print("START DATE:", self._start_dt)
        print("MATURITY DATE:", self._end_dt)
        print("MORTGAGE TYPE:", self._mortgage_type)
        print("FREQUENCY:", self._freq_type)
        print("CALENDAR:", self._cal_type)
        print("BUSDAYRULE:", self._bd_type)
        print("DATEGENRULE:", self._dg_type)

        num_flows = len(self._schedule._adjusted_dts)

        print("%15s %12s %12s %12s %12s" %
              ("PAYMENT DATE", "INTEREST", "PRINCIPAL",
               "OUTSTANDING", "TOTAL"))

        print("")
        for i in range(0, num_flows):
            print("%15s %12.2f %12.2f %12.2f %12.2f" %
                  (self._schedule._adjusted_dts[i],
                   self._interest_flows[i],
                   self._principal_flows[i],
                   self._principal_remaining[i],
                   self._total_flows[i]))

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._start_dt)
        s += label_to_string("MATURITY DATE", self._end_dt)
        s += label_to_string("MORTGAGE TYPE", self._mortgage_type)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUSDAYRULE", self._bd_type)
        s += label_to_string("DATEGENRULE", self._dg_type)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
