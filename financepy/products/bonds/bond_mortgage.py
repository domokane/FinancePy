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
    """A mortgage is a vector of dates and flows generated in order to repay
    a fixed amount given a known interest rate. Payments are all the same
    amount but with a varying mixture of interest and repayment of principal.
    """

    def __init__(
        self,
        start_dt: Date,
        end_dt: Date,
        principal: float,
        freq_type: FrequencyTypes = FrequencyTypes.MONTHLY,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
    ):
        """Create the mortgage using start and end dates and principal."""

        check_argument_types(self.__init__, locals())

        if start_dt > end_dt:
            raise FinError("Start Date after End Date")

        self.start_dt = start_dt
        self.end_dt = end_dt
        self.principal = principal
        self.freq_type = freq_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type
        self.dc_type = dc_type

        self.schedule = Schedule(
            start_dt,
            end_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

        self.mortgage_type = None
        self.interest_flows = None
        self.principal_flows = None
        self.principal_remaining = None
        self.total_flows = None

    ###########################################################################

    def repayment_amount(self, zero_rate: float):
        """Determine monthly repayment amount based on current zero rate."""

        frequency = annual_frequency(self.freq_type)

        num_flows = len(self.schedule.adjusted_dts)
        p = (1.0 + zero_rate / frequency) ** (num_flows - 1)
        m = zero_rate * p / (p - 1.0) / frequency
        m = m * self.principal
        return m

    ###########################################################################

    def generate_flows(
        self, zero_rate: float, mortgage_type: BondMortgageTypes
    ):
        """Generate the bond flow amounts."""

        self.mortgage_type = mortgage_type
        self.interest_flows = [0]
        self.principal_flows = [0]
        self.principal_remaining = [self.principal]
        self.total_flows = [0]

        num_flows = len(self.schedule.adjusted_dts)
        principal = self.principal
        frequency = annual_frequency(self.freq_type)

        if mortgage_type == BondMortgageTypes.REPAYMENT:
            monthly_flow = self.repayment_amount(zero_rate)
        elif mortgage_type == BondMortgageTypes.INTEREST_ONLY:
            monthly_flow = zero_rate * self.principal / frequency
        else:
            raise FinError("Unknown Mortgage type.")

        for _ in range(1, num_flows):
            interest_flow = principal * zero_rate / frequency
            principal_flow = monthly_flow - interest_flow
            principal = principal - principal_flow
            self.interest_flows.append(interest_flow)
            self.principal_flows.append(principal_flow)
            self.principal_remaining.append(principal)
            self.total_flows.append(monthly_flow)

    ###########################################################################

    def print_leg(self):
        print("START DATE:", self.start_dt)
        print("MATURITY DATE:", self.end_dt)
        print("MORTGAGE TYPE:", self.mortgage_type)
        print("FREQUENCY:", self.freq_type)
        print("CALENDAR:", self.cal_type)
        print("BUSDAYRULE:", self.bd_type)
        print("DATEGENRULE:", self.dg_type)

        num_flows = len(self.schedule.adjusted_dts)

        print(
            "%15s %12s %12s %12s %12s"
            % ("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING", "TOTAL")
        )

        print("")
        for i in range(0, num_flows):
            print(
                "%15s %12.2f %12.2f %12.2f %12.2f"
                % (
                    self.schedule.adjusted_dts[i],
                    self.interest_flows[i],
                    self.principal_flows[i],
                    self.principal_remaining[i],
                    self.total_flows[i],
                )
            )

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.start_dt)
        s += label_to_string("MATURITY DATE", self.end_dt)
        s += label_to_string("MORTGAGE TYPE", self.mortgage_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUSDAYRULE", self.bd_type)
        s += label_to_string("DATEGENRULE", self.dg_type)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
