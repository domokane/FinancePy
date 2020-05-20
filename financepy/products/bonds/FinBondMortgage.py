##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString

###############################################################################

from enum import Enum


class FinBondMortgageType(Enum):
    REPAYMENT = 1
    INTEREST_ONLY = 2

###############################################################################


class FinBondMortgage(object):
    ''' A mortgage is a vector of dates and flows generated in order to repay
    a fixed amount given a known interest rate. Payments are all the same
    amount but with a varying mixture of interest and repayment of principal.
    '''
###############################################################################

    def __init__(self,
                 startDate,
                 endDate,
                 principal,
                 frequencyType=FinFrequencyTypes.MONTHLY,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD,
                 dayCountConventionType=FinDayCountTypes.ACT_360):
        ''' Create the mortgage using start and end dates and principal. '''

        if startDate > endDate:
            raise ValueError("Start Date after End Date")

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayAdjustTypes:
            raise ValueError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise ValueError(
                "Unknown Date Gen Rule type " +
                str(dateGenRuleType))

        if dayCountConventionType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Day Count type " +
                str(dayCountConventionType))

        self._startDate = startDate
        self._endDate = endDate
        self._principal = principal
        self._frequencyType = frequencyType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType
        self._dayCountConventionType = dayCountConventionType

        self._schedule = FinSchedule(startDate,
                                     endDate,
                                     self._frequencyType,
                                     self._calendarType,
                                     self._busDayAdjustType,
                                     self._dateGenRuleType)

###############################################################################

    def repaymentAmount(self, zeroRate):
        ''' Determine monthly repayment amount based on current zero rate. '''

        frequency = FinFrequency(self._frequencyType)

        numFlows = len(self._schedule._adjustedDates)
        p = (1.0 + zeroRate/frequency) ** (numFlows-1)
        m = zeroRate * p / (p - 1.0) / frequency
        m = m * self._principal
        return m

###############################################################################

    def generateFlows(self, zeroRate, mortgageType):
        ''' Generate the bond flow amounts. '''

        self._mortgageType = mortgageType
        self._interestFlows = [0]
        self._principalFlows = [0]
        self._principalRemaining = [self._principal]
        self._totalFlows = [0]

        numFlows = len(self._schedule._adjustedDates)
        principal = self._principal
        frequency = FinFrequency(self._frequencyType)

        if mortgageType == FinBondMortgageType.REPAYMENT:
            monthlyFlow = self.repaymentAmount(zeroRate)
        elif mortgageType == FinBondMortgageType.INTEREST_ONLY:
            monthlyFlow = zeroRate * self._principal / frequency
        else:
            raise FinError("Unknown Mortgage type.")

        for i in range(1, numFlows):
            interestFlow = principal * zeroRate / frequency
            principalFlow = monthlyFlow - interestFlow
            principal = principal - principalFlow
            self._interestFlows.append(interestFlow)
            self._principalFlows.append(principalFlow)
            self._principalRemaining.append(principal)
            self._totalFlows.append(monthlyFlow)

###############################################################################

    def __repr__(self):
        s = labelToString("START DATE:", self._startDate)
        s += labelToString("END DATE:", self._endDate)
        s += labelToString("MORTGAGE TYPE:", self._mortgageType)
        s += labelToString("FREQUENCY:", self._frequencyType)
        s += labelToString("CALENDAR:", self._calendarType)
        s += labelToString("BUSDAYRULE:", self._busDayAdjustType)
        s += labelToString("DATEGENRULE:", self._dateGenRuleType)

        numFlows = len(self._schedule._adjustedDates)

        print("%15s %12s %12s %12s %12s" %
              ("PAYMENT DATE", "INTEREST", "PRINCIPAL",
               "OUTSTANDING", "TOTAL"))

        print("")
        for i in range(0, numFlows):
            print("%15s %12.2f %12.2f %12.2f %12.2f" %
                  (self._schedule._adjustedDates[i],
                   self._interestFlows[i],
                   self._principalFlows[i],
                   self._principalRemaining[i],
                   self._totalFlows[i]))

        return s

###############################################################################

    def print(self):
        print(self)

###############################################################################
