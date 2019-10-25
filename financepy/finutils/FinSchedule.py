# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""

from . import FinError
from .FinCalendar import (FinCalendar,FinCalendarTypes)
from .FinCalendar import (FinBusDayConventionTypes, FinDateGenRuleTypes)
from .FinFrequency import (FinFrequency, FinFrequencyTypes)

###############################################################################

class FinSchedule(object):

    ''' A Schedule is a vector of dates generated according to ISDA standard 
    rules which starts on the next date after the start date and runs up to 
    an end date. Dates are adjusted to a provided calendar. The zeroth 
    element is the PCD and the first element is the NCD '''

    def __init__(self,
                 startDate,
                 endDate,
                 frequencyType = FinFrequencyTypes.ANNUAL,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayConventionTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):

        if startDate > endDate:
            raise FinError("Start Date after End Date")

        if calendarType not in FinCalendarTypes:
            raise FinError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise FinError("Unknown Business Day Adjust type " + str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise FinError("Unknown Date Gen Rule type " + str(dateGenRuleType))

        # validation complete
        self._startDate = startDate
        self._endDate = endDate
        self._frequencyType = frequencyType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType
        self._adjustedDates = []

        self.generate()

###############################################################################

    def generate(self):
        ''' Generate schedule of dates according to specified date generation 
        rules and also adjust these dates for holidays according to 
        the business day convention and the specified calendar. '''

        self._adjustedDates = []
        calendar = FinCalendar(self._calendarType)
        frequency = FinFrequency(self._frequencyType)
        numMonths = int(12/frequency)

        unadjustedScheduleDates = []

        if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

            nextDate = self._endDate
            flowNum = 0

            while nextDate > self._startDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(-numMonths)
                flowNum += 1

            # Add on the Previous Coupon Date
            unadjustedScheduleDates.append(nextDate)
            flowNum += 1

            # reverse order and holiday adjust dates
            for i in range(0, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[flowNum-i-1], 
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

        elif self._dateGenRuleType == FinDateGenRuleTypes.FORWARD:

            # This needs checking
            nextDate = self._startDate
            flowNum = 0

            unadjustedScheduleDates.append(nextDate)
            flowNum = 1

            while nextDate < self._endDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(numMonths)
                flowNum = flowNum + 1

            for i in range(1, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[i], 
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

            self._adjustedDates.append(self._endDate)

        return self._adjustedDates

###############################################################################

    def print(self):
        ''' Print out the details of the schedule and the actual dates. This 
        can be used for providing transparency on schedule calculations. '''
        print("START DATE:", self._startDate)
        print("END DATE:", self._endDate)
        print("FREQUENCY:",self._frequencyType)
        print("CALENDAR:",self._calendarType)
        print("BUSDAYRULE:",self._busDayAdjustType)
        print("DATEGENRULE:",self._dateGenRuleType)

        s = ""
        for dt in self._adjustedDates:
            print(str(dt))
            s = s + "\n" + str(dt)

        return s

###############################################################################
