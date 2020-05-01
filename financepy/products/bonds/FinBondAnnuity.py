# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""

from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes

###############################################################################


class FinBondAnnuity(object):
    ''' An annuity is a vector of dates and flows generated according to ISDA 
    standard rules which starts on the next date after the start date 
    (effective date) and runs up to an end date. Dates are then adjusted 
    according to a specified calendar. '''

    def __init__(self,
                 startDate,
                 endDate,
                 frequencyType=FinFrequencyTypes.ANNUAL,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD,
                 dayCountConventionType=FinDayCountTypes.ACT_360):

        if startDate > endDate:
            raise ValueError("Start Date after End Date")

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinDayAdjustTypes:
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
        self._flows = []

        # now generate the schedule
        self.generate(startDate)

###############################################################################

    def generate(self, startDate):

        self._flows = []
        self._yearFractions = []

        dayCount = FinDayCount(self._dayCountConventionType)

        numFlows = len(self._schedule._adjustedDates)

        if numFlows > 0:
            yearFrac = dayCount.yearFrac(
                startDate, self._schedule._adjustedDates[0])
            self._yearFractions.append(yearFrac)

        for i in range(1, numFlows):
            prevDate = self._schedule._adjustedDates[i - 1]
            nextDate = self._schedule._adjustedDates[i]

            yearFrac = dayCount.yearFrac(prevDate, nextDate)
            self._yearFractions.append(yearFrac)

        return self._yearFractions

###############################################################################

    def dump(self):
        print("START DATE:", self._startDate)
        print("END DATE:", self._endDate)
        print("FREQUENCY:", self._frequencyType)
        print("CALENDAR:", self._calendarType)
        print("BUSDAYRULE:", self._busDayAdjustType)
        print("DATEGENRULE:", self._dateGenRuleType)

        numFlows = len(self._schedule._adjustedDates)
        for i in range(0, numFlows):
            print(self._schedule._adjustedDates[i], self._yearFractions[i])

###############################################################################
