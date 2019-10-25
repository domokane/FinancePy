# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

#from .FinMath import FinMath
from .FinDate import FinDate
from .FinError import FinError

easterMondayDay = [98, 90, 103, 95, 114, 106, 91, 111, 102, 87,
                   107, 99, 83, 103, 95, 115, 99, 91, 111, 96, 87,
                   107, 92, 112, 103, 95, 108, 100, 91,
                   111, 96, 88, 107, 92, 112, 104, 88, 108, 100,
                   85, 104, 96, 116, 101, 92, 112, 97, 89, 108,
                   100, 85, 105, 96, 109, 101, 93, 112, 97, 89,
                   109, 93, 113, 105, 90, 109, 101, 86, 106, 97,
                   89, 102, 94, 113, 105, 90, 110, 101, 86, 106,
                   98, 110, 102, 94, 114, 98, 90, 110, 95, 86,
                   106, 91, 111, 102, 94, 107, 99, 90, 103, 95,
                   115, 106, 91, 111, 103, 87, 107, 99, 84, 103,
                   95, 115, 100, 91, 111, 96, 88, 107, 92, 112,
                   104, 95, 108, 100, 92, 111, 96, 88, 108, 92,
                   112, 104, 89, 108, 100, 85, 105, 96, 116, 101,
                   93, 112, 97, 89, 109, 100, 85, 105, 97, 109,
                   101, 93, 113, 97, 89, 109, 94, 113, 105, 90,
                   110, 101, 86, 106, 98, 89, 102, 94, 114, 105,
                   90, 110, 102, 86, 106, 98, 111, 102, 94, 114,
                   99, 90, 110, 95, 87, 106, 91, 111, 103, 94,
                   107, 99, 91, 103, 95, 115, 107, 91, 111, 103,
                   88, 108, 100, 85, 105, 96, 109, 101, 93, 112,
                   97, 89, 109, 93, 113, 105, 90, 109, 101, 86,
                   106, 97, 89, 102, 94, 113, 105, 90, 110, 101,
                   86, 106, 98, 110, 102, 94, 114, 98, 90, 110,
                   95, 86, 106, 91, 111, 102, 94, 107, 99, 90,
                   103, 95, 115, 106, 91, 111, 103, 87, 107, 99,
                   84, 103, 95, 115, 100, 91, 111, 96, 88, 107,
                   92, 112, 104, 95, 108, 100, 92, 111, 96, 88,
                   108, 92, 112, 104, 89, 108, 100, 85, 105, 96,
                   116, 101, 93, 112, 97, 89, 109, 100, 85, 105]

from enum import Enum

class FinBusDayConventionTypes(Enum):
    NONE=1
    FOLLOWING=2
    MODIFIED_FOLLOWING=3
    PRECEDING=4
    MODIFIED_PRECEDING=5
    
class FinCalendarTypes(Enum):
    TARGET=1
    US=2
    UK=3
    NONE=4
    WEEKEND=5
    
class FinDateGenRuleTypes(Enum):
    FORWARD=1
    BACKWARD=2
    
################################################################################
    
class FinCalendar(object):

    ''' Class to manage designation of payment dates as holidays according to a 
    calendar convention specified by the user. '''

    def __init__(self, calendarType):
        ''' Create a calendar based on a specified calendar type. '''

        if calendarType not in FinCalendarTypes:
            raise FinError("Need to pass FinCalendarType and not " + str(calendarType))

        self._type = calendarType

    ###########################################################################
    
    def adjust(self, dt, busDayConventionType ):
        ''' Adjust a payment date if it falls on a holiday according to the 
        specified business day convention. '''

        m = dt._m

        if busDayConventionType == FinBusDayConventionTypes.NONE:
        
            return dt
        
        elif busDayConventionType == FinBusDayConventionTypes.FOLLOWING:

            # step forward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(1)

            return dt

        elif busDayConventionType == FinBusDayConventionTypes.MODIFIED_FOLLOWING:

            # step forward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(1)

            # if the business day is in a different month look back 
            # for previous first business day one day at a time
            # I could speed this up by starting it at initial date
            if dt._m != m:
                while self.isBusinessDay(dt) is False:
                    dt.addDays(-1)

            return dt

        elif busDayConventionType == FinBusDayConventionTypes.PRECEDING:

            # if the business day is in the next month look back 
            # for previous first business day one day at a time
            while self.isBusinessDay(dt) is False:
                dt.addDays(-1)

            return dt

        elif busDayConventionType == FinBusDayConventionTypes.MODIFIED_PRECEDING:

            # step backward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt.addDays(-1)

            # if the business day is in a different month look forward 
            # for previous first business day one day at a time
            # I could speed this up by starting it at initial date
            if dt._m != m:
                while self.isBusinessDay(dt) is False:
                    dt.addDays(+1)

            return dt
        else:
            raise FinError("Unknown adjustment convention",str(busDayConventionType))

        return dt

###############################################################################

    def isBusinessDay(self, dt):
        ''' Determines if a date is a business day according to the calendar. '''

        y = dt._y
        m = dt._m
        d = dt._d

        startDate = FinDate(y, 1, 1)

        dd = dt._excelDate - startDate._excelDate

        weekday = dt._weekday

        em = easterMondayDay[y-1901]

        if dt.isWeekend() is True:
            return False

        if self._type == FinCalendarTypes.NONE:
            return True
        
        if self._type == FinCalendarTypes.WEEKEND:
            return True

        if self._type == FinCalendarTypes.UK:

            if m == 1 and d == 1:  # new years day
                return False

            if dd == em:
                return False

            if dd == em - 3:  # good friday
                return False

            if m == 5 and d <= 7 and weekday == FinDate.MON:
                return False

            if m == 5 and d >= 25 and weekday == FinDate.MON:
                return False

            if m == 8 and d <= 7 and weekday == FinDate.MON:
                return False

            if m == 12 and d == 25:  # Xmas
                return False

            if m == 12 and d == 26:  # Boxing day
                return False

            if m == 12 and d == 27 and weekday == FinDate.MON:  # Xmas
                return False

            if m == 12 and d == 27 and weekday == FinDate.TUE:  # Xmas
                return False

            if m == 12 and d == 28 and weekday == FinDate.MON:  # Xmas
                return False

            if m == 12 and d == 28 and weekday == FinDate.TUE:  # Xmas
                return False

            return True

        elif self._type == FinCalendarTypes.US:

            if m == 1 and d == 1:  # NYD
                return False

            if m == 12 and d == 31 and weekday == FinDate.FRI:  # NYE
                return False

            if dd == em:
                return False

            if dd == em - 3:
                return False

            if m == 1 and d >= 15 and weekday == FinDate.MON:
                return False

            if m == 2 and d >= 15 and d <= 21 and weekday == FinDate.MON:
                return False

            if m == 5 and d >= 25 and d <= 21 and weekday == FinDate.MON:
                return False

            if m == 7 and d == 4:  # Indep day
                return False

            if m == 7 and d == 5 and weekday == FinDate.MON:  # Indep day
                return False

            if m == 7 and d == 3 and weekday == FinDate.FRI:  # Indep day
                return False

            if m == 9 and d >= 8 and d <= 14 and weekday == FinDate.MON:
                return False

            if m == 11 and d == 11:  # Veterans day
                return False

            if m == 11 and d == 12 and weekday == FinDate.MON:  # Columbus day
                return False

            if m == 11 and d == 10 and weekday == FinDate.FRI:  # Columbus day
                return False

            if m == 11 and d >= 22 and d <= 28 and weekday == FinDate.THU:
                return False

            if m == 12 and d == 25:  # Xmas holiday
                return False

            if m == 12 and d == 26 and weekday == FinDate.MON:
                return False

            if m == 12 and d == 24 and weekday == FinDate.FRI:
                return False

            return True

        elif self._type == FinCalendarTypes.TARGET:

            if m == 1 and d == 1:  # new year's day
                return False

            if m == 5 and d == 1:  # May day
                return False

            if dd == em:  # Easter monday holiday
                return False

            if m == 12 and d == 25:  # Xmas bank holiday
                return False

            if m == 12 and d == 26:  # Xmas bank holiday
                return False

            if m == 12 and d == 31:  # NYD bank holiday
                return False

            return True

###############################################################################

    def easterMonday(self, y):
        ''' Get the day in a givenm year that is Easter Monday. This is not easy 
        to compute so we rely on a pre-calculated array. '''

        if y > 2100:
            raise FinError("Unable to determine Easter monday in year " + str(y))

        easterMonday = easterMondayDay[y - 1901]
        return easterMonday

###############################################################################

    def str(self):
        s = self._type
        return s

###############################################################################
