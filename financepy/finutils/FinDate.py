# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

# TODO Vectorize the functions such as addDays and addMonths

import datetime
from .FinError import FinError
from .FinMath import isLeapYear

shortDayNames = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
longDayNames = [
    'MONDAY',
    'TUESDAY',
    'WEDNESDAY',
    'THURSDAY',
    'FRIDAY',
    'SATURDAY',
    'SUNDAY']
shortMonthNames = [
    'JAN',
    'FEB',
    'MAR',
    'APR',
    'MAY',
    'JUN',
    'JUL',
    'AUG',
    'SEP',
    'OCT',
    'NOV',
    'DEC']
longMonthNames = [
    'JANUARY',
    'FEBRUARY',
    'MARCH',
    'APRIL',
    'MAY',
    'JUNE',
    'JULY',
    'AUGUST',
    'SEPTEMBER',
    'OCTOBER',
    'NOVEMBER',
    'DECEMBER']

monthDaysNotLeapYear = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
monthDaysLeapYear = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

###############################################################################
# Date functions that are not class members but are useful
###############################################################################


def dailyWorkingDaySchedule(self, startDate, endDate):
    dateList = []

    dt = startDate
    dateList.append(dt)
    while dt < endDate:
        dt = dt.addWorkDays(1)
        dateList.append(dt)

    return dateList

###############################################################################

def datediff(d1, d2):
    ''' Calculate the number of days between two dates. '''

    return (d2._excelDate - d1._excelDate)

###############################################################################

class FinDate():

    ''' Date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance. '''

    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6

    ###########################################################################

    def __init__(self, y_or_d, m, d_or_y):
        ''' Create a date given year, month and day of month. The order is not
        enforced so 4th July 2019 can be created as FinDate(4,7,2019) or as
        FinDate(2019,7,4) so long as the middle number is the month.
        The year must be a 4-digit number greater than or equal to 1900. '''

        if y_or_d >= 1900:
            y = y_or_d
            d = d_or_y
        elif d_or_y >= 1900:
            y = d_or_y
            d = y_or_d
        else:
            raise FinError("Invalid date: Inconsistent year and month day.")

        if y < 1900 or y > 2100:
            raise FinError(
                "Date: year " + str(y) + " should be 1900 to 2100.")

        if d < 1:
            raise FinError("Date: Leap year. Day not valid.")

        leapYear = isLeapYear(y)

        if leapYear:
            if d > monthDaysLeapYear[m - 1]:
                raise FinError("Date: Leap year. Day not valid.")
        else:
            if d > monthDaysNotLeapYear[m - 1]:
                raise FinError("Date: Not Leap year. Day not valid.")

        self._y = y
        self._m = m
        self._d = d
        self._excelDate = 0

        # update the excel date used for doing lots of financial calculations
        self.refresh()

    ###########################################################################

    def fromDatetime(self, dt):
        ''' Construct a FinDate from a datetime as this is often needed if we
        receive inputs from other Python objects such as Pandas dataframes. '''

        finDate = FinDate(dt._y, dt._m, dt._d)
        return finDate

    ###########################################################################

    def refresh(self):
        ''' Update internal representation of date as number of days since the
        1st Jan 1900. This is same as Excel convention. '''

        dt = datetime.date(self._y, self._m, self._d)
        delta = dt - datetime.date(1900, 1, 1)
        self._excelDate = delta.days
        self._weekday = dt.weekday()

    ###########################################################################

    def __lt__(self, other):
        return self._excelDate < other._excelDate

    ###########################################################################

    def __gt__(self, other):
        return self._excelDate > other._excelDate

    ###########################################################################

    def __le__(self, other):
        return self._excelDate <= other._excelDate

    ###########################################################################

    def __ge__(self, other):
        return self._excelDate >= other._excelDate

    ###########################################################################

    def __sub__(self, other):
        return self._excelDate - other._excelDate

    ###########################################################################

    def __eq__(self, other):
        return self._excelDate == other._excelDate

    ###########################################################################

    def isWeekend(self):
        ''' returns True if the date falls on a weekend. '''

        if self._weekday == FinDate.SAT or self._weekday == FinDate.SUN:
            return True

        return False

    ###########################################################################

    def addDays(self, numDays):
        ''' Returns a new date that is numDays after the FinDate. '''

        dt = datetime.date(self._y, self._m, self._d)
        dt = dt + datetime.timedelta(days=numDays)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = FinDate(y, m, d)
        return newDt

    ###########################################################################

    def addWorkDays(self, numDays):
        ''' Returns a new date that is numDays working days after FinDate. '''

        if type(numDays) is not int:
            raise FinError("Num days must be an integer")

        dt = datetime.date(self._y, self._m, self._d)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = FinDate(y, m, d)

        while numDays > 0:
            dt = dt + datetime.timedelta(days=1)
            d = dt.day
            m = dt.month
            y = dt.year
            newDt = FinDate(y, m, d)

            if newDt.isWeekend() is False:
                numDays = numDays - 1

        return newDt

    ###########################################################################

    def addMonths(self, mm):
        ''' Returns a new date that is mm months after the FinDate. '''

        m = self._m + mm
        y = self._y
        d = self._d

        while m > 12:
            m = m - 12
            y += 1

        while m < 1:
            m = m + 12
            y -= 1

        newDt = FinDate(y, m, d)
        return newDt

    ##########################################################################

    def nextCDSDate(self, mm=0):
        ''' Returns a CDS date that is mm months after the FinDate. If no
        argument is supplied then the next CDS date after today is returned.'''

        nextDate = self.addMonths(mm)

        y = nextDate._y
        m = nextDate._m
        d = nextDate._d

        d_cds = 20
        y_cds = y
        m_cds = 999

        if m == 12 and d >= 20:
            m_cds = 3
            y_cds = y + 1
        elif m == 10 or m == 11 or m == 12:
            m_cds = 12
        elif m == 9 and d >= 20:
            m_cds = 12
        elif m == 7 or m == 8 or m == 9:
            m_cds = 9
        elif m == 6 and d >= 20:
            m_cds = 9
        elif m == 4 or m == 5 or m == 6:
            m_cds = 6
        elif m == 3 and d >= 20:
            m_cds = 6
        elif m == 1 or m == 2 or m == 3:
            m_cds = 3

        cdsDate = FinDate(y_cds, m_cds, d_cds)
        return cdsDate

    ##########################################################################

    def thirdWednesdayOfMonth(self, m, y):
        ''' For a specific month and year this returns the day number of the
            3rd Wednesday by scanning through dates in the third week. '''

        d_start = 14
        d_end = 21

        for d in range(d_start, d_end):
            immDate = FinDate(y, m, d)
            if immDate._weekday == self.WED:
                return d

        # Should never reach this line but just to be defensive
        raise FinError("Third Wednesday not found")

    ##########################################################################

    def nextIMMDate(self):
        ''' This function returns the next IMM date after the current date
            This is a 3rd Wednesday of Jun, March, Sep or December '''

        y = self._y
        m = self._m
        d = self._d

        y_imm = y

        if m == 12 and d >= self.thirdWednesdayOfMonth(m, y):
            m_imm = 3
            y_imm = y + 1
        elif m == 10 or m == 11 or m == 12:
            m_imm = 12
        elif m == 9 and d >= self.thirdWednesdayOfMonth(m, y):
            m_imm = 12
        elif m == 7 or m == 8 or m == 9:
            m_imm = 9
        elif m == 6 and d >= self.thirdWednesdayOfMonth(m, y):
            m_imm = 9
        elif m == 4 or m == 5 or m == 6:
            m_imm = 6
        elif m == 3 and d >= self.thirdWednesdayOfMonth(m, y):
            m_imm = 6
        elif m == 1 or m == 2 or m == 3:
            m_imm = 3

        d_imm = self.thirdWednesdayOfMonth(m_imm, y_imm)

        immDate = FinDate(y_imm, m_imm, d_imm)
        return immDate

    ###########################################################################

    def addTenor(self, tenor):
        ''' Return the date following the FinDate by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. '''

        if type(tenor) != str:
            raise ValueError("Tenor must be a string e.g. '5Y'")

        tenor = tenor.upper()
        DAYS = 1
        WEEKS = 2
        MONTHS = 3
        YEARS = 4

        periodType = 0
        numPeriods = 0

        if tenor[-1] == "D":
            periodType = DAYS
        elif tenor[-1] == "W":
            periodType = WEEKS
        elif tenor[-1] == "M":
            periodType = MONTHS
        elif tenor[-1] == "Y":
            periodType = YEARS
        else:
            raise ValueError("Unknown tenor type in " + tenor)

        numPeriods = int(tenor[0:-1])

        newDate = FinDate(self._y, self._m, self._d)

        if periodType == DAYS:
            for i in range(0, numPeriods):
                newDate = newDate.addDays(1)
        elif periodType == WEEKS:
            for i in range(0, numPeriods):
                newDate = newDate.addDays(7)
        elif periodType == MONTHS:
            for i in range(0, numPeriods):
                newDate = newDate.addMonths(1)
        elif periodType == YEARS:
            for i in range(0, numPeriods):
                newDate = newDate.addMonths(12)

        return newDate

    ###########################################################################

    def date(self):

        return datetime.date(self._y, self._m, self._d)

    ###########################################################################

    def __str__(self):
        dateStr = ""
        dateStr += shortDayNames[self._weekday]
        dateStr += " " + str(self._d) + " "
        dateStr += shortMonthNames[self._m - 1]
        dateStr += " " + str(self._y)
        return dateStr

    ###########################################################################

    def print(self):
        print(self.__str__())

    ###########################################################################
