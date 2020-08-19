##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import datetime
from .FinError import FinError
from .FinMath import isLeapYear
import numpy as np

# from numba import njit, float64, int32

ENFORCE_DAY_FIRST = True

##########################################################################

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


class FinDate():
    ''' A date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance. '''

    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6

    ###########################################################################

    def __init__(self,
                 d: int,  # Day number in month with values from 1 to 31
                 m: int,  # Month number where January = 1, ..., December = 12
                 y: int):  # Year number which must be between 1900 and 2100
        ''' Create a date given a day of month, month and year. The arguments
        must be in the order of day (of month), month number and then the year.
        The year must be a 4-digit number greater than or equal to 1900. '''

        if d >= 1900 and d < 2100 and y > 0 and y <= 31:
            tmp = y
            y = d
            d = tmp

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
        self._refresh()

    ###########################################################################

    def _refresh(self):
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

    def addDays(self,
                numDays: int):
        ''' Returns a new date that is numDays after the FinDate. '''

        dt = datetime.date(self._y, self._m, self._d)
        dt = dt + datetime.timedelta(days=numDays)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = FinDate(d, m, y)
        return newDt

    ###########################################################################

    def addWorkDays(self,
                    numDays: int):
        ''' Returns a new date that is numDays working days after FinDate. '''

        if type(numDays) is not int:
            raise FinError("Num days must be an integer")

        if numDays < 0:
            raise FinError("Num days must be positive.")

        dt = datetime.date(self._y, self._m, self._d)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = FinDate(d, m, y)

        while numDays > 0:
            dt = dt + datetime.timedelta(days=1)
            d = dt.day
            m = dt.month
            y = dt.year
            newDt = FinDate(d, m, y)

            if newDt.isWeekend() is False:
                numDays = numDays - 1

        return newDt

    ###########################################################################

    def addMonths(self,
                  mm: int):
        ''' Returns a new date that is mm months after the FinDate. If mm is an
        integer or float you get back a single date. If mm is a vector you get
        back a vector of dates.'''

        numMonths = 1
        scalarFlag = False
        if isinstance(mm, int) or isinstance(mm, float):
            mmVector = [mm]
            scalarFlag = True
        else:
            mmVector = mm

        numMonths = len(mmVector)

        dateList = []

        for i in range(0, numMonths):

            mmi = mmVector[i]

            # If I get a float I check it has no decimal places
            if int(mmi) != mmi:
                raise FinError("Must only pass integers or float integers.")

            mmi = int(mmi)

            m = self._m + mmi
            y = self._y
            d = self._d

            while m > 12:
                m = m - 12
                y += 1

            while m < 1:
                m = m + 12
                y -= 1

            newDt = FinDate(d, m, y)
            dateList.append(newDt)

        if scalarFlag is True:
            return dateList[0]
        else:
            return dateList

    ###########################################################################

    def addYears(self,
                 yy: (np.ndarray, float)):
        ''' Returns a new date that is yy years after the FinDate. If yy is an
        integer or float you get back a single date. If yy is a list you get
        back a vector of dates.'''

        numYears = 1
        scalarFlag = False

        if isinstance(yy, int) or isinstance(yy, float):
            yyVector = [yy]
            scalarFlag = True
        else:
            yyVector = yy

        numYears = len(yyVector)

        dateList = []

        for i in range(0, numYears):

            yyi = yyVector[i]

            # If yyi is not a whole month I adjust for days using average
            # number of days in a month which is 365.242/12
            daysInMonth = 365.242/12.0

            mmi = int(yyi * 12.0)
            ddi = int((yyi * 12.0 - mmi) * daysInMonth)
            newDt = self.addMonths(mmi)
            newDt = newDt.addDays(ddi)

            dateList.append(newDt)

        if scalarFlag is True:
            return dateList[0]
        else:
            return dateList

    ##########################################################################

    def nextCDSDate(self,
                    mm: int = 0):
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

        cdsDate = FinDate(d_cds, m_cds, y_cds)
        return cdsDate

    ##########################################################################

    def thirdWednesdayOfMonth(self,
                              m: int,  # Month number
                              y: int):  # Year number
        ''' For a specific month and year this returns the day number of the
            3rd Wednesday by scanning through dates in the third week. '''

        d_start = 14
        d_end = 21

        for d in range(d_start, d_end):
            immDate = FinDate(d, m, y)
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

        immDate = FinDate(d_imm, m_imm, y_imm)
        return immDate

    ###########################################################################

    def addTenor(self,
                 tenor: str):
        ''' Return the date following the FinDate by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. '''

        if isinstance(tenor, str) is False:
            raise FinError("Tenor must be a string e.g. '5Y'")

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
            raise FinError("Unknown tenor type in " + tenor)

        numPeriods = int(tenor[0:-1])

        newDate = FinDate(self._d, self._m, self._y)

        if periodType == DAYS:
            for _ in range(0, numPeriods):
                newDate = newDate.addDays(1)
        elif periodType == WEEKS:
            for _ in range(0, numPeriods):
                newDate = newDate.addDays(7)
        elif periodType == MONTHS:
            for _ in range(0, numPeriods):
                newDate = newDate.addMonths(1)
        elif periodType == YEARS:
            for _ in range(0, numPeriods):
                newDate = newDate.addMonths(12)

        return newDate

    ###########################################################################

    def datetime(self):
        ''' Returns a datetime of the date '''
        return datetime.date(self._d, self._m, self._y)

    ###########################################################################

    def __repr__(self):
        ''' returns a formatted string of the date '''
        dateStr = ""
        dateStr += shortDayNames[self._weekday]
        dateStr += " " + str(self._d) + " "
        dateStr += shortMonthNames[self._m - 1]
        dateStr += " " + str(self._y)
        return dateStr

    ###########################################################################

    def print(self):
        ''' prints formatted string of the date. '''
        print(self)

    ###########################################################################


###############################################################################
# Date functions that are not class members but are useful
###############################################################################


def dailyWorkingDaySchedule(self,
                            startDate: FinDate,
                            endDate: FinDate):
    ''' Returns a list of working dates between startDate and endDate.
    This function should be replaced by dateRange once addTenor allows
    for working days. '''
    dateList = []

    dt = startDate
    dateList.append(dt)
    while dt < endDate:
        dt = dt.addWorkDays(1)
        dateList.append(dt)

    return dateList

###############################################################################


def datediff(d1: FinDate,
             d2: FinDate):
    ''' Calculate the number of days between two Findates. '''
    dd = (d2._excelDate - d1._excelDate)
    return dd

###############################################################################


def fromDatetime(dt: FinDate):
    ''' Construct a FinDate from a datetime as this is often needed if we
    receive inputs from other Python objects such as Pandas dataframes. '''

#    finDate = FinDate(dt.year, dt.month, dt.day)
    finDate = FinDate(dt.day, dt.month, dt.year)
    return finDate

###############################################################################


def dateRange(startDate: FinDate,
              endDate: FinDate,
              tenor: str = "1D"):
    ''' Returns a list of dates between startDate (inclusive)
    and endDate (inclusive). The tenor represents the distance between two
    consecutive dates and is set to daily by default. '''

    if startDate > endDate:
        return []

    dateList = []

    dt = startDate
    while dt < endDate:
        dateList.append(dt)
        dt = dt.addTenor(tenor)
    dateList.append(endDate)

    return dateList

###############################################################################
