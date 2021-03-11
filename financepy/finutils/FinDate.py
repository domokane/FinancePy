##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import datetime
from .FinError import FinError
from numba import njit, boolean, int64
import numpy as np

###############################################################################    

from enum import Enum

class FinDateFormatTypes(Enum):
    BLOOMBERG = 1
    US_SHORT = 2
    US_MEDIUM = 3
    US_LONG = 4
    US_LONGEST = 5
    UK_SHORT = 6
    UK_MEDIUM = 7
    UK_LONG = 8
    UK_LONGEST = 9
    DATETIME = 10
    

# Set the default
gDateFormatType = FinDateFormatTypes.UK_LONG

def setDateFormatType(formatType):
    ''' Function that sets the global date format type. '''
    global gDateFormatType
    gDateFormatType = formatType

###############################################################################    

ENFORCE_DAY_FIRST = True

###############################################################################

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

# TODO: Fix this - it has stopped working
# @njit(boolean(int64), fastmath=True, cache=True)
def isLeapYear(y: int):
    ''' Test whether year y is a leap year - if so return True, else False '''
    leapYear = ((y % 4 == 0) and (y % 100 != 0) or (y % 400 == 0))
    return leapYear

###############################################################################


def parse_date(dateStr, dateFormat):
    dt_obj = datetime.datetime.strptime(dateStr, dateFormat)
    return dt_obj.day, dt_obj.month, dt_obj.year

###############################################################################
# CREATE DATE COUNTER
###############################################################################


gDateCounterList = None
gStartYear = 1900
gEndYear = 2100

def calculateList():
    ''' Calculate list of dates so that we can do quick lookup to get the
    number of dates since 1 Jan 1900 (inclusive) BUT TAKING INTO ACCOUNT THE
    FACT THAT EXCEL MISTAKENLY CALLS 1900 A LEAP YEAR. For us, agreement with
    Excel is more important than this leap year error and in any case, we will
    not usually be calculating day differences with start dates before 28 Feb
    1900. Note that Excel inherited this "BUG" from LOTUS 1-2-3. '''

    dayCounter = 0
    maxDays = 0
    
    global gDateCounterList
    global gStartYear
    global gEndYear

#    print("Calculating list between", gStartYear, "and", gEndYear)

    gDateCounterList = []

    idx = -1  # the first element will be idx=0

    for yy in range(1900, gEndYear+1):

        # DO NOT CHANGE THIS FOR AGREEMENT WITH EXCEL WHICH ASSUMES THAT 1900
        # WAS A LEAP YEAR AND THAT 29 FEB 1900 ACTUALLY HAPPENED. A LOTUS BUG.
        if yy == 1900:
            leapYear = True
        else:
            leapYear = isLeapYear(yy)

        for mm in range(1, 13):

            if leapYear is True:
                maxDays = monthDaysLeapYear[mm-1]
            else:
                maxDays = monthDaysNotLeapYear[mm-1]

            for _ in range(1, maxDays+1):
                idx += 1
                dayCounter += 1
                if yy >= gStartYear:
                    gDateCounterList.append(dayCounter)

            for _ in range(maxDays, 31):
                idx += 1
                if yy >= gStartYear:
                    gDateCounterList.append(-999)

###############################################################################
# The index in these functions is not the excel date index used as the 
# internal representation of the date but the index of that date in the
# padded date object used to store the dates in a way that allows for a
# quick lookup. Do not confuse them as you will find they are out by months
###############################################################################

@njit(fastmath=True, cache=True)
def dateIndex(d, m, y):
    idx = (y-gStartYear) * 12 * 31 + (m-1) * 31 + (d-1)
    return idx

###############################################################################

@njit(fastmath=True, cache=True)
def dateFromIndex(idx):
    ''' Reverse mapping from index to date. Take care with numba as it can do
    weird rounding on the integer. Seems OK now. '''
    y = int(gStartYear + idx/12/31)
    m = 1 + int((idx - (y-gStartYear) * 12 * 31) / 31)
    d = 1 + idx - (y-gStartYear) * 12 * 31 - (m-1) * 31
    return (d, m, y)

###############################################################################

@njit(fastmath=True, cache=True)
def weekDay(dayCount):
    weekday = (dayCount+5) % 7
    return weekday

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

    def __init__(self, d, m, y, hh=0, mm=0, ss=0):  
        ''' Create a date given a day of month, month and year. The arguments
        must be in the order of day (of month), month number and then the year.
        The year must be a 4-digit number greater than or equal to 1900. The
        user can also supply an hour, minute and second for intraday work.

        Example Input:
        start_date = FinDate(1, 1, 2018)
        '''

        global gStartYear
        global gEndYear

        # If the date has been entered as y, m, d we flip it to d, m, y
        if d >= gStartYear and d < gEndYear and y > 0 and y <= 31:
            raise FinError("Dates must be in the format FinDate(dd, mm, yyyy")

        if gDateCounterList is None:
            calculateList()

        if y < 1900:
            raise FinError("Year cannot be before 1900")

        # Resize date list dynamically if required
        if y < gStartYear:
            gStartYear = y
            calculateList()

        if y > gEndYear:
            gEndYear = y
            calculateList()

        if y < gStartYear or y > gEndYear:
            raise FinError(
                "Date: year " + str(y) + " should be " + str(gStartYear) +
                " to " + str(gEndYear))

        if d < 1:
            raise FinError("Date: Leap year. Day not valid.")

        leapYear = isLeapYear(y)

        if leapYear:
            if d > monthDaysLeapYear[m - 1]:
                print(d, m, y)
                raise FinError("Date: Leap year. Day not valid.")
        else:
            if d > monthDaysNotLeapYear[m - 1]:
                print(d, m, y)
                raise FinError("Date: Not Leap year. Day not valid.")

        if hh < 0 or hh > 23:
            raise FinError("Hours must be in range 0-23")

        if mm < 0 or mm > 59:
            raise FinError("Minutes must be in range 0-59")

        if ss < 0 or ss > 59:
            raise FinError("Seconds must be in range 0-59")

        self._y = y
        self._m = m
        self._d = d
        
        self._hh = hh
        self._mm = mm
        self._ss = ss

        self._excelDate = 0 # This is a float as it includes intraday time

        # update the excel date used for doing lots of financial calculations
        self._refresh()

        dayFraction = self._hh/24.0
        dayFraction += self._mm/24.0/60.0
        dayFraction += self._ss/24.0/60.0/60.0

        self._excelDate += dayFraction # This is a float as it includes intraday time

    ###########################################################################

    @classmethod
    def fromString(cls, dateString, formatString):
        '''  Create a FinDate from a date and format string.
        Example Input:
        start_date = FinDate('1-1-2018', '%d-%m-%Y') '''
 
        d, m, y = parse_date(dateString, formatString)
        return cls(d, m, y)

    ###########################################################################

    def _refresh(self):
        ''' Update internal representation of date as number of days since the
        1st Jan 1900. This is same as Excel convention. '''

        idx = dateIndex(self._d, self._m, self._y)
        daysSinceFirstJan1900 = gDateCounterList[idx]
        wd = weekDay(daysSinceFirstJan1900)
        self._excelDate = daysSinceFirstJan1900
        self._weekday = wd

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

    def isEOM(self):
        ''' returns True if this date falls on a month end. '''
        
        y = self._y
        m = self._m
        d = self._d

        leapYear = isLeapYear(y)

        if leapYear:
            if d == monthDaysLeapYear[m - 1]:
                return True
        else:
            if d == monthDaysNotLeapYear[m - 1]:
                return True
                    
        return False

    ###########################################################################

    def EOM(self):
        ''' returns last date of month of this date. '''

        y = self._y
        m = self._m

        leapYear = isLeapYear(y)

        if leapYear:
            lastDay = monthDaysLeapYear[m - 1]
            return FinDate(lastDay, m, y)
        else:
            lastDay = monthDaysNotLeapYear[m - 1]
            return FinDate(lastDay, m, y)

        return False

    ###########################################################################

    def addHours(self, hours):
        ''' Returns a new date that is h hours after the FinDate. '''

        if hours < 0:
            raise FinError("Number of hours must be positive")
        
        startHour = self._hh        
        finalHour = startHour + hours
        days = int(finalHour/24)
        hour = finalHour % 24

        # Move forward a specific number of days
        dt1 = self.addDays(days)

        # On that date we then move to the correct hour
        dt2 = FinDate(dt1._d, dt1._m, dt1._y, hour, dt1._mm, dt1._ss)                
        return dt2

    ###########################################################################

    def addDays(self,
                numDays: int = 1):
        ''' Returns a new date that is numDays after the FinDate. I also make
        it possible to go backwards a number of days. '''

        idx = dateIndex(self._d, self._m, self._y)

        step = +1
        if numDays < 0:
            step = -1

        while numDays != 0:
            idx += step
            if gDateCounterList[idx] > 0:
                numDays -= step

        (d, m, y) = dateFromIndex(idx)
        newDt = FinDate(d, m, y)
        return newDt

    ###########################################################################

    def addWeekDays(self,
                    numDays: int):
        ''' Returns a new date that is numDays working days after FinDate. Note
        that only weekends are taken into account. Other Holidays are not. If
        you want to include regional holidays then use addBusinessDays from
        the FinCalendar class. '''

        # TODO: REMOVE DATETIME DEPENDENCE HERE

        if isinstance(numDays, int) is False:
            raise FinError("Num days must be an integer")

        positiveNumDays = (numDays > 0)
        numDays = abs(numDays)
        
        # 5 week days make up a week
        numWeeks = int(numDays / 5)
        remainingDays = numDays % 5

        if (positiveNumDays):
            if (self._weekday + remainingDays > self.FRI):
                # add weekend
                remainingDays += 2

            return self.addDays(numWeeks * 7 + remainingDays)
        else:
            if (self._weekday - remainingDays < self.MON):
                # add weekend
                remainingDays += 2

            return self.addDays(-(numWeeks * 7 + remainingDays))

    ###########################################################################

    def addMonths(self,
                  mm: (list, int)):
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

            d = self._d
            m = self._m + mmi
            y = self._y

            while m > 12:
                m = m - 12
                y += 1

            while m < 1:
                m = m + 12
                y -= 1

            leapYear = isLeapYear(y)

            if leapYear:
                if d > monthDaysLeapYear[m - 1]:
                    d = monthDaysLeapYear[m-1]
            else:
                if d > monthDaysNotLeapYear[m - 1]:
                    d = monthDaysNotLeapYear[m-1]

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

        # Suppose 1st is Weds then 8th is Wed and 15th is 3rd Wed
        # Suppose 1st is Thur then 7th is Wed and 14th is 2nd Wed so 21 is 3rd
        # so earliest and latest dates are 15th and 21st

        d_start = 15
        d_end = 21

        for d in range(d_start, d_end+1):
            immDate = FinDate(d, m, y)
            if immDate._weekday == self.WED:
                return d

        # Should never reach this line but just to be defensive
        raise FinError("Third Wednesday not found")

    ##########################################################################

    def nextIMMDate(self):
        ''' This function returns the next IMM date after the current date
            This is a 3rd Wednesday of Jun, March, Sep or December. For an 
            IMM contract the IMM date is the First Delivery Date of the
            futures contract. '''

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
                 tenor: (list,str)):
        ''' Return the date following the FinDate by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. The date is NOT weekend or holiday calendar adjusted. This must
        be done AFTERWARDS. '''

        listFlag = False

        if isinstance(tenor, list) is True:
            listFlag = True
            for ten in tenor:
                if isinstance(ten, str) is False:
                    raise FinError("Tenor must be a string e.g. '5Y'")
        else:
            if isinstance(tenor, str) is True:
                tenor = [tenor]
            else:
                raise FinError("Tenor must be a string e.g. '5Y'")

        newDates = []

        for tenStr in tenor:            
            tenStr = tenStr.upper()
            DAYS = 1
            WEEKS = 2
            MONTHS = 3
            YEARS = 4
    
            periodType = 0
            numPeriods = 0
    
            if tenStr == "ON":   # overnight - should be used only if spot days = 0
                periodType = DAYS
                numPeriods = 1
            elif tenStr == "TN":  # overnight - should be used when spot days > 0
                periodType = DAYS
                numPeriods = 1
            elif tenStr[-1] == "D":
                periodType = DAYS
                numPeriods = int(tenStr[0:-1])
            elif tenStr[-1] == "W":
                periodType = WEEKS
                numPeriods = int(tenStr[0:-1])
            elif tenStr[-1] == "M":
                periodType = MONTHS
                numPeriods = int(tenStr[0:-1])
            elif tenStr[-1] == "Y":
                periodType = YEARS
                numPeriods = int(tenStr[0:-1])
            else:
                raise FinError("Unknown tenor type in " + tenor)
    
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

            newDates.append(newDate)

        if listFlag is True:
            return newDates
        else:
            return newDates[0]

    ###########################################################################

    def datetime(self):
        ''' Returns a datetime of the date '''
        return datetime.date(self._y, self._m, self._d)

    ###########################################################################
    # TODO: Find elegant way to return long and short strings
    ###########################################################################

    def str(self, format):
        ''' returns a formatted string of the date '''
        dateStr = ""

        if self._d < 10:
            dateStr += "0" + str(self._d) + ""
        else:
            dateStr += "" + str(self._d) + ""

        dateStr += shortMonthNames[self._m - 1]
        dateStr += "" + str(self._y)
        return dateStr

    ###########################################################################

    def __repr__(self):
        ''' returns a formatted string of the date '''

        global gDateFormatType
        
        dayNameStr = shortDayNames[self._weekday]

        if self._d < 10:
            dayStr = "0" + str(self._d)
        else:
            dayStr = "" + str(self._d)

        if self._m < 10:
            shortMonthStr = "0" + str(self._m)
        else:
            shortMonthStr = str(self._m)

        longMonthStr = shortMonthNames[self._m - 1]

        shortYearStr = str(self._y)[2:]
        longYearStr = str(self._y)

        if gDateFormatType == FinDateFormatTypes.UK_LONGEST:

            sep = " "
            dateStr = dayNameStr + " " + dayStr + sep + longMonthStr + sep + longYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.UK_LONG:

            sep = "-"
            dateStr = dayStr + sep + longMonthStr + sep + longYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.UK_MEDIUM:

            sep = "/"
            dateStr = dayStr + sep + shortMonthStr + sep + longYearStr
            return dateStr
    
        elif gDateFormatType == FinDateFormatTypes.UK_SHORT:

            sep = "/"
            dateStr = dayStr + sep + shortMonthStr + sep + shortYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.US_LONGEST:

            sep = " "
            dateStr = dayNameStr + " " + longMonthStr + sep + dayStr + sep + longYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.US_LONG:

            sep = "-"
            dateStr = longMonthStr + sep + dayStr + sep + longYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.US_MEDIUM:
            
            sep = "-"
            dateStr = shortMonthStr + sep + dayStr + sep + longYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.US_SHORT:

            sep = "-"
            dateStr = shortMonthStr + sep + dayStr + sep + shortYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.BLOOMBERG:

            sep = "/"
            dateStr = shortMonthStr + sep + dayStr + sep + shortYearStr
            return dateStr

        elif gDateFormatType == FinDateFormatTypes.DATETIME:

            sep = "/"

            if self._hh < 10:
                hourStr = "0" + str(self._hh)
            else:
                hourStr = str(self._hh)

            if self._mm < 10:
                minuteStr = "0" + str(self._mm)
            else:
                minuteStr = str(self._mm)

            if self._ss < 10:
                secondStr = "0" + str(self._ss)
            else:
                secondStr = str(self._ss)

            timeStr = hourStr + ":" + minuteStr + ":" + secondStr
            dateStr = dayStr + sep + shortMonthStr + sep + longYearStr
            dateStr = dateStr + " " + timeStr
            return dateStr

        else:
            
            raise FinError("Unknown date format")

    ###########################################################################
    # REMOVE THIS
    
    def _print(self):
        ''' prints formatted string of the date. '''
        print(self)


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
        dt = dt.addWeekDays(1)
        dateList.append(dt)

    return dateList

###############################################################################


def datediff(d1: FinDate,
             d2: FinDate):
    ''' Calculate the number of days between two Findates. '''
    dd = (d2._excelDate - d1._excelDate)
    return int(dd)

###############################################################################


def fromDatetime(dt: FinDate):
    ''' Construct a FinDate from a datetime as this is often needed if we
    receive inputs from other Python objects such as Pandas dataframes. '''

    finDate = FinDate(dt.day, dt.month, dt.year)
    return finDate

###############################################################################

def daysInMonth(m, y):
    ''' Get the number of days in the month (1-12) of a given year y. '''

    if m < 1 or m > 12:
        raise FinError("Month must be 1-12")
    
    if isLeapYear(y) is False:
        return monthDaysNotLeapYear[m-1]
    else:
        return monthDaysLeapYear[m-1]

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

def testType():
    global gDateFormatType
    print("TEST TYPE", gDateFormatType)
    
###############################################################################
