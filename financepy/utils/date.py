##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit
import numpy as np
import datetime

from financepy.utils.error import FinError

###############################################################################

from enum import Enum


class DateFormatTypes(Enum):
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
gDateFormatType = DateFormatTypes.UK_LONG


def set_date_format(format_type):
    """ Function that sets the global date format type. """
    global gDateFormatType
    gDateFormatType = format_type

###############################################################################


short_day_names = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
long_day_names = [
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


def is_leap_year(y: int):
    """ Test whether year y is a leap year - if so return True, else False """
    leap_year = ((y % 4 == 0) and (y % 100 != 0) or (y % 400 == 0))
    return leap_year

###############################################################################


def parse_date(date_str, date_format):
    dt_obj = datetime.datetime.strptime(date_str, date_format)
    return dt_obj.day, dt_obj.month, dt_obj.year

###############################################################################
# CREATE DATE COUNTER
###############################################################################


gDateCounterList = None
gStartYear = 1900
gEndYear = 2100


def calculate_list():
    """ Calculate list of dates so that we can do quick lookup to get the
    number of dates since 1 Jan 1900 (inclusive) BUT TAKING INTO ACCOUNT THE
    FACT THAT EXCEL MISTAKENLY CALLS 1900 A LEAP YEAR. For us, agreement with
    Excel is more important than this leap year error and in any case, we will
    not usually be calculating day differences with start dates before 28 Feb
    1900. Note that Excel inherited this "BUG" from LOTUS 1-2-3. """

    day_counter = 0
    max_days = 0

    global gDateCounterList
    global gStartYear
    global gEndYear

    gDateCounterList = []

    idx = -1  # the first element will be idx=0

    for yy in range(1900, gEndYear+1):

        # DO NOT CHANGE THIS FOR AGREEMENT WITH EXCEL WHICH ASSUMES THAT 1900
        # WAS A LEAP YEAR AND THAT 29 FEB 1900 ACTUALLY HAPPENED. A LOTUS BUG.
        if yy == 1900:
            leap_year = True
        else:
            leap_year = is_leap_year(yy)

        for mm in range(1, 13):

            if leap_year is True:
                max_days = monthDaysLeapYear[mm-1]
            else:
                max_days = monthDaysNotLeapYear[mm-1]

            for _ in range(1, max_days+1):
                idx += 1
                day_counter += 1
                if yy >= gStartYear:
                    gDateCounterList.append(day_counter)

            for _ in range(max_days, 31):
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
def date_index(d, m, y):
    idx = (y-gStartYear) * 12 * 31 + (m-1) * 31 + (d-1)
    return idx

###############################################################################


@njit(fastmath=True, cache=True)
def date_from_index(idx):
    """ Reverse mapping from index to date. Take care with numba as it can do
    weird rounding on the integer. Seems OK now. """
    y = int(gStartYear + idx/12/31)
    m = 1 + int((idx - (y-gStartYear) * 12 * 31) / 31)
    d = 1 + idx - (y-gStartYear) * 12 * 31 - (m-1) * 31
    return (d, m, y)

###############################################################################


@njit(fastmath=True, cache=True)
def weekday(day_count):
    weekday = (day_count+5) % 7
    return weekday

###############################################################################


class Date():
    """ A date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance. """

    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6

    ###########################################################################

    def __init__(self, d, m, y, hh=0, mm=0, ss=0):
        """ Create a date given a day of month, month and year. The arguments
        must be in the order of day (of month), month number and then the year.
        The year must be a 4-digit number greater than or equal to 1900. The
        user can also supply an hour, minute and second for intraday work.

        Example Input:
        start_date = Date(1, 1, 2018)
        """

        global gStartYear
        global gEndYear

        # If the date has been entered as y, m, d we flip it to d, m, y
        # This message should be removed after a few releases
        if d >= gStartYear and d < gEndYear and y > 0 and y <= 31:
            raise FinError(
                "Date arguments must now be in the order Date(dd, mm, yyyy)")

        if gDateCounterList is None:
            calculate_list()

        if y < 1900:
            raise FinError("Year cannot be before 1900")

        # Resize date list dynamically if required
        if y < gStartYear:
            gStartYear = y
            calculate_list()

        if y > gEndYear:
            gEndYear = y
            calculate_list()

        if y < gStartYear or y > gEndYear:
            raise FinError(
                "Date: year " + str(y) + " should be " + str(gStartYear) +
                " to " + str(gEndYear))

        if d < 1:
            raise FinError("Date: Leap year. Day not valid.")

        leap_year = is_leap_year(y)

        if leap_year:
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

        self._excel_date = 0  # This is a float as it includes intraday time

        # update the excel date used for doing lots of financial calculations
        self._refresh()

        dayFraction = self._hh/24.0
        dayFraction += self._mm/24.0/60.0
        dayFraction += self._ss/24.0/60.0/60.0

        self._excel_date += dayFraction  # This is a float as it includes intraday time

    ###########################################################################

    @classmethod
    def from_string(cls, date_string, formatString):
        """  Create a Date from a date and format string.
        Example Input:
        start_date = Date('1-1-2018', '%d-%m-%Y') """

        d, m, y = parse_date(date_string, formatString)
        return cls(d, m, y)

    ###########################################################################

    def _refresh(self):
        """ Update internal representation of date as number of days since the
        1st Jan 1900. This is same as Excel convention. """

        idx = date_index(self._d, self._m, self._y)
        daysSinceFirstJan1900 = gDateCounterList[idx]
        wd = weekday(daysSinceFirstJan1900)
        self._excel_date = daysSinceFirstJan1900
        self._weekday = wd

    ###########################################################################

    def __lt__(self, other):
        return self._excel_date < other._excel_date

    ###########################################################################

    def __gt__(self, other):
        return self._excel_date > other._excel_date

    ###########################################################################

    def __le__(self, other):
        return self._excel_date <= other._excel_date

    ###########################################################################

    def __ge__(self, other):
        return self._excel_date >= other._excel_date

    ###########################################################################

    def __sub__(self, other):
        return self._excel_date - other._excel_date

    ###########################################################################

    def __eq__(self, other):
        return self._excel_date == other._excel_date

    ###########################################################################

    def is_weekend(self):
        """ returns True if the date falls on a weekend. """

        if self._weekday == Date.SAT or self._weekday == Date.SUN:
            return True

        return False

    ###########################################################################

    def is_eom(self):
        """ returns True if this date falls on a month end. """

        y = self._y
        m = self._m
        d = self._d

        leap_year = is_leap_year(y)

        if leap_year:
            if d == monthDaysLeapYear[m - 1]:
                return True
        else:
            if d == monthDaysNotLeapYear[m - 1]:
                return True

        return False

    ###########################################################################

    def eom(self):
        """ returns last date of month of this date. """

        y = self._y
        m = self._m

        leap_year = is_leap_year(y)

        if leap_year:
            lastDay = monthDaysLeapYear[m - 1]
            return Date(lastDay, m, y)
        else:
            lastDay = monthDaysNotLeapYear[m - 1]
            return Date(lastDay, m, y)

        return False

    ###########################################################################

    def add_hours(self, hours):
        """ Returns a new date that is h hours after the Date. """

        if hours < 0:
            raise FinError("Number of hours must be positive")

        startHour = self._hh
        finalHour = startHour + hours
        days = int(finalHour/24)
        hour = finalHour % 24

        # Move forward a specific number of days
        dt1 = self.add_days(days)

        # On that date we then move to the correct hour
        dt2 = Date(dt1._d, dt1._m, dt1._y, hour, dt1._mm, dt1._ss)
        return dt2

    ###########################################################################

    def add_days(self,
                 numDays: int = 1):
        """ Returns a new date that is numDays after the Date. I also make
        it possible to go backwards a number of days. """

        idx = date_index(self._d, self._m, self._y)

        step = +1
        if numDays < 0:
            step = -1

        while numDays != 0:
            idx += step
            if gDateCounterList[idx] > 0:
                numDays -= step

        (d, m, y) = date_from_index(idx)
        newDt = Date(d, m, y)
        return newDt

    ###########################################################################

    def add_weekdays(self,
                     numDays: int):
        """ Returns a new date that is numDays working days after Date. Note
        that only weekends are taken into account. Other Holidays are not. If
        you want to include regional holidays then use add_business_days from
        the FinCalendar class. """

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

            return self.add_days(numWeeks * 7 + remainingDays)
        else:
            if (self._weekday - remainingDays < self.MON):
                # add weekend
                remainingDays += 2

            return self.add_days(-(numWeeks * 7 + remainingDays))

    ###########################################################################

    def add_months(self,
                   mm: (list, int)):
        """ Returns a new date that is mm months after the Date. If mm is an
        integer or float you get back a single date. If mm is a vector you get
        back a vector of dates."""

        num_months = 1
        scalarFlag = False

        if isinstance(mm, int) or isinstance(mm, float):
            mmVector = [mm]
            scalarFlag = True
        else:
            mmVector = mm

        num_months = len(mmVector)

        dateList = []

        for i in range(0, num_months):

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

            leap_year = is_leap_year(y)

            if leap_year:
                if d > monthDaysLeapYear[m - 1]:
                    d = monthDaysLeapYear[m-1]
            else:
                if d > monthDaysNotLeapYear[m - 1]:
                    d = monthDaysNotLeapYear[m-1]

            newDt = Date(d, m, y)
            dateList.append(newDt)

        if scalarFlag is True:
            return dateList[0]
        else:
            return dateList

    ###########################################################################

    def add_years(self,
                  yy: (np.ndarray, float)):
        """ Returns a new date that is yy years after the Date. If yy is an
        integer or float you get back a single date. If yy is a list you get
        back a vector of dates."""

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
            newDt = self.add_months(mmi)
            newDt = newDt.add_days(ddi)

            dateList.append(newDt)

        if scalarFlag is True:
            return dateList[0]
        else:
            return dateList

    ##########################################################################

    def next_cds_date(self,
                      mm: int = 0):
        """ Returns a CDS date that is mm months after the Date. If no
        argument is supplied then the next CDS date after today is returned."""

        next_date = self.add_months(mm)

        y = next_date._y
        m = next_date._m
        d = next_date._d

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

        cdsDate = Date(d_cds, m_cds, y_cds)
        return cdsDate

    ##########################################################################

    def third_wednesday_of_month(self,
                                 m: int,  # Month number
                                 y: int):  # Year number
        """ For a specific month and year this returns the day number of the
            3rd Wednesday by scanning through dates in the third week. """

        # Suppose 1st is Weds then 8th is Wed and 15th is 3rd Wed
        # Suppose 1st is Thur then 7th is Wed and 14th is 2nd Wed so 21 is 3rd
        # so earliest and latest dates are 15th and 21st

        d_start = 15
        d_end = 21

        for d in range(d_start, d_end+1):
            immDate = Date(d, m, y)
            if immDate._weekday == self.WED:
                return d

        # Should never reach this line but just to be defensive
        raise FinError("Third Wednesday not found")

    ##########################################################################

    def next_imm_date(self):
        """ This function returns the next IMM date after the current date
            This is a 3rd Wednesday of Jun, March, Sep or December. For an 
            IMM contract the IMM date is the First Delivery Date of the
            futures contract. """

        y = self._y
        m = self._m
        d = self._d

        y_imm = y

        if m == 12 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 3
            y_imm = y + 1
        elif m == 10 or m == 11 or m == 12:
            m_imm = 12
        elif m == 9 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 12
        elif m == 7 or m == 8 or m == 9:
            m_imm = 9
        elif m == 6 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 9
        elif m == 4 or m == 5 or m == 6:
            m_imm = 6
        elif m == 3 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 6
        elif m == 1 or m == 2 or m == 3:
            m_imm = 3

        d_imm = self.third_wednesday_of_month(m_imm, y_imm)

        immDate = Date(d_imm, m_imm, y_imm)
        return immDate

    ###########################################################################

    def add_tenor(self,
                  tenor: (list, str)):
        """ Return the date following the Date by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. The date is NOT weekend or holiday calendar adjusted. This must
        be done AFTERWARDS. """

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
            num_periods = 0

            if tenStr == "ON":   # overnight - should be used only if spot days = 0
                periodType = DAYS
                num_periods = 1
            elif tenStr == "TN":  # overnight - should be used when spot days > 0
                periodType = DAYS
                num_periods = 1
            elif tenStr[-1] == "D":
                periodType = DAYS
                num_periods = int(tenStr[0:-1])
            elif tenStr[-1] == "W":
                periodType = WEEKS
                num_periods = int(tenStr[0:-1])
            elif tenStr[-1] == "M":
                periodType = MONTHS
                num_periods = int(tenStr[0:-1])
            elif tenStr[-1] == "Y":
                periodType = YEARS
                num_periods = int(tenStr[0:-1])
            else:
                raise FinError("Unknown tenor type in " + tenor)

            newDate = Date(self._d, self._m, self._y)

            if periodType == DAYS:
                for _ in range(0, num_periods):
                    newDate = newDate.add_days(1)
            elif periodType == WEEKS:
                for _ in range(0, num_periods):
                    newDate = newDate.add_days(7)
            elif periodType == MONTHS:
                for _ in range(0, num_periods):
                    newDate = newDate.add_months(1)

                # in case we landed on a 28th Feb and lost the month day we add this logic
                y = newDate._y
                m = newDate._m
                d = min(self._d, newDate.eom()._d)
                newDate = Date(d, m, y)

            elif periodType == YEARS:
                for _ in range(0, num_periods):
                    newDate = newDate.add_months(12)

            newDates.append(newDate)

        if listFlag is True:
            return newDates
        else:
            return newDates[0]

    ###########################################################################

    def datetime(self):
        """ Returns a datetime of the date """

        # Remember that datetime likes inputs in opposite order
        return datetime.date(self._y, self._m, self._d)

    ###########################################################################
    # TODO: Find elegant way to return long and short strings
    ###########################################################################

    def str(self, format):
        """ returns a formatted string of the date """
        date_str = ""

        if self._d < 10:
            date_str += "0" + str(self._d) + ""
        else:
            date_str += "" + str(self._d) + ""

        date_str += shortMonthNames[self._m - 1]
        date_str += "" + str(self._y)
        return date_str

    ###########################################################################

    def __repr__(self):
        """ returns a formatted string of the date """

        global gDateFormatType

        dayNameStr = short_day_names[self._weekday]

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

        if gDateFormatType == DateFormatTypes.UK_LONGEST:

            sep = " "
            date_str = dayNameStr + " " + dayStr + sep + longMonthStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.UK_LONG:

            sep = "-"
            date_str = dayStr + sep + longMonthStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.UK_MEDIUM:

            sep = "/"
            date_str = dayStr + sep + shortMonthStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.UK_SHORT:

            sep = "/"
            date_str = dayStr + sep + shortMonthStr + sep + shortYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.US_LONGEST:

            sep = " "
            date_str = dayNameStr + " " + longMonthStr + sep + dayStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.US_LONG:

            sep = "-"
            date_str = longMonthStr + sep + dayStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.US_MEDIUM:

            sep = "-"
            date_str = shortMonthStr + sep + dayStr + sep + longYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.US_SHORT:

            sep = "-"
            date_str = shortMonthStr + sep + dayStr + sep + shortYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.BLOOMBERG:

            sep = "/"
            date_str = shortMonthStr + sep + dayStr + sep + shortYearStr
            return date_str

        elif gDateFormatType == DateFormatTypes.DATETIME:

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
            date_str = dayStr + sep + shortMonthStr + sep + longYearStr
            date_str = date_str + " " + timeStr
            return date_str

        else:

            raise FinError("Unknown date format")

    ###########################################################################
    # REMOVE THIS

    def _print(self):
        """ prints formatted string of the date. """
        print(self)


###############################################################################
# Date functions that are not class members but are useful
###############################################################################


def daily_working_day_schedule(self,
                               start_date: Date,
                               end_date: Date):
    """ Returns a list of working dates between start_date and end_date.
    This function should be replaced by dateRange once add_tenor allows
    for working days. """
    dateList = []

    dt = start_date
    dateList.append(dt)
    while dt < end_date:
        dt = dt.add_weekdays(1)
        dateList.append(dt)

    return dateList

###############################################################################


def datediff(d1: Date,
             d2: Date):
    """ Calculate the number of days between two Findates. """
    dd = (d2._excel_date - d1._excel_date)
    return int(dd)

###############################################################################


def from_datetime(dt: Date):
    """ Construct a Date from a datetime as this is often needed if we
    receive inputs from other Python objects such as Pandas dataframes. """

    finDate = Date(dt.day, dt.month, dt.year)
    return finDate

###############################################################################


def days_in_month(m, y):
    """ Get the number of days in the month (1-12) of a given year y. """

    if m < 1 or m > 12:
        raise FinError("Month must be 1-12")

    if is_leap_year(y) is False:
        return monthDaysNotLeapYear[m-1]
    else:
        return monthDaysLeapYear[m-1]

###############################################################################


def date_range(start_date: Date,
               end_date: Date,
               tenor: str = "1D"):
    """ Returns a list of dates between start_date (inclusive)
    and end_date (inclusive). The tenor represents the distance between two
    consecutive dates and is set to daily by default. """

    if start_date > end_date:
        return []

    dateList = []

    dt = start_date
    while dt < end_date:
        dateList.append(dt)
        dt = dt.add_tenor(tenor)
    dateList.append(end_date)

    return dateList

###############################################################################


def test_type():
    global gDateFormatType
    print("TEST TYPE", gDateFormatType)

###############################################################################
