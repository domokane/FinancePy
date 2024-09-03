##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from collections.abc import Iterable
from functools import partial
from enum import Enum

from typing import Union
from numba import njit

import numpy as np
import datetime

import math

# from financepy.utils.error import FinError
from .error import FinError
from .tenor import Tenor, TenorUnit


###############################################################################


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
g_date_type_format = DateFormatTypes.UK_LONG


def set_date_format(format_type):
    """Function that sets the global date format type."""
    global g_date_type_format
    g_date_type_format = format_type


###############################################################################


short_day_names = ["MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN"]
long_day_names = [
    "MONDAY",
    "TUESDAY",
    "WEDNESDAY",
    "THURSDAY",
    "FRIDAY",
    "SATURDAY",
    "SUNDAY",
]
short_month_names = [
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
]
longMonthNames = [
    "JANUARY",
    "FEBRUARY",
    "MARCH",
    "APRIL",
    "MAY",
    "JUNE",
    "JULY",
    "AUGUST",
    "SEPTEMBER",
    "OCTOBER",
    "NOVEMBER",
    "DECEMBER",
]

month_days_not_leap_year = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_days_leap_year = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

###############################################################################

# TODO: Fix this - it has stopped working
# @njit(boolean(int64), fastmath=True, cache=True)


def is_leap_year(y: int):
    """Test whether year y is a leap year - if so return True, else False"""
    leap_year = (y % 4 == 0) and (y % 100 != 0) or (y % 400 == 0)
    return leap_year


###############################################################################


def parse_dt(date_str, date_format):
    dt_obj = datetime.datetime.strptime(date_str, date_format)
    return dt_obj.day, dt_obj.month, dt_obj.year


###############################################################################
# CREATE DATE COUNTER
###############################################################################


g_dt_counter_list = None
g_start_year = 1900
g_end_year = 2100


def calculate_list():
    """Calculate list of dates so that we can do quick lookup to get the
    number of dates since 1 Jan 1900 (inclusive) BUT TAKING INTO ACCOUNT THE
    FACT THAT EXCEL MISTAKENLY CALLS 1900 A LEAP YEAR. For us, agreement with
    Excel is more important than this leap year error and in any case, we will
    not usually be calculating day differences with start dates before 28 Feb
    1900. Note that Excel inherited this "BUG" from LOTUS 1-2-3."""

    day_counter = 0
    max_days = 0

    global g_dt_counter_list
    global g_start_year
    global g_end_year

    g_dt_counter_list = []

    idx = -1  # the first element will be idx=0

    for yy in range(1900, g_end_year + 1):

        # DO NOT CHANGE THIS FOR AGREEMENT WITH EXCEL WHICH ASSUMES THAT 1900
        # WAS A LEAP YEAR AND THAT 29 FEB 1900 ACTUALLY HAPPENED. A LOTUS BUG.
        if yy == 1900:
            leap_year = True
        else:
            leap_year = is_leap_year(yy)

        for mm in range(1, 13):

            if leap_year is True:
                max_days = month_days_leap_year[mm - 1]
            else:
                max_days = month_days_not_leap_year[mm - 1]

            for _ in range(1, max_days + 1):
                idx += 1
                day_counter += 1
                if yy >= g_start_year:
                    g_dt_counter_list.append(day_counter)

            for _ in range(max_days, 31):
                idx += 1
                if yy >= g_start_year:
                    g_dt_counter_list.append(-999)


###############################################################################
# The index in these functions is not the Excel date index used as the
# internal representation of the date but the index of that date in the
# padded date object used to store the dates in a way that allows for a
# quick lookup. Do not confuse them as you will find they are out by months
###############################################################################


@njit(fastmath=True, cache=True)
def date_index(d, m, y):
    """Calculate the index of a date assuming 31 days in all months"""
    idx = (y - g_start_year) * 12 * 31 + (m - 1) * 31 + (d - 1)
    return idx


###############################################################################


@njit(fastmath=True, cache=True)
def date_from_index(idx):
    """Reverse mapping from index to date. Take care with numba as it can do
    weird rounding on the integer. Seems OK now."""
    y = int(g_start_year + idx / 12 / 31)
    m = 1 + int((idx - (y - g_start_year) * 12 * 31) / 31)
    d = 1 + idx - (y - g_start_year) * 12 * 31 - (m - 1) * 31
    return (d, m, y)


###############################################################################


@njit(fastmath=True, cache=True)
def weekday(day_count):
    """Converts day count to a weekday based on Excel date"""
    week_day = (day_count + 5) % 7
    return week_day


###############################################################################


def vectorisation_helper(func):
    def wrapper(self_, other):
        if isinstance(other, Iterable):
            # Store the type of other, then cast the output to be the same type
            output_type = type(other)
            f = partial(func, self_)
            return output_type(map(f, other))
        return func(self_, other)

    return wrapper


###############################################################################


class Date:
    """A date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance."""

    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6

    ###########################################################################

    def __init__(self, d, m, y, hh=0, mm=0, ss=0):
        """Create a date given a day of month, month and year. The arguments
        must be in the order of day (of month), month number and then the year.
        The year must be a 4-digit number greater than or equal to 1900. The
        user can also supply an hour, minute and second for intraday work.

        Example Input:
        start_dt = Date(1, 1, 2018)
        """

        global g_start_year
        global g_end_year

        # If the date has been entered as y, m, d we flip it to d, m, y
        # This message should be removed after a few releases
        if d >= g_start_year and d < g_end_year and y > 0 and y <= 31:
            raise FinError(
                "Date arguments must now be in the order Date(dd, mm, yyyy)"
            )

        if g_dt_counter_list is None:
            calculate_list()

        if y < 1900:
            raise FinError("Year cannot be before 1900")

        # Resize date list dynamically if required
        if y < g_start_year:
            g_start_year = y
            calculate_list()

        if y > g_end_year:
            g_end_year = y
            calculate_list()

        if y < g_start_year or y > g_end_year:
            raise FinError(
                "Date: year "
                + str(y)
                + " should be "
                + str(g_start_year)
                + " to "
                + str(g_end_year)
            )

        if d < 1:
            raise FinError("Date: Leap year. Day not valid.")

        leap_year = is_leap_year(y)

        if leap_year:
            if d > month_days_leap_year[m - 1]:
                print(d, m, y)
                raise FinError("Date: Leap year. Day not valid.")
        else:
            if d > month_days_not_leap_year[m - 1]:
                print(d, m, y)
                raise FinError("Date: Not Leap year. Day not valid.")

        if hh < 0 or hh > 23:
            raise FinError("Hours must be in range 0-23")

        if mm < 0 or mm > 59:
            raise FinError("Minutes must be in range 0-59")

        if ss < 0 or ss > 59:
            raise FinError("Seconds must be in range 0-59")

        self.y = y
        self.m = m
        self.d = d

        self.hh = hh
        self.mm = mm
        self.ss = ss

        self.excel_dt = 0  # This is a float as it includes intraday time

        # update the Excel date used for doing lots of financial calculations
        self._refresh()

        day_fraction = self.hh / 24.0
        day_fraction += self.mm / 24.0 / 60.0
        day_fraction += self.ss / 24.0 / 60.0 / 60.0

        self.excel_dt += day_fraction  # This is float - holds intraday time

    ###########################################################################

    @classmethod
    def from_string(cls, date_string, format_string):
        """Create a Date from a date and format string.
        Example Input:
        start_dt = Date('1-1-2018', '%d-%m-%Y')"""

        d, m, y = parse_dt(date_string, format_string)
        return cls(d, m, y)

    ###########################################################################

    @classmethod
    def from_date(cls, date: [datetime.date, np.datetime64]):
        """Create a Date from a python datetime.date object or from a
        Numpy datetime64 object.
        Example Input:
        start_dt = Date.from_dt(datetime.date(2022, 11, 8))"""

        if isinstance(date, datetime.date):
            d, m, y = date.day, date.month, date.year
            return cls(d, m, y)

        if isinstance(date, np.datetime64):
            time_stamp = (
                date - np.datetime64("1970-01-01T00:00:00")
            ) / np.timedelta64(1, "s")

            date = datetime.datetime.utcfromtime_stamp(time_stamp)
            d, m, y = date.day, date.month, date.year
            return cls(d, m, y)

    ###########################################################################

    def _refresh(self):
        """Update internal representation of date as number of days since the
        1st Jan 1900. This is same as Excel convention."""
        idx = date_index(self.d, self.m, self.y)
        days_since_first_jan_1900 = g_dt_counter_list[idx]
        wd = weekday(days_since_first_jan_1900)
        self.excel_dt = days_since_first_jan_1900
        self.weekday = wd

    ###########################################################################

    @vectorisation_helper
    def __gt__(self, other):
        return self.excel_dt > other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __lt__(self, other):
        return self.excel_dt < other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __ge__(self, other):
        return self.excel_dt >= other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __le__(self, other):
        return self.excel_dt <= other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __sub__(self, other):
        return self.excel_dt - other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __rsub__(self, other):
        return self.excel_dt - other.excel_dt

    ###########################################################################

    @vectorisation_helper
    def __eq__(self, other):
        return self.excel_dt == other.excel_dt

    ###########################################################################

    def __hash__(self):
        return hash(self.excel_dt)

    ###########################################################################

    def is_weekend(self):
        """returns True if the date falls on a weekend."""

        if self.weekday == Date.SAT or self.weekday == Date.SUN:
            return True

        return False

    ###########################################################################

    def is_eom(self):
        """returns True if this date falls on a month end."""

        y = self.y
        m = self.m
        d = self.d

        leap_year = is_leap_year(y)

        if leap_year:
            if d == month_days_leap_year[m - 1]:
                return True
        else:
            if d == month_days_not_leap_year[m - 1]:
                return True

        return False

    ###########################################################################

    def eom(self):
        """returns last date of month of this date."""

        y = self.y
        m = self.m

        leap_year = is_leap_year(y)

        if leap_year:
            last_day = month_days_leap_year[m - 1]
            return Date(last_day, m, y)
        else:
            last_day = month_days_not_leap_year[m - 1]
            return Date(last_day, m, y)

        return False

    ###########################################################################

    def add_hours(self, hours):
        """Returns a new date that is h hours after the Date."""

        if hours < 0:
            raise FinError("Number of hours must be positive")

        start_hour = self.hh
        final_hour = start_hour + hours
        days = int(final_hour / 24)
        hour = final_hour % 24

        # Move forward a specific number of days
        dt_1 = self.add_days(days)

        # On that date we then move to the correct hour
        dt_2 = Date(dt_1.d, dt_1.m, dt_1.y, hour, dt_1.mm, dt_1.ss)
        return dt_2

    ###########################################################################

    def add_days(self, num_days: int = 1):
        """Returns a new date that is num_days after the Date. I also make
        it possible to go backwards a number of days."""

        idx = date_index(self.d, self.m, self.y)

        step = +1
        if num_days < 0:
            step = -1

        while num_days != 0:
            idx += step
            if g_dt_counter_list[idx] > 0:
                num_days -= step

        (d, m, y) = date_from_index(idx)
        new_dt = Date(d, m, y)
        return new_dt

    ###########################################################################

    def add_weekdays(self, num_days: int):
        """Returns a new date that is num_days working days after Date. Note
        that only weekends are taken into account. Other Holidays are not. If
        you want to include regional holidays then use add_business_days from
        the FinCalendar class."""

        # TODO: REMOVE DATETIME DEPENDENCE HERE

        end_dt = self

        if isinstance(num_days, int) is False:
            raise FinError("Num days must be an integer")

        positive_num_days = num_days > 0
        num_days = abs(num_days)

        # 5 week days make up a week
        old_logic = False

        if old_logic is True:

            num_weeks = int(num_days / 5)
            remaining_days = num_days % 5

            if self.weekday == Date.SAT:
                weekend_adjust = 1
            elif self.weekday == Date.SUN:
                weekend_adjust = 0
            else:
                weekend_adjust = 2

            if positive_num_days is True:
                if self.weekday + remaining_days > self.FRI:
                    # add weekend
                    remaining_days += weekend_adjust

                return self.add_days(num_weeks * 7 + remaining_days)

            else:

                if self.weekday - remaining_days < self.MON:
                    # add weekend
                    remaining_days += weekend_adjust

            return self.add_days(-(num_weeks * 7 + remaining_days))

        else:  # new logic

            num_days_left = num_days
            end_dt = self

            while num_days_left > 0:

                if positive_num_days is True:
                    end_dt = end_dt.add_days(1)
                else:
                    end_dt = end_dt.add_days(-1)

                if end_dt.weekday == Date.SAT or end_dt.weekday == Date.SUN:
                    pass
                else:
                    num_days_left -= 1

            return end_dt

    ###########################################################################

    def add_months(self, mm: (list, int)):
        """Returns a new date that is mm months after the Date. If mm is an
        integer or float you get back a single date. If mm is a vector you get
        back a vector of dates."""

        num_months = 1
        scalar_flag = False

        if isinstance(mm, int) or isinstance(mm, float):
            mm_vector = [mm]
            scalar_flag = True
        else:
            mm_vector = mm

        num_months = len(mm_vector)

        date_list = []

        for i in range(0, num_months):

            mmi = mm_vector[i]

            # If I get a float I check it has no decimal places
            if int(mmi) != mmi:
                raise FinError("Must only pass integers or float integers.")

            mmi = int(mmi)

            d = self.d
            m = self.m + mmi
            y = self.y

            while m > 12:
                m = m - 12
                y += 1

            while m < 1:
                m = m + 12
                y -= 1

            leap_year = is_leap_year(y)

            if leap_year:
                if d > month_days_leap_year[m - 1]:
                    d = month_days_leap_year[m - 1]
            else:
                if d > month_days_not_leap_year[m - 1]:
                    d = month_days_not_leap_year[m - 1]

            new_dt = Date(d, m, y)
            date_list.append(new_dt)

        if scalar_flag is True:
            return date_list[0]
        else:
            return date_list

    ###########################################################################

    def add_years(self, yy: (np.ndarray, float)):
        """Returns a new date that is yy years after the Date. If yy is an
        integer or float you get back a single date. If yy is a list you get
        back a vector of dates."""

        num_years = 1
        scalar_flag = False

        if isinstance(yy, int) or isinstance(yy, float):
            yy_vector = [yy]
            scalar_flag = True
        else:
            yy_vector = yy

        num_years = len(yy_vector)

        date_list = []

        for i in range(0, num_years):

            yyi = yy_vector[i]

            # If yyi is not a whole month I adjust for days using average
            # number of days in a month which is 365.242/12
            daysInMonth = 365.242 / 12.0

            mmi = int(yyi * 12.0)
            ddi = int((yyi * 12.0 - mmi) * daysInMonth)
            new_dt = self.add_months(mmi)
            new_dt = new_dt.add_days(ddi)

            date_list.append(new_dt)

        if scalar_flag is True:
            return date_list[0]
        else:
            return date_list

    ##########################################################################

    def next_cds_date(self, mm: int = 0):
        """Returns a CDS date that is mm months after the Date. If no
        argument is supplied then the next CDS date after today is returned."""

        next_dt = self.add_months(mm)

        y = next_dt.y
        m = next_dt.m
        d = next_dt.d

        d_cds = 20
        y_cds = y
        m_cds = 999

        if m == 12 and d >= 20:
            m_cds = 3
            y_cds = y + 1
        elif m in (10, 11, 12):
            m_cds = 12
        elif m == 9 and d >= 20:
            m_cds = 12
        elif m in (7, 8, 9):
            m_cds = 9
        elif m == 6 and d >= 20:
            m_cds = 9
        elif m in (4, 5, 6):
            m_cds = 6
        elif m == 3 and d >= 20:
            m_cds = 6
        elif m in (1, 2, 3):
            m_cds = 3

        cds_dt = Date(d_cds, m_cds, y_cds)
        return cds_dt

    ##########################################################################

    def third_wednesday_of_month(self, m: int, y: int):
        """For a specific month and year this returns the day number of the
        3rd Wednesday by scanning through dates in the third week."""

        # Suppose 1st is Weds then 8th is Wed and 15th is 3rd Wednesday
        # Suppose 1st is Thur then 7th is Wed and 14th is 2nd Wednesday so 21
        # is 3rd so earliest and latest dates are 15th and 21st

        d_start = 15
        d_end = 21

        for d in range(d_start, d_end + 1):
            imm_dt = Date(d, m, y)
            if imm_dt.weekday == self.WED:
                return d

        # Should never reach this line but just to be defensive
        raise FinError("Third Wednesday not found")

    ##########################################################################

    def next_imm_date(self):
        """This function returns the next IMM date after the current date
        This is a 3rd Wednesday of Jun, March, Sep or December. For an
        IMM contract the IMM date is the First Delivery Date of the
        futures contract."""

        y = self.y
        m = self.m
        d = self.d

        y_imm = y

        if m == 12 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 3
            y_imm = y + 1
        elif m in (10, 11, 12):
            m_imm = 12
        elif m == 9 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 12
        elif m in (7, 8, 9):
            m_imm = 9
        elif m == 6 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 9
        elif m in (4, 5, 6):
            m_imm = 6
        elif m == 3 and d >= self.third_wednesday_of_month(m, y):
            m_imm = 6
        elif m in (1, 2, 3):
            m_imm = 3

        d_imm = self.third_wednesday_of_month(m_imm, y_imm)

        imm_dt = Date(d_imm, m_imm, y_imm)
        return imm_dt

    ###########################################################################

    def add_tenor(self, tenor: Union[list, str, Tenor]):
        """Return the date following the Date by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. The date is NOT weekend or holiday calendar adjusted. This must
        be done AFTERWARDS."""

        list_flag = False

        if isinstance(tenor, list) is True:
            list_flag = True
            for ten in tenor:
                if (
                    isinstance(ten, str) is False
                    and isinstance(ten, Tenor) is False
                ):
                    raise FinError(
                        "Tenor must be a string e.g. '5Y' or a Tenor object"
                    )
        else:
            if (
                isinstance(tenor, str) is True
                or isinstance(tenor, Tenor) is True
            ):
                tenor = [tenor]
            else:
                raise FinError(
                    "Tenor must be a string e.g. '5Y' or a Tenor object"
                )

        new_dts = []

        for tenor_string in tenor:

            tenor_obj = Tenor.as_tenor(str_or_tenor=tenor_string)

            new_dt = Date(self.d, self.m, self.y)

            if tenor_obj._units == TenorUnit.DAYS:
                for _ in range(0, abs(tenor_obj._num_periods)):
                    new_dt = new_dt.add_days(
                        math.copysign(1, tenor_obj._num_periods)
                    )
            elif tenor_obj._units == TenorUnit.WEEKS:
                for _ in range(0, abs(tenor_obj._num_periods)):
                    new_dt = new_dt.add_days(
                        math.copysign(7, tenor_obj._num_periods)
                    )
            elif tenor_obj._units == TenorUnit.MONTHS:
                for _ in range(0, abs(tenor_obj._num_periods)):
                    new_dt = new_dt.add_months(
                        math.copysign(1, tenor_obj._num_periods)
                    )

                # in case we landed on a 28th Feb and lost the month day
                # we add this logic
                y = new_dt.y
                m = new_dt.m
                d = min(self.d, new_dt.eom().d)
                new_dt = Date(d, m, y)

            elif tenor_obj._units == TenorUnit.YEARS:
                for _ in range(0, abs(tenor_obj._num_periods)):
                    new_dt = new_dt.add_months(
                        math.copysign(12, tenor_obj._num_periods)
                    )

            new_dts.append(new_dt)

        if list_flag is True:
            return new_dts
        else:
            return new_dts[0]

    ###########################################################################

    def datetime(self):
        """Returns a datetime of the date"""

        # Remember that datetime likes inputs in opposite order
        return datetime.date(self.y, self.m, self.d)

    ###########################################################################
    # TODO: Find elegant way to return long and short strings
    ###########################################################################

    def str(self):
        """returns a formatted string of the date"""
        date_str = ""

        if self.d < 10:
            date_str += "0" + str(self.d) + ""
        else:
            date_str += "" + str(self.d) + ""

        date_str += short_month_names[self.m - 1]
        date_str += "" + str(self.y)
        return date_str

    ###########################################################################

    def __repr__(self):
        """returns a formatted string of the date"""

        global g_date_type_format

        day_name_str = short_day_names[self.weekday]

        if self.d < 10:
            day_str = "0" + str(self.d)
        else:
            day_str = "" + str(self.d)

        if self.m < 10:
            short_month_str = "0" + str(self.m)
        else:
            short_month_str = str(self.m)

        long_month_str = short_month_names[self.m - 1]

        short_year_str = str(self.y)[2:]
        long_year_str = str(self.y)

        if g_date_type_format == DateFormatTypes.UK_LONGEST:

            sep = " "
            date_str = (
                day_name_str
                + " "
                + day_str
                + sep
                + long_month_str
                + sep
                + long_year_str
            )
            return date_str

        elif g_date_type_format == DateFormatTypes.UK_LONG:

            sep = "-"
            date_str = day_str + sep + long_month_str + sep + long_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.UK_MEDIUM:

            sep = "/"
            date_str = day_str + sep + short_month_str + sep + long_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.UK_SHORT:

            sep = "/"
            date_str = day_str + sep + short_month_str + sep + short_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.US_LONGEST:

            sep = " "
            date_str = (
                day_name_str
                + " "
                + long_month_str
                + sep
                + day_str
                + sep
                + long_year_str
            )
            return date_str

        elif g_date_type_format == DateFormatTypes.US_LONG:

            sep = "-"
            date_str = long_month_str + sep + day_str + sep + long_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.US_MEDIUM:

            sep = "-"
            date_str = short_month_str + sep + day_str + sep + long_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.US_SHORT:

            sep = "-"
            date_str = short_month_str + sep + day_str + sep + short_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.BLOOMBERG:

            sep = "/"
            date_str = short_month_str + sep + day_str + sep + short_year_str
            return date_str

        elif g_date_type_format == DateFormatTypes.DATETIME:

            sep = "/"

            if self.hh < 10:
                hour_str = "0" + str(self.hh)
            else:
                hour_str = str(self.hh)

            if self.mm < 10:
                minute_str = "0" + str(self.mm)
            else:
                minute_str = str(self.mm)

            if self.ss < 10:
                second_str = "0" + str(self.ss)
            else:
                second_str = str(self.ss)

            time_str = hour_str + ":" + minute_str + ":" + second_str
            date_str = day_str + sep + short_month_str + sep + long_year_str
            date_str = date_str + " " + time_str
            return date_str

        else:

            raise FinError("Unknown date format")

    ###########################################################################
    # REMOVE THIS

    def _print(self):
        """prints formatted string of the date."""
        print(self)


###############################################################################
# Date functions that are not class members but are useful
###############################################################################


def daily_working_day_schedule(start_dt: Date, end_dt: Date):
    """Returns a list of working dates between start_dt and end_dt.
    This function should be replaced by dateRange once add_tenor allows
    for working days."""
    date_list = []

    dt = start_dt
    date_list.append(dt)
    while dt < end_dt:
        dt = dt.add_weekdays(1)
        date_list.append(dt)

    return date_list


###############################################################################


def datediff(d1: Date, d2: Date):
    """Calculate the number of days between two Findates."""
    dd = d2.excel_dt - d1.excel_dt
    return int(dd)


###############################################################################


def from_datetime(dt: Date):
    """Construct a Date from a datetime as this is often needed if we
    receive inputs from other Python objects such as Pandas dataframes."""

    fin_dt = Date(dt.day, dt.month, dt.year)
    return fin_dt


###############################################################################


def days_in_month(m, y):
    """Get the number of days in the month (1-12) of a given year y."""

    if m < 1 or m > 12:
        raise FinError("Month must be 1-12")

    if is_leap_year(y) is False:
        return month_days_not_leap_year[m - 1]
    else:
        return month_days_leap_year[m - 1]


###############################################################################


def date_range(start_dt: Date, end_dt: Date, tenor: str = "1D"):
    """Returns a list of dates between start_dt (inclusive)
    and end_dt (inclusive). The tenor represents the distance between two
    consecutive dates and is set to daily by default."""

    if start_dt > end_dt:
        return []

    date_list = []

    dt = start_dt
    while dt < end_dt:
        date_list.append(dt)
        dt = dt.add_tenor(tenor)
    date_list.append(end_dt)

    return date_list


###############################################################################


def test_type():
    global g_date_type_format
    print("TEST TYPE", g_date_type_format)


###############################################################################
