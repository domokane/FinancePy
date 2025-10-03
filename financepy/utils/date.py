##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from collections.abc import Iterable
from functools import partial
import datetime
import math


from typing import Union
from numba import njit

import numpy as np

from .error import FinError
from .tenor import Tenor, TenorUnit
from .date_format import DateFormatTypes, get_date_format
from .date_arrays import short_day_names, short_month_names
from .date_arrays import month_days_leap_year, month_days_not_leap_year

########################################################################################
# --- Precomputed day offsets ---
START_YEAR = 1900
END_YEAR = 2200

# Precompute year offsets
year_offsets_list = []
days = 0
for y in range(START_YEAR, END_YEAR + 1):
    year_offsets_list.append(days)
    if y == 1900:
        days += 366
    else:
        days += 366 if (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0) else 365

YEAR_OFFSETS = np.array(year_offsets_list, dtype=np.int32)

# Cumulative month days
NONLEAP_CUM = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334], dtype=np.int32)
LEAP_CUM    = np.array([0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335], dtype=np.int32)

########################################################################################

@njit(cache=True, fastmath=True)
def is_leap(y: int) -> bool:
    return (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)

########################################################################################

@njit(cache=True, fastmath=True)
def excel_from_ymd(d: int, m: int, y: int) -> int:
    y_off = YEAR_OFFSETS[y - START_YEAR]
    cum = LEAP_CUM if (y == 1900 or is_leap(y)) else NONLEAP_CUM
    return y_off + cum[m - 1] + d

########################################################################################

@njit(cache=True)
def ymd_from_excel(excel_dt: int):
    """Excel day number (1 = 1900-01-01) -> (d, m, y)."""
    days = excel_dt

    # Excel bug: treat 1900 as a leap year
    if days > 59:  # i.e. after 28 Feb 1900
        days -= 1  # skip the fake Feb 29

    # Now proceed with normal proleptic Gregorian calculation
    y = 1900
    while True:
        leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
        days_in_year = 366 if leap else 365
        if days > days_in_year:
            days -= days_in_year
            y += 1
        else:
            break

    # find month/day
    leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
    mdays = month_days_leap_year if leap else month_days_not_leap_year
    m = 1
    while days > mdays[m-1]:
        days -= mdays[m-1]
        m += 1
    d = days
    return d, m, y

########################################################################################

def is_leap_year(y: int):
    """Test whether year y is a leap year - if so return True, else False"""
    leap_year = (y % 4 == 0) and (y % 100 != 0) or (y % 400 == 0)
    return leap_year


########################################################################################


def parse_dt(date_str, date_format):
    """
    Parse a date string into day, month, and year components.
    """
    dt_obj = datetime.datetime.strptime(date_str, date_format)
    return dt_obj.day, dt_obj.month, dt_obj.year

########################################################################################


@njit(fastmath=True, cache=True)
def weekday(day_count) -> int:
    """Converts day count to a weekday based on Excel date"""
    week_day = (day_count + 5) % 7
    return week_day


########################################################################################


def vectorisation_helper(func):
    """Decorator to allow a function to be vectorised over an iterable"""

    def wrapper(self_, other):
        if isinstance(other, Iterable):
            # Store the type of other, then cast the output to be the same type
            output_type = type(other)
            f = partial(func, self_)
            return output_type(map(f, other))
        return func(self_, other)

    return wrapper


########################################################################################


class Date:
    """A date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance."""

    __slots__ = ("d", "m", "y", "hh", "mm", "ss", "excel_dt", "weekday")

    MON, TUE, WED, THU, FRI, SAT, SUN = range(7)

    ####################################################################################

class Date:
    """A date class to manage dates that is simple to use and includes a
    number of useful date functions used frequently in Finance."""

    __slots__ = ("d", "m", "y", "hh", "mm", "ss", "excel_dt", "weekday")

    MON, TUE, WED, THU, FRI, SAT, SUN = range(7)

    ####################################################################################

    def __init__(self, d, m, y, hh=0, mm=0, ss=0):
        # validation checks (as you already had)
        if 1900 <= d <= 2100 and 1 <= y <= 31:
            raise FinError("Date arguments must be in order Date(day, month, year)")
        if y < 1900:
            raise FinError("Year cannot be before 1900")
        ...
        # set fields
        self.y, self.m, self.d = y, m, d
        self.hh, self.mm, self.ss = hh, mm, ss
        self._refresh()
        day_fraction = (hh / 24.0) + (mm / (24.0 * 60.0)) + (ss / (24.0 * 3600.0))
        self.excel_dt += day_fraction

    ####################################################################################

    @staticmethod
    def _make(d, m, y, hh=0, mm=0, ss=0):
        """Normal safe constructor."""
        return Date(d, m, y, hh, mm, ss)

    ####################################################################################

    @classmethod
    def _make_fast(cls, d, m, y, excel_dt):
        """Fast constructor when excel_dt is already known."""
        obj = cls.__new__(cls)
        obj.d, obj.m, obj.y = d, m, y
        obj.hh = obj.mm = obj.ss = 0
        obj.excel_dt = excel_dt
        obj.weekday = (int(excel_dt) + 5) % 7
        return obj

    ####################################################################################

    @classmethod
    def from_ymd_excel(cls, d, m, y, excel_dt):
        obj = cls.__new__(cls)  # allocate without __init__
        obj.d, obj.m, obj.y = d, m, y
        obj.hh = obj.mm = obj.ss = 0
        obj.excel_dt = excel_dt
        obj.weekday = (int(excel_dt) + 5) % 7
        return obj

    ####################################################################################

    @staticmethod
    def from_excel(excel_dt: int) -> "Date":
        d, m, y = ymd_from_excel(int(excel_dt))
        return Date(d, m, y)

    ####################################################################################

    @classmethod
    def from_string(cls, date_string, format_string):
        """Create a Date from a date and format string.
        Example Input:
        start_dt = Date('1-1-2018', '%d-%m-%Y')"""

        d, m, y = parse_dt(date_string, format_string)
        return cls(d, m, y)

    ####################################################################################

    @classmethod
    def from_date(cls, date: Union[datetime.date, np.datetime64]):
        """Create a Date from a python datetime.date object or from a
        Numpy datetime64 object.
        Example Input:
        start_dt = Date.from_dt(datetime.date(2022, 11, 8))"""

        if isinstance(date, datetime.date):
            d, m, y = date.day, date.month, date.year
            return cls(d, m, y)

        if isinstance(date, np.datetime64):
            time_stamp = (date - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(
                1, "s"
            )

            date = datetime.datetime.utcfromtimestamp(time_stamp)
            d, m, y = date.day, date.month, date.year
            return cls(d, m, y)

    ####################################################################################

    def _refresh(self):
        """Update Excel day count and weekday."""
        self.excel_dt = excel_from_ymd(self.d, self.m, self.y)
        self.weekday = (int(self.excel_dt) + 5) % 7

    ####################################################################################

    @vectorisation_helper
    def __gt__(self, other):
        return self.excel_dt > other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __lt__(self, other):
        return self.excel_dt < other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __ge__(self, other):
        return self.excel_dt >= other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __le__(self, other):
        return self.excel_dt <= other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __sub__(self, other):
        return self.excel_dt - other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __rsub__(self, other):
        return self.excel_dt - other.excel_dt

    ####################################################################################

    @vectorisation_helper
    def __eq__(self, other):
        return self.excel_dt == other.excel_dt

    ####################################################################################

    def __hash__(self):
        return hash(self.excel_dt)

    ####################################################################################

    def is_weekend(self):
        """returns True if the date falls on a weekend."""

        if self.weekday == Date.SAT or self.weekday == Date.SUN:
            return True

        return False

    ####################################################################################

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

    ####################################################################################

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

    ####################################################################################

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

    ####################################################################################

    def add_days(self, num_days: int = 1):
        new_excel_dt = int(self.excel_dt) + int(num_days)
        d, m, y = ymd_from_excel(new_excel_dt)
        new_dt = Date._make_fast(d, m, y, new_excel_dt)
        return new_dt

    ####################################################################################

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

    ####################################################################################

    def add_months(self, mm: Union[list, int]) -> Union["Date", list]:
        """Returns a new date that is mm months after the Date. If mm is an
        integer or float you get back a single date. If mm is a vector you get
        back a vector of dates."""

        scalar_flag = False

        if isinstance(mm, (int, float)):
            mm_vector = [mm]
            scalar_flag = True
        else:
            mm_vector = mm

        date_list = []

        for mmi in mm_vector:

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

            excel_dt = excel_from_ymd(d, m, y)
            new_dt = Date._make_fast(d, m, y, excel_dt)
            date_list.append(new_dt)

        if scalar_flag is True:
            return date_list[0]

        return date_list

    ####################################################################################

    def add_years(self, yy: Union[np.ndarray, float]):
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
            days_in_month_float = 365.242 / 12.0

            mmi = int(yyi * 12.0)
            ddi = int((yyi * 12.0 - mmi) * days_in_month_float)

            new_dt = self.add_months(mmi)
            new_dt = new_dt.add_days(ddi)
            date_list.append(new_dt)

        if scalar_flag is True:
            return date_list[0]
        else:
            return date_list

    ####################################################################################

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

        excel_dt = excel_from_ymd(d_cds, m_cds, y_cds)
        cds_dt = Date._make_fast(d_cds, m_cds, y_cds, excel_dt)
        return cds_dt

    ####################################################################################

    def third_wednesday_of_month(self, m: int, y: int):
        """For a specific month and year this returns the day number of the
        3rd Wednesday by scanning through dates in the third week."""

        # Suppose 1st is Weds then 8th is Wed and 15th is 3rd Wednesday
        # Suppose 1st is Thur then 7th is Wed and 14th is 2nd Wednesday so 21
        # is 3rd so earliest and latest dates are 15th and 21st

        d_start = 15
        d_end = 21

        for d in range(d_start, d_end + 1):
            excel_dt = excel_from_ymd(d, m, y)
            if (excel_dt + 5) % 7 == self.WED:
                return d

        # Should never reach this line but just to be defensive
        raise FinError("Third Wednesday not found")

    ####################################################################################

    def next_imm_date(self):
        """This function returns the next IMM date after the current date
        This is a 3rd Wednesday of Jun, March, Sep or December. For an
        IMM contract the IMM date is the First Delivery Date of the
        futures contract."""

        y = self.y
        m = self.m
        d = self.d

        m_imm = None
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
        excel_dt = excel_from_ymd(d_imm, m_imm, y_imm)
        imm_dt = Date._make_fast(d_imm, m_imm, y_imm, excel_dt)
        return imm_dt

    ####################################################################################

    def add_tenor(self, tenor: Union[list, str, Tenor]):
        """Return the date following the Date by a period given by the
        tenor which is a string consisting of a number and a letter, the
        letter being d, w, m , y for day, week, month or year. This is case
        independent. For example 10Y means 10 years while 120m also means 10
        years. The date is NOT weekend or holiday calendar adjusted. This must
        be done AFTERWARDS."""

        list_flag = isinstance(tenor, list)
        if not list_flag:
            tenor = [tenor]

        new_dts = []

        for tenor_item in tenor:
            tenor_obj = Tenor.as_tenor(str_or_tenor=tenor_item)

            # start with current date
            d, m, y = self.d, self.m, self.y
            excel_dt = int(self.excel_dt)

            if tenor_obj.units == TenorUnit.DAYS:
                new_excel_dt = excel_dt + int(tenor_obj.num_periods)
                d, m, y = ymd_from_excel(new_excel_dt)
                new_dt = Date._make_fast(d, m, y, new_excel_dt)

            elif tenor_obj.units == TenorUnit.WEEKS:
                new_excel_dt = excel_dt + int(tenor_obj.num_periods) * 7
                d, m, y = ymd_from_excel(new_excel_dt)
                new_dt = Date._make_fast(d, m, y, new_excel_dt)

            elif tenor_obj.units == TenorUnit.MONTHS:
                # reuse add_months (already _make_fast optimized)
                new_dt = self.add_months(int(tenor_obj.num_periods))

                # handle case when Feb truncates to 28/29 (preserve original d if possible)
                y, m = new_dt.y, new_dt.m
                d = min(self.d, new_dt.eom().d)
                excel_dt = excel_from_ymd(d, m, y)
                new_dt = Date._make_fast(d, m, y, excel_dt)

            elif tenor_obj.units == TenorUnit.YEARS:
                # just months × 12
                new_dt = self.add_months(int(tenor_obj.num_periods) * 12)

            else:
                raise FinError("Unsupported tenor unit: " + str(tenor_obj.units))

            new_dts.append(new_dt)

        return new_dts if list_flag else new_dts[0]

    ####################################################################################

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

    ####################################################################################

    def __repr__(self):
        """returns a formatted string of the date"""

        global G_DATE_TYPE_FORMAT

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

        if get_date_format() == DateFormatTypes.UK_LONGEST:

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

        elif get_date_format() == DateFormatTypes.UK_LONG:

            sep = "-"
            date_str = day_str + sep + long_month_str + sep + long_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.UK_MEDIUM:

            sep = "/"
            date_str = day_str + sep + short_month_str + sep + long_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.UK_SHORT:

            sep = "/"
            date_str = day_str + sep + short_month_str + sep + short_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.US_LONGEST:

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

        elif get_date_format() == DateFormatTypes.US_LONG:

            sep = "-"
            date_str = long_month_str + sep + day_str + sep + long_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.US_MEDIUM:

            sep = "-"
            date_str = short_month_str + sep + day_str + sep + long_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.US_SHORT:

            sep = "-"
            date_str = short_month_str + sep + day_str + sep + short_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.BLOOMBERG:

            sep = "/"
            date_str = short_month_str + sep + day_str + sep + short_year_str
            return date_str

        elif get_date_format() == DateFormatTypes.DATETIME:

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

        raise FinError("Unknown date format")

    ###########################################################################
    # REMOVE THIS

    def _print(self):
        """prints formatted string of the date."""
        print(self)


########################################################################################
# Date functions that are not class members but are useful
########################################################################################


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


########################################################################################


def datediff(d1: Date, d2: Date):
    """Calculate the number of days between two Findates."""
    dd = d2.excel_dt - d1.excel_dt
    return int(dd)


########################################################################################


def from_datetime(dt: Date):
    """Construct a Date from a datetime as this is often needed if we
    receive inputs from other Python objects such as Pandas dataframes."""

    fin_dt = Date(dt.day, dt.month, dt.year)
    return fin_dt


########################################################################################


def days_in_month(m, y):
    """Get the number of days in the month (1-12) of a given year y."""

    if m < 1 or m > 12:
        raise FinError("Month must be 1-12")

    if is_leap_year(y) is False:
        return month_days_not_leap_year[m - 1]
    else:
        return month_days_leap_year[m - 1]


########################################################################################


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


########################################################################################


