########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from functools import lru_cache
import warnings

########################################################################################
# TODO: Do some timings and tidy up logic in adjustment function
########################################################################################

import datetime
from enum import Enum
from .date import Date
from .error import FinError

MIN_CALENDAR_YEAR = 1901
MAX_CALENDAR_YEAR = 2100

easter_monday_day = [
    98,
    90,
    103,
    95,
    114,
    106,
    91,
    111,
    102,
    87,
    107,
    99,
    83,
    103,
    95,
    115,
    99,
    91,
    111,
    96,
    87,
    107,
    92,
    112,
    103,
    95,
    108,
    100,
    91,
    111,
    96,
    88,
    107,
    92,
    112,
    104,
    88,
    108,
    100,
    85,
    104,
    96,
    116,
    101,
    92,
    112,
    97,
    89,
    108,
    100,
    85,
    105,
    96,
    109,
    101,
    93,
    112,
    97,
    89,
    109,
    93,
    113,
    105,
    90,
    109,
    101,
    86,
    106,
    97,
    89,
    102,
    94,
    113,
    105,
    90,
    110,
    101,
    86,
    106,
    98,
    110,
    102,
    94,
    114,
    98,
    90,
    110,
    95,
    86,
    106,
    91,
    111,
    102,
    94,
    107,
    99,
    90,
    103,
    95,
    115,
    106,
    91,
    111,
    103,
    87,
    107,
    99,
    84,
    103,
    95,
    115,
    100,
    91,
    111,
    96,
    88,
    107,
    92,
    112,
    104,
    95,
    108,
    100,
    92,
    111,
    96,
    88,
    108,
    92,
    112,
    104,
    89,
    108,
    100,
    85,
    105,
    96,
    116,
    101,
    93,
    112,
    97,
    89,
    109,
    100,
    85,
    105,
    97,
    109,
    101,
    93,
    113,
    97,
    89,
    109,
    94,
    113,
    105,
    90,
    110,
    101,
    86,
    106,
    98,
    89,
    102,
    94,
    114,
    105,
    90,
    110,
    102,
    86,
    106,
    98,
    111,
    102,
    94,
    114,
    99,
    90,
    110,
    95,
    87,
    106,
    91,
    111,
    103,
    94,
    107,
    99,
    91,
    103,
    95,
    115,
    107,
    91,
    111,
    103,
    88,
    108,
    100,
    85,
    105,
    96,
    109,
    101,
    93,
    112,
    97,
    89,
    109,
    93,
    113,
    105,
    90,
    109,
    101,
    86,
    106,
    97,
    89,
    102,
    94,
    113,
    105,
    90,
    110,
    101,
    86,
    106,
    98,
    110,
    102,
    94,
    114,
    98,
    90,
    110,
    95,
    86,
    106,
    91,
    111,
    102,
    94,
    107,
    99,
    90,
    103,
    95,
    115,
    106,
    91,
    111,
    103,
    87,
    107,
    99,
    84,
    103,
    95,
    115,
    100,
    91,
    111,
    96,
    88,
    107,
    92,
    112,
    104,
    95,
    108,
    100,
    92,
    111,
    96,
    88,
    108,
    92,
    112,
    104,
    89,
    108,
    100,
    85,
    105,
    96,
    116,
    101,
    93,
    112,
    97,
    89,
    109,
    100,
    85,
    105,
]

@lru_cache(maxsize=None)
def _warn_incomplete_calendar(calendar_name: str, year: int):
    warnings.warn(
        f"{calendar_name} calendar has incomplete holiday data "
        f"for year {year}.",
        RuntimeWarning,
        stacklevel=3,
    )

class BusDayAdjustTypes(Enum):
    """Enum for business day adjustment types."""

    NONE = 1
    FOLLOWING = 2
    MODIFIED_FOLLOWING = 3
    PRECEDING = 4
    MODIFIED_PRECEDING = 5


class CalendarTypes(Enum):
    """Enum for calendar types."""

    NONE = 1
    WEEKEND = 2

    # Legacy / country calendars
    AUSTRALIA = 3
    CANADA = 4
    FRANCE = 5
    GERMANY = 6
    ITALY = 7
    JAPAN = 8
    NEW_ZEALAND = 9
    NORWAY = 10
    SWEDEN = 11
    SWITZERLAND = 12

    # Core financial-centre calendars
    TARGET = 13
    UNITED_STATES = 14
    NEW_YORK = 15
    LONDON = 16
    UNITED_KINGDOM = 17

    # Derivatives / benchmark calendars
    US_GOVERNMENT_SECURITIES = 18
    US_FEDERAL_RESERVE = 19

    SYDNEY = 20
    AUSTRALIA_RITS = 21
    AUD_FX_SETTLEMENT = 22

    TORONTO = 23
    TOKYO = 24
    ZURICH = 25
    HONG_KONG = 26
    SINGAPORE = 27


class DateGenRuleTypes(Enum):
    """Enum for date generation rule types."""

    FORWARD = 1
    BACKWARD = 2


########################################################################################


class Calendar:
    """Class to manage designation of payment dates as holidays according to
    a regional or country-specific calendar convention specified by the user.
    It also supplies an adjustment method which takes in an adjustment
    convention and then applies that to any date that falls on a holiday in the
    specified calendar."""

    def __init__(self, cal_type: CalendarTypes | list | tuple):
        """Create a calendar based on a specified calendar type.

        A list or tuple of calendar types creates a joint calendar whose
        holidays are the union of the holidays of the constituent calendars,
        so a date is a business day only if it is a business day in every
        one of them. This is needed for cross-market products, e.g. a swap
        with a USD SOFR leg and a EUR fixed leg would use
        Calendar([CalendarTypes.US_GOVERNMENT_SECURITIES, CalendarTypes.TARGET])."""

        if isinstance(cal_type, (list, tuple)):

            cal_types = []
            for ct in cal_type:
                if not isinstance(ct, CalendarTypes):
                    raise FinError("Invalid calendar type " + str(ct))
                if ct not in cal_types:
                    cal_types.append(ct)

            if len(cal_types) == 0:
                raise FinError("Calendar type list is empty")

            if len(cal_types) == 1:
                # A single-entry list is just that calendar
                cal_type = cal_types[0]
            else:
                # Joint calendar: cal_type holds the constituent types as a
                # tuple so it never compares equal to a CalendarTypes member
                self.cal_type = tuple(cal_types)
                self._joint_cals = [Calendar(ct) for ct in cal_types]
                return

        if cal_type not in CalendarTypes:
            raise FinError("Invalid calendar type " + str(cal_type))

        self.cal_type = cal_type

    ####################################################################################

    def _validate_year(self, dt):
        if not MIN_CALENDAR_YEAR <= dt.y <= MAX_CALENDAR_YEAR:
            raise FinError(
                f"Calendar supports years "
                f"{MIN_CALENDAR_YEAR}-{MAX_CALENDAR_YEAR}; got {dt.y}"
            )

    ####################################################################################

    @lru_cache(maxsize=None)
    def adjust_cached(self, dt: "Date", bd_type: "BusDayAdjustTypes") -> "Date":
        return self.adjust(dt, bd_type)

    ####################################################################################

    def adjust(self, dt: Date, bd_type: BusDayAdjustTypes):
        """Adjust a payment date if it falls on a holiday according to the
        specified business day convention."""

        if isinstance(bd_type, BusDayAdjustTypes) is False:
            raise FinError("Invalid type passed. Need Finbd_type")

        # If calendar type is NONE then every day is a business day
        if self.cal_type == CalendarTypes.NONE:
            return dt

        if bd_type == BusDayAdjustTypes.NONE:
            return dt

        if bd_type == BusDayAdjustTypes.FOLLOWING:

            # step forward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(1)

            return dt

        if bd_type == BusDayAdjustTypes.MODIFIED_FOLLOWING:

            d_start = dt.d
            m_start = dt.m
            y_start = dt.y

            # step forward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(1)

            # if the business day is in a different month look back
            # for previous first business day one day at a time
            # TODO: I could speed this up by starting it at initial date
            if dt.m != m_start:
                dt = Date(d_start, m_start, y_start)
                while self.is_business_day(dt) is False:
                    dt = dt.add_days(-1)

            return dt

        if bd_type == BusDayAdjustTypes.PRECEDING:

            # if the business day is in the next month look back
            # for previous first business day one day at a time
            while self.is_business_day(dt) is False:
                dt = dt.add_days(-1)

            return dt

        if bd_type == BusDayAdjustTypes.MODIFIED_PRECEDING:

            d_start = dt.d
            m_start = dt.m
            y_start = dt.y

            # step backward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(-1)

            # if the business day is in a different month look forward
            # for previous first business day one day at a time
            # I could speed this up by starting it at initial date
            if dt.m != m_start:
                dt = Date(d_start, m_start, y_start)
                while self.is_business_day(dt) is False:
                    dt = dt.add_days(+1)

            return dt

        raise FinError("Unknown adjustment convention" + str(bd_type))

    ####################################################################################

    def fast_adjust(self, dt: Date, bd_type: BusDayAdjustTypes):
        """Fast adjust a payment date using business day conventions."""

        if not isinstance(bd_type, BusDayAdjustTypes):
            raise FinError("Invalid type passed. Need BusDayAdjustTypes")

        # If no calendar or no adjustment, nothing to do
        if self.cal_type == CalendarTypes.NONE or bd_type == BusDayAdjustTypes.NONE:
            return dt

        # FOLLOWING convention
        if bd_type == BusDayAdjustTypes.FOLLOWING:
            if dt.is_weekend():
                # jump directly to Monday
                if dt.weekday == Date.SAT:
                    dt = dt.add_days(2)
                elif dt.weekday == Date.SUN:
                    dt = dt.add_days(1)
            # if still a holiday (rare), walk forward
            while not self.is_business_day(dt):
                dt = dt.add_days(1)
            return dt

        # MODIFIED FOLLOWING convention
        if bd_type == BusDayAdjustTypes.MODIFIED_FOLLOWING:
            m_start = dt.m
            orig_dt = dt

            if dt.is_weekend():
                if dt.weekday == Date.SAT:
                    dt = dt.add_days(2)
                elif dt.weekday == Date.SUN:
                    dt = dt.add_days(1)
            while not self.is_business_day(dt):
                dt = dt.add_days(1)

            # if moved into a different month → go backwards
            if dt.m != m_start:
                dt = orig_dt
                if dt.is_weekend():
                    if dt.weekday == Date.SAT:
                        dt = dt.add_days(-1)
                    elif dt.weekday == Date.SUN:
                        dt = dt.add_days(-2)
                while not self.is_business_day(dt):
                    dt = dt.add_days(-1)
            return dt

        # PRECEDING convention
        if bd_type == BusDayAdjustTypes.PRECEDING:
            if dt.is_weekend():
                if dt.weekday == Date.SAT:
                    dt = dt.add_days(-1)
                elif dt.weekday == Date.SUN:
                    dt = dt.add_days(-2)
            while not self.is_business_day(dt):
                dt = dt.add_days(-1)
            return dt

        # MODIFIED PRECEDING convention
        if bd_type == BusDayAdjustTypes.MODIFIED_PRECEDING:
            m_start = dt.m
            orig_dt = dt

            if dt.is_weekend():
                if dt.weekday == Date.SAT:
                    dt = dt.add_days(-1)
                elif dt.weekday == Date.SUN:
                    dt = dt.add_days(-2)
            while not self.is_business_day(dt):
                dt = dt.add_days(-1)

            # if moved into a different month → go forward
            if dt.m != m_start:
                dt = orig_dt
                if dt.is_weekend():
                    if dt.weekday == Date.SAT:
                        dt = dt.add_days(2)
                    elif dt.weekday == Date.SUN:
                        dt = dt.add_days(1)
                while not self.is_business_day(dt):
                    dt = dt.add_days(1)
            return dt

        raise FinError("Unknown adjustment convention: " + str(bd_type))

    ####################################################################################

    def add_business_days(self, start_dt: Date, num_days: int):
        """Return the date num_days business days from start_dt.

        The start date itself is not counted.
        """

        if not isinstance(num_days, int):
            raise FinError("Num days must be an integer")

        if num_days == 0:
            return start_dt

        if self.cal_type == CalendarTypes.NONE:
            return start_dt.add_days(num_days)

        step = 1 if num_days > 0 else -1
        days_left = abs(num_days)
        dt = start_dt

        # Safe fast path for weekend-only calendar
        if self.cal_type == CalendarTypes.WEEKEND:
            weeks, days_left = divmod(days_left, 5)
            dt = dt.add_days(step * 7 * weeks)

            while days_left > 0:
                dt = dt.add_days(step)
                if self.is_business_day(dt):
                    days_left -= 1

            return dt

        # General calendar path: correct for holidays
        while days_left > 0:
            dt = dt.add_days(step)
            if self.is_business_day(dt):
                days_left -= 1

        return dt

    ###########################################################################

    def add_business_days_old(self, start_dt: Date, num_days: int):
        """Returns a new date that is num_days business days after Date.
        All holidays in the chosen calendar are assumed not business days."""

        # TODO: REMOVE DATETIME DEPENDENCE HERE ???

        if isinstance(num_days, int) is False:
            raise FinError("Num days must be an integer")

        dt = datetime.date(start_dt.y, start_dt.m, start_dt.d)
        d = dt.day
        m = dt.month
        y = dt.year
        new_dt = Date(d, m, y)

        s = +1
        if num_days < 0:
            num_days = -1 * num_days
            s = -1

        while num_days > 0:
            dt = dt + s * datetime.timedelta(days=1)
            d = dt.day
            m = dt.month
            y = dt.year
            new_dt = Date(d, m, y)

            if self.is_business_day(new_dt) is True:
                num_days -= 1

        return new_dt

    ####################################################################################

    def is_business_day(self, dt: Date):
        """Determines if a date is a business day according to the specified
        calendar. If it is it returns True, otherwise False."""

        self._validate_year(dt)

        # Every day is a business day when the Calendar is NONE
        if self.cal_type == CalendarTypes.NONE:
            return True

        # For all calendars so far, SAT and SUN are not business days
        # If this ever changes I will need to add a filter here.
        if dt.is_weekend():
            return False

        if self.is_holiday(dt) is True:
            return False

        return True

    ####################################################################################

    def is_holiday(self, dt: Date):
        """Determines if a date is a Holiday according to the specified
        calendar. Weekends are not holidays unless the holiday falls on a
        weekend date."""

        self._validate_year(dt)

        if isinstance(self.cal_type, tuple):
            # Joint calendar: a holiday in any constituent calendar
            for cal in self._joint_cals:
                if cal.is_holiday(dt):
                    return True
            return False

        start_dt = Date(1, 1, dt.y)
        day_in_year = dt.excel_dt - start_dt.excel_dt + 1
        weekday = dt.weekday

        if self.cal_type == CalendarTypes.NONE:
            return self.holiday_none(dt)

        if self.cal_type == CalendarTypes.WEEKEND:
            return self.holiday_weekend(dt)

        if self.cal_type == CalendarTypes.AUSTRALIA:
            return self.holiday_australia(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.AUSTRALIA_RITS:
            return self.holiday_australia_rits(dt, weekday)

        if self.cal_type == CalendarTypes.AUD_FX_SETTLEMENT:
            return self.holiday_aud_fx_settlement(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.CANADA:
            return self.holiday_canada(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.FRANCE:
            return self.holiday_france(dt, day_in_year)

        if self.cal_type == CalendarTypes.GERMANY:
            return self.holiday_germany(dt, day_in_year)

        if self.cal_type == CalendarTypes.HONG_KONG:
            return self.holiday_hong_kong(dt, weekday)

        if self.cal_type == CalendarTypes.ITALY:
            return self.holiday_italy(dt, day_in_year)

        if self.cal_type == CalendarTypes.JAPAN:
            return self.holiday_japan(dt, weekday)

        if self.cal_type == CalendarTypes.LONDON:
            return self.holiday_london(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.NEW_YORK:
            return self.holiday_new_york(dt, weekday)

        if self.cal_type == CalendarTypes.NEW_ZEALAND:
            return self.holiday_new_zealand(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.NORWAY:
            return self.holiday_norway(dt, day_in_year)

        if self.cal_type == CalendarTypes.SWEDEN:
            return self.holiday_sweden(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.SWITZERLAND:
            return self.holiday_switzerland(dt, day_in_year)

        if self.cal_type == CalendarTypes.SYDNEY:
            return self.holiday_sydney(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.SINGAPORE:
            return self.holiday_singapore(dt, weekday)

        if self.cal_type == CalendarTypes.TARGET:
            return self.holiday_target(dt, day_in_year)

        if self.cal_type == CalendarTypes.TOKYO:
            return self.holiday_tokyo(dt, weekday)

        if self.cal_type == CalendarTypes.TORONTO:
            return self.holiday_toronto(dt, weekday)

        if self.cal_type == CalendarTypes.UNITED_KINGDOM:
            return self.holiday_united_kingdom(dt, day_in_year, weekday)

        if self.cal_type == CalendarTypes.UNITED_STATES:
            return self.holiday_united_states(dt, weekday)

        if self.cal_type == CalendarTypes.US_FEDERAL_RESERVE:
            return self.holiday_us_federal_reserve(dt, weekday)

        if self.cal_type == CalendarTypes.US_GOVERNMENT_SECURITIES:
            return self.holiday_us_government_securities(dt, weekday)

        if self.cal_type == CalendarTypes.ZURICH:
            return self.holiday_zurich(dt, weekday)

        print(self.cal_type)
        raise FinError("Unknown calendar:" + str(self.cal_type))

    ####################################################################################

    def holiday_weekend(self, dt: Date):
        """Weekends by themselves are a holiday."""

        if dt.is_weekend():
            return True

        return False

    ####################################################################################

    def holiday_australia(self, dt: Date, day_in_year: int, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 26:  # Australia day
            return True

        if m == 1 and d == 27 and weekday == Date.MON:  # Australia day
            return True

        if m == 1 and d == 28 and weekday == Date.MON:  # Australia day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # Australia day
            return True

        if m == 4 and d == 26 and weekday == Date.MON:  # Australia day
            return True

        if m == 6 and d > 7 and d < 15 and weekday == Date.MON:  # Queen
            return True

        # REMOVED AS BEING SYDNEY ONLY
        #if m == 8 and d < 8 and weekday == Date.MON:  # BANK holiday
        #    return True

        #if m == 10 and d < 8 and weekday == Date.MON:  # BANK holiday
        #    return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26 and weekday == Date.MON:  # Xmas
            return True

        if m == 12 and d == 27 and weekday == Date.MON:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday == Date.MON:  # Boxing
            return True

        if m == 12 and d == 28 and weekday == Date.MON:  # Boxing
            return True

        return False

    ###########################################################################

    def holiday_sydney(self, dt: Date, day_in_year: int, weekday: int):
        """
        Sydney / NSW banking calendar.
    
        Weekends are handled separately by is_business_day().
    
        Includes NSW public holidays relevant to Sydney banking,
        plus the NSW Bank Holiday on the first Monday in August.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # New Year's Day
        # --------------------------------------------------------------
    
        if m == 1 and d == 1:
            return True
    
        # Additional day when Jan 1 falls on weekend
        if m == 1 and d == 2 and weekday == Date.MON:
            return True
    
        if m == 1 and d == 3 and weekday == Date.MON:
            return True
    
        # --------------------------------------------------------------
        # Australia Day
        # Jan 26; when weekend -> following Monday
        # --------------------------------------------------------------
    
        if m == 1 and d == 26:
            return True
    
        if (
            m == 1
            and d == 27
            and weekday == Date.MON
        ):
            return True
    
        if (
            m == 1
            and d == 28
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Easter
        # --------------------------------------------------------------
    
        em = self.easter_monday(y)
    
        if dt == em.add_days(-3):  # Good Friday
            return True
    
        if dt == em:  # Easter Monday
            return True
    
        # --------------------------------------------------------------
        # ANZAC Day
        #
        # NSW normally observes 25 April itself.
        # There was an additional NSW holiday in 2026 and 2027.
        # --------------------------------------------------------------
    
        if m == 4 and d == 25:
            return True
    
        if y == 2026 and m == 4 and d == 27:
            return True
    
        if y == 2027 and m == 4 and d == 26:
            return True
    
        # --------------------------------------------------------------
        # King's Birthday
        # Second Monday in June
        # --------------------------------------------------------------
    
        if (
            m == 6
            and 8 <= d <= 14
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # NSW Bank Holiday
        # First Monday in August
        #
        # This is not a general NSW public holiday, but bank branches
        # and certain financial institutions are closed.
        # Relevant for a Sydney banking calendar.
        # --------------------------------------------------------------
    
        if (
            m == 8
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Labour Day
        # First Monday in October
        # --------------------------------------------------------------
    
        if (
            m == 10
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Christmas Day
        # --------------------------------------------------------------
    
        if m == 12 and d == 25:
            return True
    
        # Additional day when Christmas falls on weekend
        if (
            m == 12
            and d == 27
            and weekday in (Date.MON, Date.TUE)
        ):
            return True
    
        # --------------------------------------------------------------
        # Boxing Day
        # --------------------------------------------------------------
    
        if m == 12 and d == 26:
            return True
    
        # Additional Boxing Day holiday
        if (
            m == 12
            and d == 28
            and weekday in (Date.MON, Date.TUE)
        ):
            return True
    
        return False

    ###########################################################################

    def holiday_australia_rits(self, dt: Date, weekday: int):
        """
        RITS settlement calendar.
    
        RITS is open when banks are generally open in either Sydney
        or Melbourne. It is closed on weekends and on holidays
        observed in both NSW and Victoria.
    
        Weekends are handled separately by is_business_day().
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # New Year's Day
        # --------------------------------------------------------------
    
        if m == 1 and d == 1:
            return True
    
        if m == 1 and d == 2 and weekday == Date.MON:
            return True
    
        if m == 1 and d == 3 and weekday == Date.MON:
            return True
    
        # --------------------------------------------------------------
        # Australia Day
        # --------------------------------------------------------------
    
        if m == 1 and d == 26:
            return True
    
        if m == 1 and d == 27 and weekday == Date.MON:
            return True
    
        if m == 1 and d == 28 and weekday == Date.MON:
            return True
    
        # --------------------------------------------------------------
        # Easter
        # --------------------------------------------------------------
    
        em = self.easter_monday(y)
    
        if dt == em.add_days(-3):  # Good Friday
            return True
    
        if dt == em:  # Easter Monday
            return True
    
        # --------------------------------------------------------------
        # ANZAC Day
        #
        # Core RITS closure is Apr 25 itself.
        # Sydney-only additional holidays must NOT automatically
        # close RITS.
        # --------------------------------------------------------------
    
        if m == 4 and d == 25:
            return True
    
        # --------------------------------------------------------------
        # King's Birthday
        # Second Monday in June
        # --------------------------------------------------------------
    
        if (
            m == 6
            and 8 <= d <= 14
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Christmas
        # --------------------------------------------------------------
    
        if m == 12 and d == 25:
            return True
    
        if (
            m == 12
            and d == 27
            and weekday in (Date.MON, Date.TUE)
        ):
            return True
    
        # --------------------------------------------------------------
        # Boxing Day
        # --------------------------------------------------------------
    
        if m == 12 and d == 26:
            return True
    
        if (
            m == 12
            and d == 28
            and weekday in (Date.MON, Date.TUE)
        ):
            return True
    
        return False

    ###########################################################################

    def holiday_aud_fx_settlement(
        self,
        dt: Date,
        day_in_year: int,
        weekday: int,
    ):
        """
        AUD FX settlement calendar.
    
        An AUD FX value date must be a Sydney banking business day.
        Weekends are handled separately by is_business_day().
        """
    
        return self.holiday_sydney(
            dt,
            day_in_year,
            weekday,
        )

    ###########################################################################

    def holiday_united_kingdom(self, dt: Date, day_in_year: int, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        return self.holiday_london(dt, day_in_year, weekday)

    ###########################################################################

    def holiday_london(self, dt: Date, day_in_year: int, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        # One-off or moved bank holidays.
        special_holidays = {
                (2020, 5, 8),   # Early May bank holiday moved for VE Day
                (2022, 6, 2),   # Spring bank holiday moved from 30 May
                (2022, 6, 3),   # Platinum Jubilee
                (2022, 9, 19),  # State funeral of Queen Elizabeth II
                (2023, 5, 8),   # Coronation of King Charles III
            }

        if (y, m, d) in special_holidays:
            return True

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # new years day
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # new years day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 5 and d <= 7 and weekday == Date.MON:
            return True

        if m == 5 and d >= 25 and weekday == Date.MON:
            return True

        if m == 6 and d == 2 and y == 2022:  # SPRING BANK HOLIDAY
            return True

        if m == 6 and d == 3 and y == 2022:  # QUEEN PLAT JUB
            return True

        # Summer bank holiday: last Monday in August.
        if m == 8 and 25 <= d <= 31 and weekday == Date.MON:
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday in (Date.MON, Date.TUE):  # Xmas
            return True

        if m == 12 and d == 28 and weekday in (Date.MON, Date.TUE):  # Xmas
            return True

        return False

    ###########################################################################

    def holiday_france(self, dt: Date, day_in_year: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 5 and d == 1:  # LABOUR DAY
            return True

        if m == 5 and d == 8:  # VICTORY DAY
            return True

        if dt == em.add_days(39 - 1):  # Ascension
            return True

        if dt == em.add_days(50 - 1):  # pentecost
            return True

        if m == 7 and d == 14:  # BASTILLE DAY
            return True

        if m == 8 and d == 15:  # ASSUMPTION
            return True

        if m == 11 and d == 1:  # ALL SAINTS
            return True

        if m == 11 and d == 11:  # ARMISTICE
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        return False

    ###########################################################################

    def holiday_singapore(self, dt: Date, weekday: int):
        """
        Singapore banking / public holiday calendar.
    
        Intended for:
          - SGD derivatives
          - Singapore banking business days
          - SGD settlement
    
        Weekends are handled separately by is_business_day().
    
        Moving religious holidays are maintained from the official
        Singapore Ministry of Manpower annual holiday calendar.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # Fixed-date holidays
        # --------------------------------------------------------------
    
        # New Year's Day
        if m == 1 and d == 1:
            return True
    
        # Labour Day
        if m == 5 and d == 1:
            return True
    
        # National Day
        if m == 8 and d == 9:
            return True
    
        # Christmas Day
        if m == 12 and d == 25:
            return True
    
        # --------------------------------------------------------------
        # Good Friday
        # --------------------------------------------------------------
    
        if dt == self.easter_monday(y).add_days(-3):
            return True
    
        # --------------------------------------------------------------
        # Annual moving holidays
        # --------------------------------------------------------------
    
        special_holidays = {
            2025: {
                (1, 29), (1, 30),   # Chinese New Year
                (3, 31),            # Hari Raya Puasa
                (5, 12),            # Vesak Day
                (6, 7),             # Hari Raya Haji
                (10, 20),           # Deepavali
            },
    
            2026: {
                (2, 17), (2, 18),   # Chinese New Year
                (3, 21),            # Hari Raya Puasa
                (5, 27),            # Hari Raya Haji
                (5, 31),            # Vesak Day
                (6, 1),             # Vesak observed
                (8, 10),            # National Day observed
                (11, 8),            # Deepavali
                (11, 9),            # Deepavali observed
            },
    
            2027: {
                (2, 6), (2, 7),     # Chinese New Year
                (2, 8),             # CNY observed
                (3, 10),            # Hari Raya Puasa
                (5, 17),            # Hari Raya Haji
                (5, 20),            # Vesak Day
                (10, 28),           # Deepavali
            },
        }
    
        if y not in special_holidays:
            pass
#            _warn_incomplete_calendar("Hong Kong", y)
        else:
            if (m, d) in special_holidays[y]:
                return True

    ###########################################################################

    def holiday_sweden(self, dt: Date, day_in_year: int, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if dt == em.add_days(39 - 1):  # Ascension
            return True

        if m == 5 and d == 1:  # labour day
            return True

        if m == 6 and d == 6:  # June
            return True

        if m == 6 and d > 18 and d < 26 and weekday == Date.FRI:  # Midsummer
            return True

        if m == 12 and d == 24:  # Xmas eve
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 31:  # NYE
            return True

        return False

    ###########################################################################

    def holiday_germany(self, dt: Date, day_in_year: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 5 and d == 1:  # LABOUR DAY
            return True

        if dt == em.add_days(39 - 1):  # Ascension
            return True

        if dt == em.add_days(50 - 1):  # pentecost
            return True

        if m == 10 and d == 3:  # GERMAN UNITY DAY
            return True

        if m == 12 and d == 24:  # Xmas eve
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        return False

    ###########################################################################

    def holiday_switzerland(self, dt: Date, day_in_year: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2:  # berchtoldstag
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if dt == em.add_days(39 - 1):  # Ascension
            return True

        if dt == em.add_days(50 - 1):  # pentecost
            return True

        if m == 5 and d == 1:  # Labour day
            return True

        if m == 8 and d == 1:  # National day
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        return False

    ###########################################################################

    def holiday_zurich(self, dt: Date, weekday: int):
        """
        Zurich / CHF currency holiday calendar.
    
        Intended for:
          - SARON
          - CHF OIS / swaps
          - CHF money-market settlement
    
        Weekends are handled separately by is_business_day().
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # New Year's Day
        if m == 1 and d == 1:
            return True
    
        # Berchtold's Day
        if m == 1 and d == 2:
            return True
    
        # Good Friday
        if dt == self.easter_monday(y).add_days(-3):
            return True
    
        # Easter Monday
        if dt == self.easter_monday(y):
            return True
    
        # Labour Day
        if m == 5 and d == 1:
            return True
    
        # Ascension Day
        if dt == self.easter_monday(y).add_days(38):
            return True
    
        # Whit Monday / Pentecost Monday
        if dt == self.easter_monday(y).add_days(49):
            return True
    
        # Swiss National Day
        if m == 8 and d == 1:
            return True
    
        # Christmas Day
        if m == 12 and d == 25:
            return True
    
        # St Stephen's Day
        if m == 12 and d == 26:
            return True
    
        return False

    ###########################################################################

    @lru_cache(maxsize=None)
    def _japan_holiday_set(self, y: int):
    
        if y < 1980 or y > 2099:
            raise FinError(
                f"Japan calendar supports years 1980-2099; got {y}"
            )
    
        national = set()
    
        def add(m, d):
            national.add((m, d))
    
        def nth_monday(month, n):
            first = datetime.date(y, month, 1)
            offset = (0 - first.weekday()) % 7
            return 1 + offset + 7 * (n - 1)
    
        def vernal_equinox_day():
            return int(
                20.8431
                + 0.242194 * (y - 1980)
                - int((y - 1980) / 4)
            )
    
        def autumnal_equinox_day():
            return int(
                23.2488
                + 0.242194 * (y - 1980)
                - int((y - 1980) / 4)
            )
    
        # --------------------------------------------------------------
        # National holidays
        # --------------------------------------------------------------
    
        add(1, 1)
    
        # Coming of Age Day
        if y >= 2000:
            add(1, nth_monday(1, 2))
        else:
            add(1, 15)
    
        # National Foundation Day
        add(2, 11)
    
        # Emperor's Birthday
        if 1989 <= y <= 2018:
            add(12, 23)
        elif y >= 2020:
            add(2, 23)
    
        # Vernal Equinox
        add(3, vernal_equinox_day())
    
        # Showa / Greenery / former Emperor's Birthday
        add(4, 29)
    
        # Golden Week
        add(5, 3)
    
        if y >= 2007:
            add(5, 4)
    
        add(5, 5)
    
        # Marine Day
        if y == 2020:
            add(7, 23)
        elif y == 2021:
            add(7, 22)
        elif y >= 2003:
            add(7, nth_monday(7, 3))
        elif y >= 1996:
            add(7, 20)
    
        # Sports Day special moves
        if y == 2020:
            add(7, 24)
        elif y == 2021:
            add(7, 23)
    
        # Mountain Day
        if y == 2020:
            add(8, 10)
        elif y == 2021:
            add(8, 8)
        elif y >= 2016:
            add(8, 11)
    
        # Respect for the Aged Day
        if y >= 2003:
            add(9, nth_monday(9, 3))
        else:
            add(9, 15)
    
        # Autumnal Equinox
        add(9, autumnal_equinox_day())
    
        # Sports / Health and Sports Day
        if y not in (2020, 2021):
            if y >= 2000:
                add(10, nth_monday(10, 2))
            else:
                add(10, 10)
    
        # Culture Day
        add(11, 3)
    
        # Labour Thanksgiving Day
        add(11, 23)
    
        # --------------------------------------------------------------
        # 2019 Imperial one-offs
        # --------------------------------------------------------------
    
        if y == 2019:
            add(5, 1)
            add(10, 22)
    
        # --------------------------------------------------------------
        # Citizens' holidays
        #
        # A weekday between two national holidays becomes a holiday.
        # Example: 22-Sep-2026.
        # --------------------------------------------------------------
    
        day = datetime.date(y, 1, 2)
        end = datetime.date(y, 12, 31)
    
        citizens = set()
    
        while day < end:
            key = (day.month, day.day)
    
            if key not in national:
                prev_day = day - datetime.timedelta(days=1)
                next_day = day + datetime.timedelta(days=1)
    
                if (
                    (prev_day.month, prev_day.day) in national
                    and (next_day.month, next_day.day) in national
                ):
                    citizens.add(key)
    
            day += datetime.timedelta(days=1)
    
        national.update(citizens)
    
        # --------------------------------------------------------------
        # Substitute holidays
        #
        # For a Sunday national holiday, move forward to the first
        # non-national-holiday date.
        #
        # Example: 3-May-2026 -> 6-May-2026.
        # --------------------------------------------------------------
    
        original_holidays = sorted(
            datetime.date(y, m, d)
            for m, d in national
        )
    
        substitutes = set()
    
        for holiday in original_holidays:
            if holiday.weekday() != 6:  # Sunday
                continue
    
            substitute = holiday + datetime.timedelta(days=1)
    
            while (
                substitute.year == y
                and (
                    (substitute.month, substitute.day) in national
                    or (substitute.month, substitute.day) in substitutes
                )
            ):
                substitute += datetime.timedelta(days=1)
    
            if substitute.year == y:
                substitutes.add(
                    (substitute.month, substitute.day)
                )
    
        national.update(substitutes)
    
        # --------------------------------------------------------------
        # BOJ-specific closures
        # --------------------------------------------------------------
    
        return frozenset(national)

    
    def holiday_japan(self, dt: Date, weekday: int):
        """Japanese statutory public holidays. Weekends handled separately."""
        return (dt.m, dt.d) in self._japan_holiday_set(dt.y)

    ###########################################################################

    @lru_cache(maxsize=None)
    def _tokyo_holiday_set(self, y: int):
    
        tokyo = set(self._japan_holiday_set(y))
    
        # BOJ-specific year-end closures
        tokyo.add((1, 2))
        tokyo.add((1, 3))
        tokyo.add((12, 31))
    
        return frozenset(tokyo)
    
    
    def holiday_tokyo(self, dt: Date, weekday: int):
        return (dt.m, dt.d) in self._tokyo_holiday_set(dt.y)

    ###########################################################################

    def holiday_new_zealand(self, dt: Date, day_in_year: int, weekday: int):
        """
        New Zealand public / banking holiday calendar.
    
        Weekends are handled separately by is_business_day().
    
        Includes:
          - New Year's Day and day after New Year's Day
          - Waitangi Day
          - Good Friday
          - Easter Monday
          - ANZAC Day
          - King's / Queen's Birthday
          - Matariki
          - Labour Day
          - Christmas Day
          - Boxing Day
    
        Provincial anniversary days are deliberately excluded.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # New Year
        #
        # Jan 1 and Jan 2 are holidays.
        # If they fall on Sat/Sun, observed on Mon/Tue.
        # --------------------------------------------------------------
    
        if m == 1 and d in (1, 2):
            return True
    
        # Jan 1 Saturday -> Jan 3 Monday
        # Jan 2 Sunday   -> Jan 4 Tuesday
        if m == 1 and d == 3 and weekday == Date.MON:
            return True
    
        if m == 1 and d == 4 and weekday in (Date.MON, Date.TUE):
            return True
    
        # --------------------------------------------------------------
        # Waitangi Day
        # February 6
        #
        # Since 2014, Mondayised when falling on weekend.
        # --------------------------------------------------------------
    
        if m == 2 and d == 6:
            return True
    
        if y >= 2014:
            if m == 2 and d == 7 and weekday == Date.MON:
                return True
    
            if m == 2 and d == 8 and weekday == Date.MON:
                return True
    
        # --------------------------------------------------------------
        # Easter
        # --------------------------------------------------------------
    
        em = self.easter_monday(y)
    
        if dt == em.add_days(-3):  # Good Friday
            return True
    
        if dt == em:  # Easter Monday
            return True
    
        # --------------------------------------------------------------
        # ANZAC Day
        # April 25
        #
        # Since 2014, Mondayised when falling on weekend.
        # --------------------------------------------------------------
    
        if m == 4 and d == 25:
            return True
    
        if y >= 2014:
            if m == 4 and d == 26 and weekday == Date.MON:
                return True
    
            if m == 4 and d == 27 and weekday == Date.MON:
                return True
    
        # --------------------------------------------------------------
        # King's / Queen's Birthday
        # First Monday in June
        # --------------------------------------------------------------
    
        if (
            m == 6
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Matariki
        #
        # Dates are prescribed annually and do not follow a simple
        # Gregorian rule, so use an explicit table.
        # --------------------------------------------------------------
    
        matariki = {
            2022: (6, 24),
            2023: (7, 14),
            2024: (6, 28),
            2025: (6, 20),
            2026: (7, 10),
            2027: (6, 25),
            2028: (7, 14),
            2029: (7, 6),
            2030: (6, 21),
            2031: (7, 11),
            2032: (7, 2),
            2033: (6, 24),
            2034: (7, 7),
            2035: (6, 29),
            2036: (7, 18),
            2037: (7, 10),
            2038: (6, 25),
            2039: (7, 15),
            2040: (7, 6),
            2041: (7, 19),
            2042: (7, 11),
            2043: (7, 3),
            2044: (6, 24),
            2045: (7, 7),
            2046: (6, 29),
            2047: (7, 19),
            2048: (7, 3),
            2049: (6, 25),
            2050: (7, 15),
            2051: (6, 30),
            2052: (6, 21),
        }
    
        if y in matariki:
            mm, dd = matariki[y]
        
            if m == mm and d == dd:
                return True


        # --------------------------------------------------------------
        # Labour Day
        # Fourth Monday in October
        # --------------------------------------------------------------
    
        if (
            m == 10
            and 22 <= d <= 28
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Christmas / Boxing Day
        #
        # Dec 25 and Dec 26.
        # Weekend holidays are observed on Mon/Tue.
        # --------------------------------------------------------------
    
        if m == 12 and d in (25, 26):
            return True
    
        if m == 12 and d == 27 and weekday in (Date.MON, Date.TUE):
            return True
    
        if m == 12 and d == 28 and weekday in (Date.MON, Date.TUE):
            return True
    
        return False

    ###########################################################################

    def holiday_norway(self, dt: Date, day_in_year: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-4):  # holy thursday
            return True

        if dt == em.add_days(-3):  # good friday
            return True

        if dt == em:  # Easter Monday
            return True

        if dt == em.add_days(38):  # Ascension
            return True

        if dt == em.add_days(49):  # Pentecost
            return True

        if m == 5 and d == 1:  # May day
            return True

        if m == 5 and d == 17:  # Independence day
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        return False

    ###########################################################################

    def holiday_united_states(self, dt: Date, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday.
        This is a legacy generic US calendar that contains the superset of
        holidays for bond markets, NYSE, and public holidays. For each of
        these and other categories there will be some variations."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # NYD
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # NYD
            return True

        if m == 1 and d >= 15 and d < 22 and weekday == Date.MON:  # MLK
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == Date.MON:  # GW
            return True

        if m == 5 and d >= 25 and d <= 31 and weekday == Date.MON:  # MD
            return True

        # Juneteenth
        if y >= 2021:
            if m == 6 and d == 19:
                return True
            if m == 6 and d == 18 and weekday == Date.FRI:
                return True
            if m == 6 and d == 20 and weekday == Date.MON:
                return True

        if m == 7 and d == 4:  # Indep day
            return True

        if m == 7 and d == 5 and weekday == Date.MON:  # Indep day
            return True

        if m == 7 and d == 3 and weekday == Date.FRI:  # Indep day
            return True

        if m == 9 and d >= 1 and d < 8 and weekday == Date.MON:  # Lab
            return True

        if m == 10 and d >= 8 and d < 15 and weekday == Date.MON:  # CD
            return True

        if m == 11 and d == 11:  # Veterans day
            return True

        if m == 11 and d == 12 and weekday == Date.MON:  # Vets
            return True

        if m == 11 and d == 10 and weekday == Date.FRI:  # Vets
            return True

        if m == 11 and d >= 22 and d < 29 and weekday == Date.THU:  # TG
            return True

        if m == 12 and d == 24 and weekday == Date.FRI:  # Xmas holiday
            return True

        if m == 12 and d == 25:  # Xmas holiday
            return True

        if m == 12 and d == 26 and weekday == Date.MON:  # Xmas holiday
            return True

        if m == 12 and d == 31 and weekday == Date.FRI:
            return True

        return False

    ###########################################################################

    def holiday_us_government_securities(self, dt: Date, weekday: int):
        """
        U.S. Government Securities / SIFMA calendar.
    
        Intended for:
          - U.S. Treasury / government securities market
          - SOFR-related business-day logic
    
        Weekends are handled separately by is_business_day().
    
        NOTE:
        SIFMA recommendations can change and exceptional dates occur.
        Keep year-specific overrides for authoritative market closures.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # Explicit full-close / SOFR holiday overrides
        # --------------------------------------------------------------
    
        special_holidays = {
            # 2026
            (2026, 4, 3),   # Good Friday - no SOFR publication
            (2026, 7, 3),   # Independence Day observed / SIFMA close
    
            # Add future exceptional dates from SIFMA / NY Fed here.
        }
    
        if (y, m, d) in special_holidays:
            return True
    
        # --------------------------------------------------------------
        # New Year's Day
        # --------------------------------------------------------------
    
        if m == 1 and d == 1:
            return True
    
        # Friday observation when Jan 1 falls Saturday
        if m == 12 and d == 31 and weekday == Date.FRI:
            return True
    
        # Monday observation when Jan 1 falls Sunday
        if m == 1 and d == 2 and weekday == Date.MON:
            return True
    
        # --------------------------------------------------------------
        # Martin Luther King Jr. Day
        # Third Monday in January
        # --------------------------------------------------------------
    
        if (
            m == 1
            and 15 <= d <= 21
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Presidents Day
        # Third Monday in February
        # --------------------------------------------------------------
    
        if (
            m == 2
            and 15 <= d <= 21
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Good Friday
        #
        # Normally a government-securities / settlement holiday.
        # Compute from Easter Monday:
        # Good Friday = Easter Monday - 3 days.
        #
        # Exceptional years should remain in special_holidays.
        # --------------------------------------------------------------

        if dt == self.easter_monday(y).add_days(-3):
            return True

        # --------------------------------------------------------------
        # Memorial Day
        # Last Monday in May
        # --------------------------------------------------------------
    
        if (
            m == 5
            and 25 <= d <= 31
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Juneteenth
        # --------------------------------------------------------------
    
        if y >= 2022:
            if m == 6 and d == 19:
                return True
    
            # June 19 Saturday -> Friday June 18
            if (
                m == 6
                and d == 18
                and weekday == Date.FRI
            ):
                return True
    
            # June 19 Sunday -> Monday June 20
            if (
                m == 6
                and d == 20
                and weekday == Date.MON
            ):
                return True
    
        # --------------------------------------------------------------
        # Independence Day
        # --------------------------------------------------------------
    
        if m == 7 and d == 4:
            return True
    
        if (
            m == 7
            and d == 3
            and weekday == Date.FRI
        ):
            return True
    
        if (
            m == 7
            and d == 5
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Labor Day
        # First Monday in September
        # --------------------------------------------------------------
    
        if (
            m == 9
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Columbus Day / Indigenous Peoples' Day
        # Second Monday in October
        #
        # U.S. government securities market closes.
        # --------------------------------------------------------------
    
        if (
            m == 10
            and 8 <= d <= 14
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Veterans Day
        # --------------------------------------------------------------
    
        if m == 11 and d == 11:
            return True
    
        if (
            m == 11
            and d == 10
            and weekday == Date.FRI
        ):
            return True
    
        if (
            m == 11
            and d == 12
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Thanksgiving
        # Fourth Thursday in November
        # --------------------------------------------------------------
    
        if (
            m == 11
            and 22 <= d <= 28
            and weekday == Date.THU
        ):
            return True
    
        # --------------------------------------------------------------
        # Christmas
        # --------------------------------------------------------------
    
        if m == 12 and d == 25:
            return True
    
        if (
            m == 12
            and d == 24
            and weekday == Date.FRI
        ):
            return True
    
        if (
            m == 12
            and d == 26
            and weekday == Date.MON
        ):
            return True
    
        return False

    ###########################################################################

    def holiday_us_federal_reserve(self, dt: Date, weekday: int):
        """
        Federal Reserve Bank business-day calendar.
    
        Intended for:
          - Federal Reserve Bank business days
          - Fedwire / National Settlement Service
          - EFFR / OBFR publication calendar
    
        Weekends are handled separately by is_business_day().
    
        Important:
          - Saturday holidays are NOT observed on Friday by Federal Reserve Banks.
          - Sunday holidays ARE observed on Monday.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # New Year's Day
        # Jan 1.
        #
        # If Sunday -> Monday Jan 2.
        # If Saturday -> no Friday observation for Federal Reserve Banks.
        # --------------------------------------------------------------
    
        if m == 1 and d == 1:
            return True
    
        if (
            m == 1
            and d == 2
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Martin Luther King Jr. Day
        # Third Monday in January
        # --------------------------------------------------------------
    
        if (
            m == 1
            and 15 <= d <= 21
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Washington's Birthday
        # Third Monday in February
        # --------------------------------------------------------------
    
        if (
            m == 2
            and 15 <= d <= 21
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Memorial Day
        # Last Monday in May
        # --------------------------------------------------------------
    
        if (
            m == 5
            and 25 <= d <= 31
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Juneteenth
        # Federal holiday since 2021.
        #
        # Sunday -> Monday.
        # Saturday -> no Friday observation for Federal Reserve Banks.
        # --------------------------------------------------------------
    
        if y >= 2021:
            if m == 6 and d == 19:
                return True
    
            if (
                m == 6
                and d == 20
                and weekday == Date.MON
            ):
                return True
    
        # --------------------------------------------------------------
        # Independence Day
        # --------------------------------------------------------------
    
        if m == 7 and d == 4:
            return True
    
        if (
            m == 7
            and d == 5
            and weekday == Date.MON
        ):
            return True
    
        # Do NOT close Friday July 3 when July 4 is Saturday.
    
        # --------------------------------------------------------------
        # Labor Day
        # First Monday in September
        # --------------------------------------------------------------
    
        if (
            m == 9
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Columbus Day
        # Second Monday in October
        # --------------------------------------------------------------
    
        if (
            m == 10
            and 8 <= d <= 14
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Veterans Day
        # Nov 11.
        #
        # Sunday -> Monday.
        # Saturday -> no Friday observation for Federal Reserve Banks.
        # --------------------------------------------------------------
    
        if m == 11 and d == 11:
            return True
    
        if (
            m == 11
            and d == 12
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Thanksgiving
        # Fourth Thursday in November
        # --------------------------------------------------------------
    
        if (
            m == 11
            and 22 <= d <= 28
            and weekday == Date.THU
        ):
            return True
    
        # --------------------------------------------------------------
        # Christmas
        # --------------------------------------------------------------
    
        if m == 12 and d == 25:
            return True
    
        if (
            m == 12
            and d == 26
            and weekday == Date.MON
        ):
            return True
    
        # Again: no Friday Dec 24 observation when Dec 25 is Saturday.
    
        return False

    ###########################################################################

    def holiday_new_york(self, dt: Date, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday.
        This is a NY specific US calendar close to Federal Reserve calendar."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # NYD
            return True

        if m == 1 and d >= 15 and d < 22 and weekday == Date.MON:  # MLK
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == Date.MON:  # GW
            return True


        if m == 5 and d >= 25 and d <= 31 and weekday == Date.MON:  # MD
            return True

        # Juneteenth
        if y >= 2021:
            if m == 6 and d == 19:
                return True
    
            if m == 6 and d == 18 and weekday == Date.FRI:
                return True
    
            if m == 6 and d == 20 and weekday == Date.MON:
                return True

        if m == 7 and d == 4:  # Indep day
            return True

        if m == 7 and d == 5 and weekday == Date.MON:  # Indep day
            return True

        if m == 7 and d == 3 and weekday == Date.FRI:  # Indep day
            return True

        if m == 9 and d >= 1 and d < 8 and weekday == Date.MON:  # Lab
            return True

        if m == 10 and d >= 8 and d < 15 and weekday == Date.MON:  # CD
            return True

        if m == 11 and d == 11:  # Veterans day
            return True

        if m == 11 and d == 12 and weekday == Date.MON:  # Vets
            return True

        if m == 11 and d == 10 and weekday == Date.FRI:  # Vets
            return True

        if m == 11 and d >= 22 and d < 29 and weekday == Date.THU:  # TG
            return True

        if m == 12 and d == 24 and weekday == Date.FRI:  # Xmas holiday
            return True

        if m == 12 and d == 25:  # Xmas holiday
            return True

        if m == 12 and d == 26 and weekday == Date.MON:  # Xmas holiday
            return True

        return False

    ###########################################################################

    def holiday_canada(self, dt: Date, day_in_year: int, weekday: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # NYD
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # NYD
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == Date.MON:  # FAMILY
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # good friday
            return True

        if m == 5 and d >= 18 and d < 25 and weekday == Date.MON:  # VICTORIA
            return True

        if m == 7 and d == 1:  # Canada day
            return True

        if m == 7 and d == 2 and weekday == Date.MON:  # Canada day
            return True

        if m == 7 and d == 3 and weekday == Date.MON:  # Canada day
            return True

        if m == 8 and d < 8 and weekday == Date.MON:  # Provincial
            return True

        if m == 9 and d < 8 and weekday == Date.MON:  # Labor
            return True

        if m == 10 and d >= 8 and d < 15 and weekday == Date.MON:  # THANKS
            return True

        if m == 11 and d == 11:  # Veterans day
            return True

        if m == 11 and d == 12 and weekday == Date.MON:  # Vets
            return True

        if m == 11 and d == 13 and weekday == Date.MON:  # Vets
            return True

        if m == 12 and d == 25:  # Xmas holiday
            return True

        if m == 12 and d == 26 and weekday == Date.MON:  # Xmas holiday
            return True

        if m == 12 and d == 27 and weekday == Date.MON:  # Xmas holiday
            return True

        if m == 12 and d == 26:  # Boxing holiday
            return True

        if m == 12 and d == 27 and weekday == Date.MON:  # Boxing holiday
            return True

        if m == 12 and d == 28 and weekday == Date.TUE:  # Boxing holiday
            return True

        return False

    ###########################################################################

    def holiday_toronto(self, dt: Date, weekday: int):
        """
        Toronto banking calendar.
    
        Intended for:
          - CORRA
          - CAD OIS / swaps
          - Toronto banking business-day conventions
    
        Weekends are handled separately by is_business_day().
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # New Year's Day
        # --------------------------------------------------------------
    
        if m == 1 and d == 1:
            return True
    
        if m == 1 and d == 2 and weekday == Date.MON:
            return True
    
        # --------------------------------------------------------------
        # Family Day
        # Third Monday in February (Ontario)
        # --------------------------------------------------------------
    
        if (
            y >= 2008
            and m == 2
            and 15 <= d <= 21
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Good Friday
        # --------------------------------------------------------------
    
        if dt == self.easter_monday(y).add_days(-3):
            return True
    
        # --------------------------------------------------------------
        # Victoria Day
        #
        # Monday preceding May 25
        # --------------------------------------------------------------
    
        if (
            m == 5
            and 18 <= d <= 24
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Canada Day
        # --------------------------------------------------------------
    
        if m == 7 and d == 1:
            return True
    
        # July 1 Saturday -> Monday July 3
        if (
            m == 7
            and d == 3
            and weekday == Date.MON
        ):
            return True
    
        # July 1 Sunday -> Monday July 2
        if (
            m == 7
            and d == 2
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Civic Holiday
        # First Monday in August
        #
        # Ontario / Toronto banking holiday.
        # --------------------------------------------------------------
    
        if (
            m == 8
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Labour Day
        # First Monday in September
        # --------------------------------------------------------------
    
        if (
            m == 9
            and 1 <= d <= 7
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # National Day for Truth and Reconciliation
        # September 30
        #
        # Include if your Toronto/CORRA convention treats Schedule I
        # banks as closed on this date.
        # --------------------------------------------------------------
    
        if y >= 2021 and m == 9 and d == 30:
            return True
    
        # --------------------------------------------------------------
        # Thanksgiving
        # Second Monday in October
        # --------------------------------------------------------------
    
        if (
            m == 10
            and 8 <= d <= 14
            and weekday == Date.MON
        ):
            return True
    
        # --------------------------------------------------------------
        # Remembrance Day
        # --------------------------------------------------------------
    
        if m == 11 and d == 11:
            return True
    
        # --------------------------------------------------------------
        # Christmas Day
        # --------------------------------------------------------------
    
        if m == 12 and d == 25:
            return True
    
        # --------------------------------------------------------------
        # Boxing Day
        # --------------------------------------------------------------
    
        if m == 12 and d == 26:
            return True
    
        # Observed Christmas / Boxing Day
        if (
            m == 12
            and d in (27, 28)
            and weekday in (Date.MON, Date.TUE)
        ):
            return True
    
        return False

    ###########################################################################

    def holiday_italy(self, dt: Date, day_in_year: int):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # LIBERATION DAY
            return True

        if m == 5 and d == 1:  # LABOUR DAY
            return True

        if m == 6 and d == 2 and y > 1999:  # REPUBLIC DAY
            return True

        if m == 8 and d == 15:  # ASSUMPTION
            return True

        if m == 11 and d == 1:  # ALL SAINTS
            return True

        if m == 12 and d == 8:  # IMMAC CONC
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        return False

    ###########################################################################

    def holiday_hong_kong(self, dt: Date, weekday: int):
        """
        Hong Kong banking / general holiday calendar.
    
        Intended for:
          - HKD derivatives
          - Hong Kong banking business days
          - HKD settlement
    
        Weekends are handled separately by is_business_day().
    
        Lunar and substitution holidays should be maintained
        from the official gazetted Hong Kong holiday calendar.
        """
    
        m = dt.m
        d = dt.d
        y = dt.y
    
        # --------------------------------------------------------------
        # Fixed-date holidays
        # --------------------------------------------------------------
    
        # New Year's Day
        if m == 1 and d == 1:
            return True
    
        # Labour Day
        if m == 5 and d == 1:
            return True
    
        # HKSAR Establishment Day
        if m == 7 and d == 1:
            return True
    
        # National Day
        if m == 10 and d == 1:
            return True
    
        # Christmas Day
        if m == 12 and d == 25:
            return True
    
        # --------------------------------------------------------------
        # Easter holidays
        # --------------------------------------------------------------
    
        em = self.easter_monday(y)
    
        # Good Friday
        if dt == em.add_days(-3):
            return True
    
        # Day following Good Friday
        if dt == em.add_days(-2):
            return True
    
        # Easter Monday
        if dt == em:
            return True
    
        # --------------------------------------------------------------
        # Year-specific Hong Kong holidays
        #
        # Lunar New Year, Ching Ming, Buddha's Birthday,
        # Tuen Ng, Mid-Autumn, Chung Yeung and substitution
        # days are best taken from gazetted annual calendars.
        # --------------------------------------------------------------
    
        special_holidays = {
    
            2026: {
                # Lunar New Year
                (2, 17),
                (2, 18),
                (2, 19),
    
                # Ching Ming substitution / Easter interaction
                (4, 6),
                (4, 7),
    
                # Birthday of the Buddha substitution
                (5, 25),
    
                # Tuen Ng Festival
                (6, 19),
    
                # Day following Mid-Autumn Festival
                (9, 26),
    
                # Chung Yeung substitution
                (10, 19),
    
                # First weekday after Christmas
                (12, 26),
            },
    
            2027: {
                # Lunar New Year
                (2, 6),
                (2, 8),
                (2, 9),
    
                # Ching Ming
                (4, 5),
    
                # Birthday of the Buddha
                (5, 13),
    
                # Tuen Ng Festival
                (6, 9),
    
                # Day following Mid-Autumn Festival
                (9, 16),
    
                # Chung Yeung Festival
                (10, 8),
    
                # First weekday after Christmas
                (12, 27),
            },
        }

        if y not in special_holidays:
            pass
#            _warn_incomplete_calendar("Singapore", y)
        else:
            if (m, d) in special_holidays[y]:
                return True

    ###########################################################################

    def holiday_target(self, dt: Date, day_in_year):
        """Only bank holidays. Weekends by themselves are not a holiday."""

        m = dt.m
        d = dt.d
        y = dt.y

        if m == 1 and d == 1:  # new year's day
            return True

        if m == 5 and d == 1:  # May day
            return True

        em = self.easter_monday(y)

        if dt == em.add_days(-3):  # Good Friday
            return True
        
        if dt == em:  # Easter Monday
            return True

        if m == 12 and d == 25:  # Xmas bank holiday
            return True

        if m == 12 and d == 26:  # Xmas bank holiday
            return True

        return False

    ###########################################################################

    def holiday_none(self, dt: Date = None):
        """No day is a holiday."""
        return False

    ###########################################################################

    def get_holiday_list(self, year: int):
        """Generate a list of declared calendar holidays in a given year.

        This includes declared holidays that fall on weekends, but does not
        include ordinary weekends unless the calendar type is WEEKEND.
        """

        start_dt = Date(1, 1, year)
        end_dt = Date(1, 1, year + 1)

        holiday_list = []

        while start_dt < end_dt:
            if self.is_holiday(start_dt):
                holiday_list.append(str(start_dt))

            start_dt = start_dt.add_days(1)

        return holiday_list

    ###########################################################################

    def easter_monday(self, year: float):
        """Get the day in a given year that is Easter Monday. This is not
        easy to compute, so we rely on a pre-calculated array."""

        if year < 1901 or year > 2100:
            raise FinError(
            f"Easter calculation only supported for years 1901-2100: {year}")

        em_days = easter_monday_day[year - 1901]
        start_dt = Date(1, 1, year)
        em = start_dt.add_days(em_days - 1)
        return em

    ###########################################################################

    def __str__(self):
        if isinstance(self.cal_type, tuple):
            s = "JOINT(" + ",".join(ct.name for ct in self.cal_type) + ")"
            return s
        s = self.cal_type.name
        return s

    ###########################################################################

    def __repr__(self):
        return self.__str__()


########################################################################################
