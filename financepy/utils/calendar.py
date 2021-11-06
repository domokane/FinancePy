###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


###############################################################################
# TODO: Do some timings and tidy up logic in adjustment function
###############################################################################

import datetime
from enum import Enum
from .date import Date
from .error import FinError

# from numba import njit, jit, int64, boolean

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


class BusDayAdjustTypes(Enum):
    NONE = 1
    FOLLOWING = 2
    MODIFIED_FOLLOWING = 3
    PRECEDING = 4
    MODIFIED_PRECEDING = 5


class CalendarTypes(Enum):
    NONE = 1
    WEEKEND = 2
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
    TARGET = 13
    UNITED_STATES = 14
    UNITED_KINGDOM = 15


class DateGenRuleTypes(Enum):
    FORWARD = 1
    BACKWARD = 2

###############################################################################


class Calendar:
    """ Class to manage designation of payment dates as holidays according to
    a regional or country-specific calendar convention specified by the user.
    It also supplies an adjustment method which takes in an adjustment
    convention and then applies that to any date that falls on a holiday in the
    specified calendar. """

    def __init__(self,
                 calendar_type: CalendarTypes):
        """ Create a calendar based on a specified calendar type. """

        if calendar_type not in CalendarTypes:
            raise FinError(
                "Need to pass FinCalendarType and not " +
                str(calendar_type))

        self._type = calendar_type

    ###########################################################################

    def adjust(self,
               dt: Date,
               busDayConventionType: BusDayAdjustTypes):
        """ Adjust a payment date if it falls on a holiday according to the
        specified business day convention. """

        if type(busDayConventionType) != BusDayAdjustTypes:
            raise FinError("Invalid type passed. Need FinBusDayConventionType")

        if busDayConventionType == BusDayAdjustTypes.NONE:
            return dt

        elif busDayConventionType == BusDayAdjustTypes.FOLLOWING:

            # step forward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(1)

            return dt

        elif busDayConventionType == BusDayAdjustTypes.MODIFIED_FOLLOWING:

            d_start = dt._d
            m_start = dt._m
            y_start = dt._y

            # step forward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(1)

            # if the business day is in a different month look back
            # for previous first business day one day at a time
            # TODO: I could speed this up by starting it at initial date
            if dt._m != m_start:
                dt = Date(d_start, m_start, y_start)
                while self.is_business_day(dt) is False:
                    dt = dt.add_days(-1)

            return dt

        elif busDayConventionType == BusDayAdjustTypes.PRECEDING:

            # if the business day is in the next month look back
            # for previous first business day one day at a time
            while self.is_business_day(dt) is False:
                dt = dt.add_days(-1)

            return dt

        elif busDayConventionType == BusDayAdjustTypes.MODIFIED_PRECEDING:

            d_start = dt._d
            m_start = dt._m
            y_start = dt._y

            # step backward until we find a business day
            while self.is_business_day(dt) is False:
                dt = dt.add_days(-1)

            # if the business day is in a different month look forward
            # for previous first business day one day at a time
            # I could speed this up by starting it at initial date
            if dt._m != m_start:
                dt = Date(d_start, m_start, y_start)
                while self.is_business_day(dt) is False:
                    dt = dt.add_days(+1)

            return dt

        else:

            raise FinError("Unknown adjustment convention" +
                           str(busDayConventionType))

        return dt

###############################################################################

    def add_business_days(self,
                          start_date: Date,
                          numDays: int):
        """ Returns a new date that is numDays business days after Date.
        All holidays in the chosen calendar are assumed not business days. """

        # TODO: REMOVE DATETIME DEPENDENCE HERE ???

        if isinstance(numDays, int) is False:
            raise FinError("Num days must be an integer")

        dt = datetime.date(start_date._y, start_date._m, start_date._d)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = Date(d, m, y)

        s = +1
        if numDays < 0:
            numDays = -1 * numDays
            s = -1

        while numDays > 0:
            dt = dt + s * datetime.timedelta(days=1)
            d = dt.day
            m = dt.month
            y = dt.year
            newDt = Date(d, m, y)

            if self.is_business_day(newDt) is True:
                numDays -= 1

        return newDt

###############################################################################

    def is_business_day(self,
                        dt: Date):
        """ Determines if a date is a business day according to the specified
        calendar. If it is it returns True, otherwise False. """

        # For all calendars so far, SAT and SUN are not business days
        # If this ever changes I will need to add a filter here.
        if dt.is_weekend():
            return False

        if self.is_holiday(dt) is True:
            return False
        else:
            return True

###############################################################################

    def is_holiday(self,
                   dt: Date):
        """ Determines if a date is a Holiday according to the specified
        calendar. Weekends are not holidays unless the holiday falls on a 
        weekend date. """

        start_date = Date(1, 1, dt._y)
        day_in_year = dt._excel_date - start_date._excel_date + 1
        weekday = dt._weekday

        self._y = dt._y
        self._m = dt._m
        self._d = dt._d
        self._day_in_year = day_in_year
        self._weekday = weekday
        self._dt = dt

        if self._type == CalendarTypes.NONE:
            return self.holiday_none()
        elif self._type == CalendarTypes.WEEKEND:
            return self.holiday_weekend()
        elif self._type == CalendarTypes.AUSTRALIA:
            return self.holiday_australia()
        elif self._type == CalendarTypes.CANADA:
            return self.holiday_canada()
        elif self._type == CalendarTypes.FRANCE:
            return self.holiday_france()
        elif self._type == CalendarTypes.GERMANY:
            return self.holiday_germany()
        elif self._type == CalendarTypes.ITALY:
            return self.holiday_italy()
        elif self._type == CalendarTypes.JAPAN:
            return self.holiday_japan()
        elif self._type == CalendarTypes.NEW_ZEALAND:
            return self.holiday_new_zealand()
        elif self._type == CalendarTypes.NORWAY:
            return self.holiday_norway()
        elif self._type == CalendarTypes.SWEDEN:
            return self.holiday_sweden()
        elif self._type == CalendarTypes.SWITZERLAND:
            return self.holiday_switzerland()
        elif self._type == CalendarTypes.TARGET:
            return self.holiday_target()
        elif self._type == CalendarTypes.UNITED_KINGDOM:
            return self.holiday_united_kingdom()
        elif self._type == CalendarTypes.UNITED_STATES:
            return self.holiday_united_states()
        else:
            print(self._type)
            raise FinError("Unknown calendar")

###############################################################################

    def holiday_weekend(self):
        """ Weekends by themselves are a holiday. """

        if self._dt.is_weekend():
            return True
        else:
            return False

###############################################################################

    def holiday_australia(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 26:  # Australia day
            return True

        if m == 1 and d == 27 and weekday == Date.MON:  # Australia day
            return True

        if m == 1 and d == 28 and weekday == Date.MON:  # Australia day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 3:  # good friday
            return True

        if day_in_year == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # Australia day
            return True

        if m == 4 and d == 26 and weekday == Date.MON:  # Australia day
            return True

        if m == 6 and d > 7 and d < 15 and weekday == Date.MON:  # Queen
            return True

        if m == 8 and d < 8 and weekday == Date.MON:  # BANK holiday
            return True

        if m == 10 and d < 8 and weekday == Date.MON:  # BANK holiday
            return True

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

###############################################################################

    def holiday_united_kingdom(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        weekday = self._weekday
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # new years day
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em - 3:  # good friday
            return True

        if m == 5 and d <= 7 and weekday == Date.MON:
            return True

        if m == 5 and d >= 25 and weekday == Date.MON:
            return True

        if m == 6 and d == 2 and y == 2022:  # SPRING BANK HOLIDAY
            return True

        if m == 6 and d == 3 and y == 2022:  # QUEEN PLAT JUB
            return True

        if m == 8 and d > 24 and weekday == Date.MON:  # Late Summer
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday == Date.MON:  # Xmas
            return True

        if m == 12 and d == 27 and weekday == Date.TUE:  # Xmas
            return True

        if m == 12 and d == 28 and weekday == Date.MON:  # Xmas
            return True

        if m == 12 and d == 28 and weekday == Date.TUE:  # Xmas
            return True

        return False

###############################################################################

    def holiday_france(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em - 3:  # good friday
            return True

        if m == 5 and d == 1:  # LABOUR DAY
            return True

        if m == 5 and d == 8:  # VICTORY DAY
            return True

        if day_in_year == em + 39 - 1:  # Ascension
            return True

        if day_in_year == em + 50 - 1:  # pentecost
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

###############################################################################

    def holiday_sweden(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 3:  # good friday
            return True

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em + 39 - 1:  # Ascension
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

###############################################################################

    def holiday_germany(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em - 3:  # good friday
            return True

        if m == 5 and d == 1:  # LABOUR DAY
            return True

        if day_in_year == em + 39 - 1:  # Ascension
            return True

        if day_in_year == em + 50 - 1:  # pentecost
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

###############################################################################

    def holiday_switzerland(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2:  # berchtoldstag
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em - 3:  # good friday
            return True

        if day_in_year == em + 39 - 1:  # Ascension
            return True

        if day_in_year == em + 50 - 1:  # pentecost / whit
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

###############################################################################

    def holiday_japan(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # bank holiday
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # bank holiday
            return True

        if m == 1 and d > 7 and d < 15 and weekday == Date.MON:  # coa
            return True

        if m == 2 and d == 11:  # nfd
            return True

        if m == 2 and d == 12 and weekday == Date.MON:  # nfd
            return True

        if m == 2 and d == 23:  # emperor's birthday
            return True

        if m == 2 and d == 24 and weekday == Date.MON:  # emperor's birthday
            return True

        if m == 3 and d == 20:  # vernal equinox - NOT EXACT
            return True

        if m == 3 and d == 21 and weekday == Date.MON:
            return True

        if m == 4 and d == 29:  # SHOWA greenery
            return True

        if m == 4 and d == 30 and weekday == Date.MON:  # SHOWA greenery
            return True

        if m == 5 and d == 3:  # Memorial Day
            return True

        if m == 5 and d == 4:  # nation
            return True

        if m == 5 and d == 5:  # children
            return True

        if m == 5 and d == 6 and weekday == Date.MON:  # children
            return True

        if m == 7 and d > 14 and d < 22 and y != 2021 and weekday == Date.MON:
            return True

        if m == 7 and d == 22 and y == 2021:  # OLYMPICS
            return True

        if m == 7 and d == 23 and y == 2021:  # OLYMPICS HEALTH AND SPORTS HERE
            return True

        # Mountain day
        if m == 8 and d == 11 and y != 2021:
            return True

        if m == 8 and d == 12 and y != 2021 and weekday == Date.MON:
            return True

        if m == 8 and d == 9 and y == 2021 and weekday == Date.MON:
            return True

        # Respect for aged
        if m == 9 and d > 14 and d < 22 and weekday == Date.MON:
            return True

        # Equinox - APPROXIMATE
        if m == 9 and d == 23:
            return True

        if m == 9 and d == 24 and weekday == Date.MON:
            return True

        if m == 10 and d > 7 and d <= 14 and y != 2021 and weekday == Date.MON:  # HS
            return True

        if m == 11 and d == 3:  # Culture
            return True

        if m == 11 and d == 4 and weekday == Date.MON:  # Culture
            return True

        if m == 11 and d == 23:  # Thanksgiving
            return True

        return False

###############################################################################

    def holiday_new_zealand(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # new years day
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # new years day
            return True

        if m == 1 and d > 18 and d < 26 and weekday == Date.MON:  # Anniversary
            return True

        if m == 2 and d == 6:  # Waitanga day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 3:  # good friday
            return True

        if day_in_year == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # ANZAC day
            return True

        if m == 6 and d < 8 and weekday == Date.MON:  # Queen
            return True

        if m == 10 and d > 21 and d < 29 and weekday == Date.MON:  # LABOR DAY
            return True

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

###############################################################################

    def holiday_norway(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 4:  # holy thursday
            return True

        if day_in_year == em - 3:  # good friday
            return True

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em + 38:  # Ascension
            return True

        if day_in_year == em + 49:  # Pentecost
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

###############################################################################

    def holiday_united_states(self):
        """ Only bank holidays. Weekends by themselves are not a holiday.
        This is a generic US calendar that contains the superset of
        holidays for bond markets, NYSE, and public holidays. For each of
        these and other categories there will be some variations. """

        m = self._m
        d = self._d
        weekday = self._weekday

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

###############################################################################

    def holiday_canada(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        weekday = self._weekday
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == Date.MON:  # NYD
            return True

        if m == 1 and d == 3 and weekday == Date.MON:  # NYD
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == Date.MON:  # FAMILY
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 3:  # good friday
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

###############################################################################

    def holiday_italy(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em:  # Easter Monday
            return True

        if day_in_year == em - 3:  # good friday
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

###############################################################################

    def holiday_target(self):
        """ Only bank holidays. Weekends by themselves are not a holiday. """

        m = self._m
        d = self._d
        y = self._y
        day_in_year = self._day_in_year

        if m == 1 and d == 1:  # new year's day
            return True

        if m == 5 and d == 1:  # May day
            return True

        em = easterMondayDay[y - 1901]

        if day_in_year == em - 3:  # Easter Friday holiday
            return True

        if day_in_year == em:  # Easter monday holiday
            return True

        if m == 12 and d == 25:  # Xmas bank holiday
            return True

        if m == 12 and d == 26:  # Xmas bank holiday
            return True

        return False

###############################################################################

    def holiday_none(self):
        """ No day is a holiday. """
        return False

###############################################################################

    def get_holiday_list(self,
                         year: float):
        """ generates a list of holidays in a specific year for the specified
        calendar. Useful for diagnostics. """
        start_date = Date(1, 1, year)
        end_date = Date(1, 1, year + 1)
        holidayList = []
        while start_date < end_date:
            if self.is_business_day(start_date) is False and \
                    start_date.is_weekend() is False:
                holidayList.append(start_date.__str__())

            start_date = start_date.add_days(1)

        return holidayList

###############################################################################

    def easter_monday(self,
                      year: float):
        """ Get the day in a given year that is Easter Monday. This is not
        easy to compute so we rely on a pre-calculated array. """

        if year > 2100:
            raise FinError(
                "Unable to determine Easter monday in year " + str(year))

        emDays = easterMondayDay[year - 1901]
        start_date = Date(1, 1, year)
        em = start_date.add_days(emDays-1)
        return em

###############################################################################

    def __str__(self):
        s = self._type.name
        return s

###############################################################################

    def __repr__(self):
        s = self._type
        return s

###############################################################################
