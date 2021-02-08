###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


###############################################################################
# TODO: Do some timings and tidy up logic in adjustment function
###############################################################################

import datetime
from enum import Enum
from .FinDate import FinDate
from .FinError import FinError

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


class FinBusDayAdjustTypes(Enum):
    NONE = 1
    FOLLOWING = 2
    MODIFIED_FOLLOWING = 3
    PRECEDING = 4
    MODIFIED_PRECEDING = 5


class FinCalendarTypes(Enum):
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


class FinDateGenRuleTypes(Enum):
    FORWARD = 1
    BACKWARD = 2

###############################################################################


class FinCalendar(object):
    ''' Class to manage designation of payment dates as holidays according to
    a regional or country-specific calendar convention specified by the user.
    It also supplies an adjustment method which takes in an adjustment
    convention and then applies that to any date that falls on a holiday in the
    specified calendar. '''

    def __init__(self,
                 calendarType: FinCalendarTypes):
        ''' Create a calendar based on a specified calendar type. '''

        if calendarType not in FinCalendarTypes:
            raise FinError(
                "Need to pass FinCalendarType and not " +
                str(calendarType))

        self._type = calendarType

    ###########################################################################

    def adjust(self,
               dt: FinDate,
               busDayConventionType: FinBusDayAdjustTypes):
        ''' Adjust a payment date if it falls on a holiday according to the
        specified business day convention. '''

        if type(busDayConventionType) != FinBusDayAdjustTypes:
            raise FinError("Invalid type passed. Need FinBusDayConventionType")

        if busDayConventionType == FinBusDayAdjustTypes.NONE:
            return dt

        elif busDayConventionType == FinBusDayAdjustTypes.FOLLOWING:

            # step forward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(1)

            return dt

        elif busDayConventionType == FinBusDayAdjustTypes.MODIFIED_FOLLOWING:

            d_start = dt._d
            m_start = dt._m
            y_start = dt._y

            # step forward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(1)

            # if the business day is in a different month look back
            # for previous first business day one day at a time
            # TODO: I could speed this up by starting it at initial date
            if dt._m != m_start:
                dt = FinDate(d_start, m_start, y_start)
                while self.isBusinessDay(dt) is False:
                    dt = dt.addDays(-1)

            return dt

        elif busDayConventionType == FinBusDayAdjustTypes.PRECEDING:

            # if the business day is in the next month look back
            # for previous first business day one day at a time
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(-1)

            return dt

        elif busDayConventionType == FinBusDayAdjustTypes.MODIFIED_PRECEDING:

            d_start = dt._d
            m_start = dt._m
            y_start = dt._y

            # step backward until we find a business day
            while self.isBusinessDay(dt) is False:
                dt = dt.addDays(-1)

            # if the business day is in a different month look forward
            # for previous first business day one day at a time
            # I could speed this up by starting it at initial date
            if dt._m != m_start:
                dt = FinDate(d_start, m_start, y_start)
                while self.isBusinessDay(dt) is False:
                    dt = dt.addDays(+1)

            return dt

        else:

            raise FinError("Unknown adjustment convention" +
                           str(busDayConventionType))

        return dt

###############################################################################

    def addBusinessDays(self,
                        startDate: FinDate,
                        numDays: int):
        ''' Returns a new date that is numDays business days after FinDate. 
        All holidays in the chosen calendar are assumed not business days. '''

        # TODO: REMOVE DATETIME DEPENDENCE HERE ???

        if isinstance(numDays, int) is False:
            raise FinError("Num days must be an integer")

        dt = datetime.date(startDate._y, startDate._m, startDate._d)
        d = dt.day
        m = dt.month
        y = dt.year
        newDt = FinDate(d, m, y)

        s = +1
        if numDays < 0:
            numDays = -1 * numDays
            s = -1

        while numDays > 0:
            dt = dt + s * datetime.timedelta(days=1)
            d = dt.day
            m = dt.month
            y = dt.year
            newDt = FinDate(d, m, y)

            if self.isBusinessDay(newDt) is True:
                numDays -= 1

        return newDt

###############################################################################

    def isBusinessDay(self,
                      dt: FinDate):
        ''' Determines if a date is a business day according to the specified
        calendar. If it is it returns True, otherwise False. '''

        # For all calendars so far, SAT and SUN are not business days
        # If this ever changes I will need to add a filter here.
        if dt.isWeekend():
            return False

        if self.isHoliday(dt) is True:
            return False
        else:
            return True

###############################################################################

    def isHoliday(self,
                      dt: FinDate):
        ''' Determines if a date is a Holiday according to the specified
        calendar. Weekends are not holidays unless the holiday falls on a 
        weekend date. '''

        startDate = FinDate(1, 1, dt._y)
        dayInYear = dt._excelDate - startDate._excelDate + 1
        weekday = dt._weekday

        self._y = dt._y
        self._m = dt._m
        self._d = dt._d
        self._dayInYear = dayInYear
        self._weekday = weekday
        self._dt = dt

        if self._type == FinCalendarTypes.NONE:
            return self.HOLIDAY_NONE()
        elif self._type == FinCalendarTypes.WEEKEND:
            return self.HOLIDAY_WEEKEND()
        elif self._type == FinCalendarTypes.AUSTRALIA:
            return self.HOLIDAY_AUSTRALIA()
        elif self._type == FinCalendarTypes.CANADA:
            return self.HOLIDAY_CANADA()
        elif self._type == FinCalendarTypes.FRANCE:
            return self.HOLIDAY_FRANCE()
        elif self._type == FinCalendarTypes.GERMANY:
            return self.HOLIDAY_GERMANY()
        elif self._type == FinCalendarTypes.ITALY:
            return self.HOLIDAY_ITALY()
        elif self._type == FinCalendarTypes.JAPAN:
            return self.HOLIDAY_JAPAN()
        elif self._type == FinCalendarTypes.NEW_ZEALAND:
            return self.HOLIDAY_NEW_ZEALAND()
        elif self._type == FinCalendarTypes.NORWAY:
            return self.HOLIDAY_NORWAY()
        elif self._type == FinCalendarTypes.SWEDEN:
            return self.HOLIDAY_SWEDEN()
        elif self._type == FinCalendarTypes.SWITZERLAND:
            return self.HOLIDAY_SWITZERLAND()
        elif self._type == FinCalendarTypes.TARGET:
            return self.HOLIDAY_TARGET()
        elif self._type == FinCalendarTypes.UNITED_KINGDOM:
            return self.HOLIDAY_UNITED_KINGDOM()
        elif self._type == FinCalendarTypes.UNITED_STATES:
            return self.HOLIDAY_UNITED_STATES()
        else:
            print(self._type)
            raise FinError("Unknown calendar")

###############################################################################

    def HOLIDAY_WEEKEND(self):
        ''' Weekends by themselves are a holiday. '''

        if self._dt.isWeekend():
            return True
        else:
            return False

###############################################################################

    def HOLIDAY_AUSTRALIA(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 26:  # Australia day
            return True

        if m == 1 and d == 27 and weekday == FinDate.MON:  # Australia day
            return True

        if m == 1 and d == 28 and weekday == FinDate.MON:  # Australia day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 3:  # good friday
            return True

        if dayInYear == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # Australia day
            return True

        if m == 4 and d == 26 and weekday == FinDate.MON:  # Australia day
            return True

        if m == 6 and d > 7 and d < 15 and weekday == FinDate.MON:  # Queen 
            return True

        if m == 8 and d < 8 and weekday == FinDate.MON:  # BANK holiday 
            return True

        if m == 10 and d < 8 and weekday == FinDate.MON:  # BANK holiday 
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Boxing
            return True

        if m == 12 and d == 28 and weekday == FinDate.MON:  # Boxing
            return True

        return False
   
###############################################################################

    def HOLIDAY_UNITED_KINGDOM(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        weekday = self._weekday ; dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == FinDate.MON:  # new years day
            return True

        if m == 1 and d == 3 and weekday == FinDate.MON:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if m == 5 and d <= 7 and weekday == FinDate.MON:
            return True

        if m == 5 and d >= 25 and weekday == FinDate.MON:
            return True

        if m == 6 and d == 2 and y == 2022: # SPRING BANK HOLIDAY
            return True

        if m == 6 and d == 3 and y == 2022: # QUEEN PLAT JUB
            return True

        if m == 8 and d > 24 and weekday == FinDate.MON:  # Late Summer
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 27 and weekday == FinDate.TUE:  # Xmas
            return True

        if m == 12 and d == 28 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 28 and weekday == FinDate.TUE:  # Xmas
            return True

        return False

###############################################################################

    def HOLIDAY_FRANCE(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if m == 5 and d == 1: # LABOUR DAY
            return True

        if m == 5 and d == 8: # VICTORY DAY
            return True

        if dayInYear == em + 39 - 1:  # Ascension
            return True

        if dayInYear == em + 50 - 1:  # pentecost
            return True

        if m == 7 and d == 14: # BASTILLE DAY
            return True

        if m == 8 and d  == 15:  # ASSUMPTION
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

    def HOLIDAY_SWEDEN(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 3:  # good friday
            return True

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em + 39 - 1:  # Ascension
            return True

        if m == 5 and d == 1:  # labour day
            return True

        if m == 6 and d == 6: # June
            return True

        if m == 6 and d > 18 and d < 26 and weekday == FinDate.FRI: # Midsummer
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

    def HOLIDAY_GERMANY(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if m == 5 and d == 1: # LABOUR DAY
            return True

        if dayInYear == em + 39 - 1:  # Ascension
            return True

        if dayInYear == em + 50 - 1:  # pentecost
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

    def HOLIDAY_SWITZERLAND(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2:  # berchtoldstag
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if dayInYear == em + 39 - 1:  # Ascension
            return True

        if dayInYear == em + 50 - 1:  # pentecost / whit
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

    def HOLIDAY_JAPAN(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y;
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == FinDate.MON:  # bank holiday
            return True

        if m == 1 and d == 3 and weekday == FinDate.MON:  # bank holiday
            return True

        if m == 1 and d > 7 and d < 15 and weekday == FinDate.MON:  # coa
            return True

        if m == 2 and d == 11:  # nfd
            return True

        if m == 2 and d == 12 and weekday == FinDate.MON:  # nfd
            return True

        if m == 2 and d == 23:  # emperor's birthday
            return True

        if m == 2 and d == 24 and weekday == FinDate.MON:  # emperor's birthday
            return True

        if m == 3 and d == 20:  # vernal equinox - NOT EXACT
            return True

        if m == 3 and d == 21 and weekday == FinDate.MON:  
            return True

        if m == 4 and d == 29:  # SHOWA greenery
            return True

        if m == 4 and d == 30 and weekday == FinDate.MON:  # SHOWA greenery
            return True

        if m == 5 and d == 3:  # Memorial Day
            return True

        if m == 5 and d == 4:  # nation
            return True

        if m == 5 and d == 5:  # children
            return True

        if m == 5 and d == 6 and weekday == FinDate.MON:  # children
            return True

        if m == 7 and d > 14 and d < 22 and y != 2021 and weekday == FinDate.MON:
            return True

        if m == 7 and d == 22 and y == 2021: # OLYMPICS
            return True

        if m == 7 and d == 23 and y == 2021: # OLYMPICS HEALTH AND SPORTS HERE
            return True

        # Mountain day
        if m == 8 and d == 11 and y != 2021:
            return True

        if m == 8 and d == 12 and y != 2021 and weekday == FinDate.MON:
            return True

        if m == 8 and d == 9 and y == 2021 and weekday == FinDate.MON:
            return True

        # Respect for aged
        if m == 9 and d > 14 and d < 22 and weekday == FinDate.MON:
            return True

        # Equinox - APPROXIMATE
        if m == 9 and d == 23:
            return True

        if m == 9 and d == 24 and weekday == FinDate.MON:
            return True

        if m == 10 and d > 7 and d <= 14 and y != 2021 and weekday == FinDate.MON:  # HS
            return True

        if m == 11 and d == 3:  # Culture
            return True

        if m == 11 and d == 4 and weekday == FinDate.MON:  # Culture
            return True

        if m == 11 and d == 23:  # Thanksgiving
            return True

        return False

###############################################################################

    def HOLIDAY_NEW_ZEALAND(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear
        weekday = self._weekday

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 2 and weekday == FinDate.MON:  # new years day
            return True

        if m == 1 and d == 3 and weekday == FinDate.MON:  # new years day
            return True

        if m == 1 and d > 18 and d < 26 and weekday == FinDate.MON:  # Anniversary 
            return True

        if m == 2 and d == 6:  # Waitanga day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 3:  # good friday
            return True

        if dayInYear == em:  # Easter Monday
            return True

        if m == 4 and d == 25:  # ANZAC day
            return True

        if m == 6 and d < 8 and weekday == FinDate.MON:  # Queen 
            return True

        if m == 10 and d > 21 and d < 29 and weekday == FinDate.MON:  # LABOR DAY 
            return True

        if m == 12 and d == 25:  # Xmas
            return True

        if m == 12 and d == 26 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Xmas
            return True

        if m == 12 and d == 26:  # Boxing day
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Boxing
            return True

        if m == 12 and d == 28 and weekday == FinDate.MON:  # Boxing
            return True

        return False

###############################################################################

    def HOLIDAY_NORWAY(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 4:  # holy thursday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em + 38:  # Ascension
            return True

        if dayInYear == em + 49:  # Pentecost
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

    def HOLIDAY_UNITED_STATES(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday.
        This is a generic US calendar that contains the superset of
        holidays for bond markets, NYSE, and public holidays. For each of
        these and other categories there will be some variations. '''

        m = self._m; d = self._d; 
        weekday = self._weekday

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == FinDate.MON:  # NYD
            return True

        if m == 1 and d == 3 and weekday == FinDate.MON:  # NYD
            return True

        if m == 1 and d >= 15 and d < 22 and weekday == FinDate.MON:  # MLK
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == FinDate.MON:  # GW
            return True

        if m == 5 and d >= 25 and d <= 31 and weekday == FinDate.MON:  # MD
            return True

        if m == 7 and d == 4:  # Indep day
            return True

        if m == 7 and d == 5 and weekday == FinDate.MON:  # Indep day
            return True

        if m == 7 and d == 3 and weekday == FinDate.FRI:  # Indep day
            return True

        if m == 9 and d >= 1 and d < 8 and weekday == FinDate.MON:  # Lab
            return True

        if m == 10 and d >= 8 and d < 15 and weekday == FinDate.MON:  # CD
            return True

        if m == 11 and d == 11:  # Veterans day
            return True

        if m == 11 and d == 12 and weekday == FinDate.MON:  # Vets
            return True

        if m == 11 and d == 10 and weekday == FinDate.FRI:  # Vets
            return True

        if m == 11 and d >= 22 and d < 29 and weekday == FinDate.THU:  # TG
            return True

        if m == 12 and d == 24 and weekday == FinDate.FRI:  # Xmas holiday
            return True

        if m == 12 and d == 25:  # Xmas holiday
            return True

        if m == 12 and d == 26 and weekday == FinDate.MON:  # Xmas holiday
            return True

        if m == 12 and d == 31 and weekday == FinDate.FRI:
            return True

        return False

###############################################################################

    def HOLIDAY_CANADA(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        weekday = self._weekday; dayInYear = self._dayInYear

        if m == 1 and d == 1:  # NYD
            return True

        if m == 1 and d == 2 and weekday == FinDate.MON:  # NYD
            return True

        if m == 1 and d == 3 and weekday == FinDate.MON:  # NYD
            return True

        if m == 2 and d >= 15 and d < 22 and weekday == FinDate.MON:  # FAMILY
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 3:  # good friday
            return True

        if m == 5 and d >= 18 and d < 25 and weekday == FinDate.MON:  # VICTORIA
            return True

        if m == 7 and d == 1:  # Canada day
            return True

        if m == 7 and d == 2 and weekday == FinDate.MON:  # Canada day
            return True

        if m == 7 and d == 3 and weekday == FinDate.MON:  # Canada day
            return True

        if m == 8 and d < 8 and weekday == FinDate.MON:  # Provincial
            return True

        if m == 9 and d < 8 and weekday == FinDate.MON:  # Labor
            return True

        if m == 10 and d >= 8 and d < 15 and weekday == FinDate.MON:  # THANKS
            return True

        if m == 11 and d == 11:  # Veterans day
            return True

        if m == 11 and d == 12 and weekday == FinDate.MON:  # Vets
            return True

        if m == 11 and d == 13 and weekday == FinDate.MON:  # Vets
            return True

        if m == 12 and d == 25:  # Xmas holiday
            return True

        if m == 12 and d == 26 and weekday == FinDate.MON:  # Xmas holiday
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Xmas holiday
            return True

        if m == 12 and d == 26:  # Boxing holiday
            return True

        if m == 12 and d == 27 and weekday == FinDate.MON:  # Boxing holiday
            return True

        if m == 12 and d == 28 and weekday == FinDate.TUE:  # Boxing holiday
            return True

        return False

###############################################################################

    def HOLIDAY_ITALY(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new years day
            return True

        if m == 1 and d == 6:  # epiphany
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em:  # Easter Monday
            return True

        if dayInYear == em - 3:  # good friday
            return True

        if m == 4 and d == 25: # LIBERATION DAY
            return True

        if m == 5 and d == 1: # LABOUR DAY
            return True

        if m == 6 and d == 2 and y > 1999: # REPUBLIC DAY
            return True

        if m == 8 and d  == 15:  # ASSUMPTION
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

    def HOLIDAY_TARGET(self):
        ''' Only bank holidays. Weekends by themselves are not a holiday. '''

        m = self._m; d = self._d; y = self._y
        dayInYear = self._dayInYear

        if m == 1 and d == 1:  # new year's day
            return True

        if m == 5 and d == 1:  # May day
            return True

        em = easterMondayDay[y - 1901]

        if dayInYear == em - 3:  # Easter Friday holiday
            return True

        if dayInYear == em:  # Easter monday holiday
            return True

        if m == 12 and d == 25:  # Xmas bank holiday
            return True

        if m == 12 and d == 26:  # Xmas bank holiday
            return True

        return False

###############################################################################

    def HOLIDAY_NONE(self):
        ''' No day is a holiday. '''
        return False

###############################################################################

    def getHolidayList(self,
                       year: float):
        ''' generates a list of holidays in a specific year for the specified
        calendar. Useful for diagnostics. '''
        startDate = FinDate(1, 1, year)
        endDate = FinDate(1, 1, year+1)
        holidayList = []
        while startDate < endDate:
            if self.isBusinessDay(startDate) is False and \
              startDate.isWeekend() is False:
                holidayList.append(startDate.__str__())

            startDate = startDate.addDays(1)

        return holidayList

###############################################################################

    def easterMonday(self,
                     year: float):
        ''' Get the day in a given year that is Easter Monday. This is not
        easy to compute so we rely on a pre-calculated array. '''

        if year > 2100:
            raise FinError(
                "Unable to determine Easter monday in year " + str(year))

        emDays = easterMondayDay[year - 1901]
        startDate = FinDate(1, 1, year)
        em = startDate.addDays(emDays-1)
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
