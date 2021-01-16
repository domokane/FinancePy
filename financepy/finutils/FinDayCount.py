##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from .FinDate import FinDate, monthDaysLeapYear, monthDaysNotLeapYear, datediff
from .FinDate import isLeapYear
from .FinError import FinError
from .FinFrequency import FinFrequencyTypes, FinFrequency
from .FinGlobalVariables import gDaysInYear

from enum import Enum

# A useful source for these definitions can be found at
# https://developers.opengamma.com/quantitative-research/Interest-Rate-Instruments-and-Market-Conventions.pdf
# and https://en.wikipedia.org/wiki/Day_count_convention
# http://www.fairmat.com/documentation/usermanual/topics/download/mediawiki/index.php/Day_Count_Conventions.htm
# http://www.eclipsesoftware.biz/DayCountConventions.html

###############################################################################


def isLastDayOfFeb(dt: FinDate):
    # Return true if we are on the last day of February
    if dt._m == 2:
        isLeap = isLeapYear(dt._y)
        if isLeap is True and dt._d == 29:
            return True
        if isLeap is False and dt._d == 28:
            return True
    else:
        return False

###############################################################################

###############################################################################
#    THIRTY_360_BOND = 1  # 30E/360 ISDA 2006 4.16f, German, Eurobond(ISDA 2000)
#    THIRTY_E_360 = 2  # ISDA 2006 4.16(g) 30/360 ISMA, ICMA
#    THIRTY_E_360_ISDA = 3  # ISDA 2006 4.16(h)
#    THIRTY_E_PLUS_360 = 4  # ROLLS D2 TO NEXT MONTH IF D2 = 31
#    ACT_ACT_ISDA = 5  # SPLITS ACCRUAL PERIOD INTO LEAP YEAR AND NON LEAP YEAR
#    ACT_ACT_ICMA = 6  # METHOD FOR ALL US TREASURY NOTES AND BONDS
#    ACT_365F = 7  # Denominator is always Fixed at 365, even in a leap year
#    ACT_360 = 8
#    ACT_365L = 9  # the 29 Feb is counted if it is in the date range
###############################################################################
        
class FinDayCountTypes(Enum):
    THIRTY_360_BOND = 1  
    THIRTY_E_360 = 2  
    THIRTY_E_360_ISDA = 3  
    THIRTY_E_PLUS_360 = 4  
    ACT_ACT_ISDA = 5  
    ACT_ACT_ICMA = 6  
    ACT_365F = 7  
    ACT_360 = 8
    ACT_365L = 9
    SIMPLE = 10 # actual divided by gDaysInYear

###############################################################################


class FinDayCount(object):
    ''' Calculate the fractional day count between two dates according to a
    specified day count convention. '''

    def __init__(self,
                 dccType: FinDayCountTypes):
        ''' Create Day Count convention by passing in the Day Count Type. '''

        if dccType not in FinDayCountTypes:
            raise FinError("Need to pass FinDayCountType")

        self._type = dccType

###############################################################################

    def yearFrac(self,
                 dt1: FinDate,    # Start of coupon period
                 dt2: FinDate,    # Settlement (for bonds) or period end(swaps)
                 dt3: FinDate = None,   # End of coupon period for accrued
                 freqType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 isTerminationDate: bool = False):  # Is dt2 a termination date
        ''' This method performs two functions:

        1) It calculates the year fraction between dates dt1 and dt2 using the
        specified day count convention which is useful for calculating year
        fractions for Libor products whose flows are day count adjusted. In
        this case we will set dt3 to be None
        2) This function is also for calculating bond accrued where dt1 is the
        last coupon date, dt2 is the settlement date of the bond and date dt3
        must be set to the next coupon date. You will also need to provide a
        coupon frequency for some conventions.

        Note that if the date is intraday, i.e. hh,mm and ss do not equal zero
        then that is used in the calculation of the year frac. This avoids 
        discontinuities for short dated intra day products. It should not
        affect normal dates for which hh=mm=ss=0.
        
        This seems like a useful source:
        https://www.eclipsesoftware.biz/DayCountConventions.html
        Wikipedia also has a decent survey of the conventions
        https://en.wikipedia.org/wiki/Day_count_convention
        and
        http://data.cbonds.info/files/cbondscalc/Calculator.pdf
        '''

        d1 = dt1._d
        m1 = dt1._m
        y1 = dt1._y

        d2 = dt2._d
        m2 = dt2._m
        y2 = dt2._y

        num = 0
        den = 0

        if self._type == FinDayCountTypes.THIRTY_360_BOND:
            # It is in accordance with section 4.16(f) of ISDA 2006 Definitions
            # Also known as 30/360, Bond Basis, 30A/360, 30-360 US Municipal
            # This method does not consider February as a special case.

            if d1 == 31:
                d1 = 30

            if d2 == 31 and d1 == 30: 
                d2 = 30

            num = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            den = 360
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.THIRTY_E_360:
            # This is in section 4.16(g) of ISDA 2006 Definitions
            # Also known as 30/360 Eurobond, 30/360 ISMA, 30/360 ICMA,
            # 30/360 European, Special German, Eurobond basis (ISDA 2006)
            # Unlike 30/360 BOND the adjustment to dt2 does not depend on dt1

            if d1 == 31:
                d1 = 30

            if d2 == 31:
                d2 = 30

            num = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            den = 360
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.THIRTY_E_360_ISDA:
            # This is 30E/360 (ISDA 2000), 30E/360 (ISDA) section 4.16(h)
            # of ISDA 2006 Definitions, German, Eurobond basis (ISDA 2000)

            if d1 == 31:
                d1 = 30

            lastDayOfFeb1 = isLastDayOfFeb(dt1)
            if lastDayOfFeb1 is True:
                d1 = 30

            if d2 == 31:
                d2 = 30

            lastDayOfFeb2 = isLastDayOfFeb(dt2)
            if lastDayOfFeb2 is True and isTerminationDate is False:
                d2 = 30

            num = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            den = 360
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.THIRTY_E_PLUS_360:

            if d1 == 31:
                d1 = 30

            if d2 == 31:
                m2 = m2 + 1  # May roll to 13 but we are doing a difference
                d2 = 1

            num = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            den = 360
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.ACT_ACT_ISDA:

            if isLeapYear(y1):
                denom1 = 366
            else:
                denom1 = 365

            if isLeapYear(y2):
                denom2 = 366
            else:
                denom2 = 365

            if y1 == y2:
                num = dt2 - dt1
                den = denom1
                accFactor = (dt2 - dt1) / denom1
                return (accFactor, num, den)
            else:
                daysYear1 = datediff(dt1, FinDate(1, 1, y1+1))
                daysYear2 = datediff(FinDate(1, 1, y2), dt2)
                accFactor1 = daysYear1 / denom1
                accFactor2 = daysYear2 / denom2
                yearDiff = y2 - y1 - 1.0
                # Note that num/den does not equal accFactor
                # I do need to pass num back
                num = daysYear1 + daysYear2
                den = denom1 + denom2
                accFactor = accFactor1 + accFactor2 + yearDiff
                return (accFactor, num, den)

        elif self._type == FinDayCountTypes.ACT_ACT_ICMA:

            freq = FinFrequency(freqType)

            if dt3 is None or freq is None:
                raise FinError("ACT_ACT_ICMA requires three dates and a freq")

            num = dt2 - dt1
            den = freq * (dt3 - dt1)
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.ACT_365F:

            num = dt2 - dt1
            den = 365
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.ACT_360:

            num = dt2 - dt1
            den = 360
            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.ACT_365L:

            # The ISDA calculator sheet appears to split this across the
            # non-leap and the leap year which I do not see in any conventions.

            freq = FinFrequency(freqType)

            if dt3 is None:
                y3 = y2
            else:
                y3 = dt3._y

            num = dt2 - dt1
            den = 365

            if isLeapYear(y1):
                feb29 = FinDate(29, 2, y1)
            elif isLeapYear(y3):
                feb29 = FinDate(29, 2, y3)
            else:
                feb29 = FinDate(1, 1, 1900)

            if freq == 1:
                if feb29 > dt1 and feb29 <= dt3:
                    den = 366
            else:
                if isLeapYear(y3) is True:
                    den = 366

            accFactor = num / den
            return (accFactor, num, den)

        elif self._type == FinDayCountTypes.SIMPLE:
            
            num = dt2 - dt1
            den = gDaysInYear
            accFactor = num / den
            return (accFactor, num, den)

        else:

            raise FinError(str(self._type) +
                           " is not one of FinDayCountTypes")

###############################################################################

    def __repr__(self):
        ''' Returns the calendar type as a string. '''
        return str(self._type)

###############################################################################
