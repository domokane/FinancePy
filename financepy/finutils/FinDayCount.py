# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:10:19 2018

@author: Dominic O'Kane
"""

from .FinDate import (FinDate, monthDaysLeapYear, monthDaysNotLeapYear)
from .FinMath import isLeapYear

from enum import Enum


# A useful source for these definitions can be found at 
# https://developers.opengamma.com/quantitative-research/Interest-Rate-Instruments-and-Market-Conventions.pdf
# and https://en.wikipedia.org/wiki/Day_count_convention
# http://www.fairmat.com/documentation/usermanual/topics/download/mediawiki/index.php/Day_Count_Conventions.htm

class FinDayCountTypes(Enum):
    THIRTY_E_360_ISDA=1 # AKA 30E/360 ISDA, EUROBOND (ISDA 2010)
    THIRTY_E_360_PLUS_ISDA=2 # NEEDS SOURCE
    ACT_ACT_ISDA=3 # DOES LEAP YEAR SPLIT
    ACT_365_ISDA=4 # SPLITS ACROSS LEAP YEAR
    THIRTY_360=5
    THIRTY_360_BOND=6
    THIRTY_E_360=7 # AKA 30/360 ICMA, EUROBOND (ISDA 2006)
    ACT_360=8
    ACT_365_FIXED=9
    ACT_365_LEAP=10

################################################################################
################################################################################

class FinDayCount(object):
    ''' Calculate the fractional day count between two dates according to a 
    specified convention. '''

    def __init__(self,dccType):
        ''' Create Day Count convention by passing in the Day Count Type. '''

        if dccType not in FinDayCountTypes:
            raise ValueError("Need to pass FinDayCountType")
            
        self._type = dccType

################################################################################

    def yearFrac(self, dt1, dt2):
        ''' Calculate the year fraction between dates dt1 and dt2 using the 
        specified day count convention. '''

        d1 = dt1._d
        d2 = dt2._d
        m1 = dt1._m
        m2 = dt2._m
        y1 = dt1._y
        y2 = dt2._y

        if self._type == FinDayCountTypes.THIRTY_360:
        
            dayDiff = 360.0 * (y2-y1) + 30.0 * (m2-m1) + (d2-d1)    
            accFactor = dayDiff/360.0    
            return accFactor

        elif self._type == FinDayCountTypes.THIRTY_360_BOND:
        
            d1 = min(d1,30)
                
            if d1 == 31 or d1 == 30:
                d2 = min(d2,30)
                    
            dayDiff = 360.0 * (y2-y1) + 30.0 * (m2-m1) + (d2-d1)    
            accFactor = dayDiff/360.0    
            return accFactor
    
        elif self._type == FinDayCountTypes.THIRTY_E_360:
        
            d1 = min(d1,30)
            d2 = min(d2,30)
            dayDiff = 360.0 * (y2-y1) + 30.0 * (m2-m1) + (d2-d1)    
            accFactor = dayDiff/360.0    
            return accFactor

        elif self._type == FinDayCountTypes.THIRTY_E_360_ISDA: 

            if isLeapYear(y1) == True:
                if d1 == monthDaysLeapYear[m1-1]:
                    d1 = 30

            if isLeapYear(y1) == False:
                if d1 == monthDaysNotLeapYear[m1-1]:
                    d1 = 30

            if isLeapYear(y2) == True:
                if d2 == monthDaysLeapYear[m2-1] and m2 != 2:
                    d2 = 30

            if isLeapYear(y2) == False:
                if d2 == monthDaysNotLeapYear[m2-1] and m2 != 2:
                    d2 = 30
                            
            # Need to exclude termination date - to check this
                                 
            dayDiff = 360.0 * (y2-y1) + 30.0 * (m2-m1) + (d2-d1)    
            accFactor = dayDiff/360.0    
            return accFactor

        elif self._type == FinDayCountTypes.THIRTY_E_360_PLUS_ISDA:

            # CHECK THIS CODE
            d1 = min(d1,30)
                
            if d2 == 31:
                d2 = 1
                m2 = m2 + 1
                                 
            dayDiff = 360.0 * (y2-y1) + 30.0 * (m2-m1) + (d2-d1)    
            accFactor = dayDiff/360.0    
            return accFactor
        
        elif self._type == FinDayCountTypes.ACT_ACT_ISDA:    

            if isLeapYear(y1):
                denom1 = 366.0
            else:
                denom1 = 365.0
    
            if isLeapYear(y2):
                denom2 = 366.0
            else:
                denom2 = 365.0

            if y1 == y2:
                accFactor = (dt2 - dt1) / denom1
                return accFactor
            else:
                daysYear1 = FinDate.datediff(dt1, FinDate(y1+1,1,1))
                daysYear2 = FinDate.datediff(FinDate(y1+1,1,1),dt2)    
                accFactor = daysYear1 / denom1
                accFactor += daysYear2 / denom2    
                return accFactor

        elif self._type == FinDayCountTypes.ACT_360:    

            accFactor = (dt2 - dt1) / 360.0    
            return accFactor    

        elif self._type == FinDayCountTypes.ACT_365_FIXED:

            accFactor = (dt2 - dt1) / 365.0    
            return accFactor    
            
        elif self._type == FinDayCountTypes.ACT_365_LEAP:    

            # The ISDA calculator sheet appears to split this across the non-leap
            # and the leap year which I do not see in any conventions.
            
            denom = 365.0
            
            if isLeapYear(y1) and dt1 <= FinDate(y1,2,28) and dt2 > FinDate(y1,2,28):
                denom = 366.0

            if isLeapYear(y2) and dt1 <= FinDate(y2,2,28) and dt2 > FinDate(y2,2,28):
                denom = 366.0
                    
            # handle case in which period straddles year end
            accFactor = (dt2 - dt1)/denom
            return accFactor

        elif self._type == FinDayCountTypes.ACT_365_ISDA:    

            if isLeapYear(y1):
                denom1 = 366.0
            else:
                denom1 = 365.0
    
            if isLeapYear(y2):
                denom2 = 366.0
            else:
                denom2 = 365.0
    
            if y1 == y2:
                accFactor = (dt2 - dt1) / denom1
                return accFactor
            else:
                daysYear1 = FinDate.datediff(dt1, FinDate(y1+1,1,1))
                daysYear2 = FinDate.datediff(FinDate(y1+1,1,1),dt2)    
                accFactor = daysYear1 / denom1
                accFactor += daysYear2 / denom2    
                return accFactor
            
        else:    
            raise ValueError(str(self._type) + " is not one of FinDayCountTypes")
            
    ###########################################################################
    
    def __str__(self):
        return str(self._type)

    ###########################################################################



