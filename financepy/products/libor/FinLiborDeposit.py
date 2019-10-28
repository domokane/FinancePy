# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinDayCount import FinDayCountTypes

###############################################################################

class FinLiborDeposit(object):

    def __init__(self, 
                 settlementDate, 
                 maturityDate, 
                 depositRate, 
                 dayCountType, 
                 calendarType = FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayConventionTypes.MODIFIED_FOLLOWING):
        
        if settlementDate > maturityDate:
            raise ValueError("Settlement date after maturity date")

#        for dcType in FinDayCountTypes:
#           print(type(dcType),type(dayCountType),dcType == dayCountType)

        if dayCountType not in FinDayCountTypes:
            raise ValueError("Unknown Day Count Rule type " + str(dayCountType))
        
        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise ValueError("Unknown Business Day Adjust type " + str(busDayAdjustType))

        self._settlementDate = settlementDate        
        self._calendarType = calendarType  
        self._dayCountType = dayCountType
        self._depositRate = depositRate

        calendar = FinCalendar(self._calendarType)
        maturityDate = calendar.adjust(maturityDate, busDayAdjustType)
        self._maturityDate = maturityDate
        
    def df(self):
        dc = FinDayCount(self._dayCountType)
        accFactor = dc.year_frac(self._settlementDate, self._maturityDate)
        discountFactor = 1.0 / (1.0 + accFactor * self._depositRate)
        return discountFactor

    def value(self, valueDate, rate):

        if valueDate > self._maturityDate:
            raise ValueError("Start date after maturity date")

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.year_frac(self._settlementDate, self._maturityDate)
        value = (1.0 + self._depositRate) / (1.0 + accFactor * rate)
        return value

    def dump(self):

        print(self._settlementDate)
        print(self._maturityDate)
        print(self._dayCountType)
        print(self._depositRate)
                        
########################################################################
########################################################################
