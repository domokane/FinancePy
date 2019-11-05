# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinDate import FinDate
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinMath import ONE_MILLION

###############################################################################


class FinLiborDeposit(object):

    def __init__(self,
                 settlementDate,
                 maturityDateOrTenor,
                 depositRate,
                 dayCountType,
                 notional=ONE_MILLION,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayConventionTypes.MODIFIED_FOLLOWING):
        ''' Create a Libor deposit object. '''

        if type(settlementDate) != FinDate:
            raise ValueError("Settlement date must be a FinDate.")

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = settlementDate.addTenor(maturityDateOrTenor)
            calendar = FinCalendar(self._calendarType)
            maturityDate = calendar.adjust(maturityDate, busDayAdjustType)

        if settlementDate > maturityDate:
            raise ValueError("Settlement date after maturity date")

        if dayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Day Count Rule type " +
                str(dayCountType))

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise ValueError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        self._settlementDate = settlementDate
        self._calendarType = calendarType
        self._dayCountType = dayCountType
        self._depositRate = depositRate
        self._notional = notional
        self._maturityDate = maturityDate

    ###########################################################################

    def maturityDf(self):
        ''' Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date.'''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        discountFactor = 1.0 / (1.0 + accFactor * self._depositRate)
        return discountFactor

    ###########################################################################

    def value(self, valueDate, liborCurve):
        ''' Determine the value of the Deposit given a Libor curve. '''

        if valueDate > self._maturityDate:
            raise ValueError("Start date after maturity date")

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        df = liborCurve.df(self._maturityDate)
        value = (1.0 + accFactor * self._depositRate) * df * self._notional

        df_settlement = liborCurve.df(self._settlementDate)
        value = value / df_settlement

        return value

    ###########################################################################

    def print(self):
        ''' Print the contractual details of the Libor deposit. '''

        print("SETTLEMENT DATE:", self._settlementDate)
        print("MATURITY DATE:", self._maturityDate)
        print("DAY COUNT TYPE:", self._dayCountType)
        print("DEPOSIT RATE:", self._depositRate)

###############################################################################
