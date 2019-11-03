# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import ONE_MILLION

###############################################################################


class FinLiborFRA(object):
    ''' Class for managing LIBOR forward rate agreements. A forward rate
    agreement is an agreement to exchange a fixed pre-agreed rate for a
    floating rate linked to LIBOR that is not known until some specified
    future fixing date. The FRA payment occurs on or soon after this date
    on the FRA settlement date. Typically the timing gap is two days.

    A FRA is used to hedge a Libor quality loan or lend of some agreed
    notional amount. This period starts on the settlement date of the
    FRA and ends on the maturity date of the FRA. For example a 1x4 FRA
    relates to a Libor starting in 1 month for a loan period ending in 4
    months. Hence it linkes to 3-month Libor rate.

    The amount received by a payer of fixed rate at settlement is

        acc(1,2) * (Libor(1,2) - FRA RATE) / (1 + acc(0,1) x Libor(0,1))

    So the value at time 0 is

        acc(1,2) * (FWD Libor(1,2) - FRA RATE) x df(0,2)

    If the base date of the curve is before the value date then we
    forward adjust this amount to that value date.

    For simplicity I have assumed that the fixing date and the settlement
    date are the same date. This should be amended later. '''

    def __init__(self,
                 settlementDate,  # The date on which the floating rate fixes
                 maturityDate,  # The end of the Libor rate period
                 fraRate,  # The fixed contractual FRA rate
                 payFixedRate,  # True if the FRA rate is being paid
                 dayCountType,  # For interest period between the fixing and maturity dates
                 notional=ONE_MILLION,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayConventionTypes.MODIFIED_FOLLOWING):
        ''' Create FRA object. '''

        if settlementDate > maturityDate:
            raise ValueError("Settlement date after maturity date")

        if type(payFixedRate) != bool:
            raise ValueError("PayFixedRate must be true or false.")

        if dayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Day Count Rule type " +
                str(dayCountType))

        self._calendarType = calendarType

        calendar = FinCalendar(self._calendarType)
        self._settlementDate = calendar.adjust(
            settlementDate, busDayAdjustType)
        self._maturityDate = calendar.adjust(maturityDate, busDayAdjustType)

        self._fraRate = fraRate
        self._payFixedRate = payFixedRate
        self._dayCountType = dayCountType
        self._notional = notional

    ###########################################################################

    def value(self, valueDate, liborCurve):
        ''' Determine mark to market value of a FRA contract based on the
        market FRA rate. The same curve is used for calculating the forward
        Libor and for doing discounting on the expected forward payment. '''

        dc = FinDayCount(self._dayCountType)
        accFactor0 = dc.yearFrac(self._settlementDate, self._maturityDate)
        accFactor1 = dc.yearFrac(valueDate, self._maturityDate)
        df2 = liborCurve.df(self._maturityDate)
        df1 = liborCurve.df(self._settlementDate)
        liborFwd = (df1 / df2 - 1.0) / accFactor0
        v = accFactor0 * (liborFwd - self._fraRate) * df2

        # Forward value the FRA to the value date
        df_to_valueDate = liborCurve.df(valueDate)
        v = v * self._notional / df_to_valueDate

        if self._payFixedRate:
            v *= -1.0
        return v

    ##########################################################################

    def maturityDf(self, liborCurve):
        ''' Determine the maturity date discount factor needed to refit
        the FRA given the libor curve anbd the contract FRA rate. '''

        df1 = liborCurve.df(self._settlementDate)
        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        df2 = df1 / (1.0 + accFactor * self._fraRate)
        return df2

    ##########################################################################

    def print(self):
        print("SETTLEMENT DATE:", self._settlementDate)
        print("MATURITY DATE  :", self._maturityDate)
        print("FRA RATE       :", self._fraRate)
        print("PAY FIXED LEG  :", self._payFixedRate)
        print("DAY COUNT TYPE :", self._dayCountType)

########################################################################
