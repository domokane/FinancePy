##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinInflationSwap(object):
    ''' Class for managing LIBOR forward rate agreements. A forward rate
    agreement is an agreement to exchange a fixed pre-agreed rate for a
    floating rate linked to LIBOR that is not known until some specified
    future fixing date. The FRA payment occurs on or soon after this date
    on the FRA settlement date. Typically the timing gap is two days.

    A FRA is used to hedge a Ibor quality loan or lend of some agreed
    notional amount. This period starts on the settlement date of the
    FRA and ends on the maturity date of the FRA. For example a 1x4 FRA
    relates to a Ibor starting in 1 month for a loan period ending in 4
    months. Hence it links to 3-month Ibor rate. The amount received by a 
    payer of fixed rate at settlement is:

        acc(1,2) * (Ibor(1,2) - FRA RATE) / (1 + acc(0,1) x Ibor(0,1))

    So the value at time 0 is

        acc(1,2) * (FWD Ibor(1,2) - FRA RATE) x df(0,2)

    If the base date of the curve is before the value date then we
    forward adjust this amount to that value date. For simplicity I have
    assumed that the fixing date and the settlement date are the same date.
    This should be amended later. '''

    def __init__(self,
                 startDate: FinDate,  # The date the FRA starts to accrue
                 maturityDateOrTenor: (FinDate, str),  # End of the Ibor rate period
                 fraRate: float,  # The fixed contractual FRA rate
                 dayCountType: FinDayCountTypes,  # For interest period
                 notional: float = 100.0,
                 payFixedRate: bool = True,  # True if the FRA rate is being paid
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        ''' Create a Forward Rate Agreeement object. '''

        print("DO NOT USE")
        raise FinError("DO NOT USE")
        
        checkArgumentTypes(self.__init__, locals())

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = startDate.addTenor(maturityDateOrTenor)
            calendar = FinCalendar(self._calendarType)
            maturityDate = calendar.adjust(maturityDate,
                                           self._busDayAdjustType)

        if startDate > maturityDate:
            raise FinError("Settlement date after maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._fraRate = fraRate
        self._payFixedRate = payFixedRate
        self._dayCountType = dayCountType
        self._notional = notional

    ###########################################################################

    def value(self, valuationDate, liborCurve):
        ''' Determine mark to market value of a FRA contract based on the
        market FRA rate. The same curve is used for calculating the forward
        Ibor and for doing discounting on the expected forward payment. '''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        df1 = liborCurve.df(self._startDate)
        df2 = liborCurve.df(self._maturityDate)
        liborFwd = (df1 / df2 - 1.0) / accFactor
        v = accFactor * (liborFwd - self._fraRate) * df2

#        print(df1, df2, accFactor, liborFwd, v)
        # Forward value the FRA to the value date
        df_to_valueDate = liborCurve.df(valuationDate)
        v = v * self._notional / df_to_valueDate

        if self._payFixedRate:
            v *= -1.0
        return v

    ##########################################################################

    def maturityDf(self, liborCurve):
        ''' Determine the maturity date discount factor needed to refit
        the FRA given the libor curve anbd the contract FRA rate. '''

        dc = FinDayCount(self._dayCountType)
        df1 = liborCurve.df(self._startDate)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        df2 = df1 / (1.0 + accFactor * self._fraRate)
        return df2

    ###########################################################################

    def printFlows(self, valuationDate):
        ''' Determine the value of the Deposit given a Ibor curve. '''

        flow_settle = self._notional
        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        flow_maturity = (1.0 + accFactor * self._fraRate) * self._notional

        if self._payFixedRate is True:
            print(self._startDate, -flow_settle)
            print(self._maturityDate, flow_maturity)
        else:
            print(self._startDate, flow_settle)
            print(self._maturityDate, -flow_maturity)

    ##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START ACCD DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("FRA RATE", self._fraRate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("PAY FIXED RATE", self._payFixedRate)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        s += labelToString("BUS DAY ADJUST TYPE", self._busDayAdjustType)
        s += labelToString("CALENDAR", self._calendarType)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
