##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.Calendar import Calendar
from ...utils.Calendar import FinCalendarTypes
from ...utils.Calendar import FinBusDayAdjustTypes
from ...utils.DayCount import DayCount, FinDayCountTypes
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinInflationSwap(object):
    """ Class for managing LIBOR forward rate agreements. A forward rate
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
    This should be amended later. """

    def __init__(self,
                 start_date: Date,  # The date the FRA starts to accrue
                 maturity_date_or_tenor: (Date, str),  # End of the Ibor rate period
                 fraRate: float,  # The fixed contractual FRA rate
                 day_count_type: FinDayCountTypes,  # For interest period
                 notional: float = 100.0,
                 payFixedRate: bool = True,  # True if the FRA rate is being paid
                 calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 bus_day_adjust_type: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        """ Create a Forward Rate Agreeement object. """

        print("DO NOT USE")
        raise FinError("DO NOT USE")
        
        checkArgumentTypes(self.__init__, locals())

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            maturity_date = start_date.addTenor(maturity_date_or_tenor)
            calendar = Calendar(self._calendar_type)
            maturity_date = calendar.adjust(maturity_date,
                                            self._bus_day_adjust_type)

        if start_date > maturity_date:
            raise FinError("Settlement date after maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._fraRate = fraRate
        self._payFixedRate = payFixedRate
        self._day_count_type = day_count_type
        self._notional = notional

    ###########################################################################

    def value(self, valuation_date, libor_curve):
        """ Determine mark to market value of a FRA contract based on the
        market FRA rate. The same curve is used for calculating the forward
        Ibor and for doing discounting on the expected forward payment. """

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        df1 = libor_curve.df(self._start_date)
        df2 = libor_curve.df(self._maturity_date)
        liborFwd = (df1 / df2 - 1.0) / acc_factor
        v = acc_factor * (liborFwd - self._fraRate) * df2

#        print(df1, df2, acc_factor, liborFwd, v)
        # Forward value the FRA to the value date
        df_to_valuation_date = libor_curve.df(valuation_date)
        v = v * self._notional / df_to_valuation_date

        if self._payFixedRate:
            v *= -1.0
        return v

    ##########################################################################

    def maturityDf(self, libor_curve):
        """ Determine the maturity date discount factor needed to refit
        the FRA given the libor curve anbd the contract FRA rate. """

        dc = DayCount(self._day_count_type)
        df1 = libor_curve.df(self._start_date)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        df2 = df1 / (1.0 + acc_factor * self._fraRate)
        return df2

    ###########################################################################

    def printFlows(self, valuation_date):
        """ Determine the value of the Deposit given a Ibor curve. """

        flow_settle = self._notional
        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        flow_maturity = (1.0 + acc_factor * self._fraRate) * self._notional

        if self._payFixedRate is True:
            print(self._start_date, -flow_settle)
            print(self._maturity_date, flow_maturity)
        else:
            print(self._start_date, flow_settle)
            print(self._maturity_date, -flow_maturity)

    ##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START ACCD DATE", self._start_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("FRA RATE", self._fraRate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("PAY FIXED RATE", self._payFixedRate)
        s += labelToString("DAY COUNT TYPE", self._day_count_type)
        s += labelToString("BUS DAY ADJUST TYPE", self._bus_day_adjust_type)
        s += labelToString("CALENDAR", self._calendar_type)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
