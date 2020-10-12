##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

##########################################################################
# THIS IS UNFINISHED
##########################################################################

from enum import Enum

from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes

##########################################################################


class FinIborConventions():

    def __init__(self,
                 currencyName: str,
                 indexName: str = "LIBOR"):

        if currencyName == "USD" and indexName == "LIBOR":
            self._spotLag = 2
            self._dayCountType=FinDayCountTypes.THIRTY_E_360_ISDA
            self._calendarType=FinCalendarTypes.TARGET
        elif currencyName == "EUR"and indexName == "EURIBOR":
            self._spotLag = 2
            self._dayCountType=FinDayCountTypes.THIRTY_E_360_ISDA
            self._calendarType=FinCalendarTypes.TARGET
        else:
            pass

###############################################################################
