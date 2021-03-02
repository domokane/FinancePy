##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

##########################################################################
# THIS IS UNFINISHED
##########################################################################

from enum import Enum

from ...utils.DayCount import DayCount, FinDayCountTypes
from ...utils.Frequency import FinFrequencyTypes
from ...utils.Calendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...utils.Calendar import FinBusDayAdjustTypes

##########################################################################


class FinIborConventions():

    def __init__(self,
                 currencyName: str,
                 indexName: str = "LIBOR"):

        if currencyName == "USD" and indexName == "LIBOR":
            self._spotLag = 2
            self._day_count_type=FinDayCountTypes.THIRTY_E_360_ISDA
            self._calendar_type=FinCalendarTypes.TARGET
        elif currencyName == "EUR"and indexName == "EURIBOR":
            self._spotLag = 2
            self._day_count_type=FinDayCountTypes.THIRTY_E_360_ISDA
            self._calendar_type=FinCalendarTypes.TARGET
        else:
            pass

###############################################################################
