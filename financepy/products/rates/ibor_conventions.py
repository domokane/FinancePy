##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

##########################################################################
# THIS IS UNFINISHED
##########################################################################

from ...utils.day_count import DayCountTypes
from ...utils.calendar import CalendarTypes

##########################################################################


class IborConventions():

    def __init__(self,
                 currencyName: str,
                 indexName: str = "LIBOR"):

        if currencyName == "USD" and indexName == "LIBOR":
            self._spotLag = 2
            self._day_count_type = DayCountTypes.THIRTY_E_360_ISDA
            self._calendar_type = CalendarTypes.TARGET
        elif currencyName == "EUR" and indexName == "EURIBOR":
            self._spotLag = 2
            self._day_count_type = DayCountTypes.THIRTY_E_360_ISDA
            self._calendar_type = CalendarTypes.TARGET
        else:
            pass

###############################################################################
