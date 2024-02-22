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
                 currency_name: str,
                 index_name: str = "LIBOR"):

        if currency_name == "USD" and index_name == "LIBOR":
            self._spot_lag = 2
            self._dc_type = DayCountTypes.THIRTY_E_360_ISDA
            self._cal_type = CalendarTypes.TARGET
        elif currency_name == "EUR" and index_name == "EURIBOR":
            self._spot_lag = 2
            self._dc_type = DayCountTypes.THIRTY_E_360_ISDA
            self._cal_type = CalendarTypes.TARGET
        else:
            pass

###############################################################################
