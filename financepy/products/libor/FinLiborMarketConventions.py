# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 14:23:13 2019

@author: Dominic O'Kane
"""

##########################################################################
# THIS IS UNFINISHED
##########################################################################

from enum import Enum

from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinDayAdjustTypes

##########################################################################


class FinLiborMarketConventions():

    def __init__(self,
                 currencyName,
                 indexName = "LIBOR"):

        if currencyName == "USD" and indexName == "LIBOR":
            spotLag = 2
            dayCountType=FinDayCountTypes.THIRTY_E_360_ISDA
            calendarType=FinCalendarTypes.TARGET
        elif currencyName == "EUR"and indexName = "EURIBOR":
            spotLag = 2
            dayCountType=FinDayCountTypes.THIRTY_E_360_ISDA
            calendarType=FinCalendarTypes.TARGET
        endif

###############################################################################
