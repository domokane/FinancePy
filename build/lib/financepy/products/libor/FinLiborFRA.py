# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

import sys
sys.path.append("..")

from finutils.FinError import FinError
from finutils.FinDayCount import FinDayCount, FinDayCountTypes

###############################################################################

class FinLiborFRA(object):
    ''' Class for managing LIBOR forward rate agreements. '''

    def __init__(self, 
                 startDate, 
                 endDate, 
                 fraRate, 
                 dayCountType ):
        ''' Create FRA object. '''
        if startDate > endDate:
            raise FinError("Settlement date after maturity date")

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Day Count Rule type " + str(dayCountType))

        self._startDate = startDate
        self._endDate = endDate
        self._fraRate = fraRate
        self._dayCountType = dayCountType

    def fwdDf(self):
        ''' Calculate the FRA rate implied forward discount factor. '''
        dayCount = FinDayCount(self._dayCountType)
        accFactor = dayCount.yearFrac(self.startDate, self.endDate)
        df = 1.0 / (1.0 + self._fraRate * accFactor)
        return df

    def value(self, valueDate, rate):
        ''' Determine value of a FRA contract based on the market FRA rate. '''
        dayCount = FinDayCount(self._dayCountType)
        accFactor0 = dayCount.yearFrac(self._startDate, self._endDate)
        accFactor1 = dayCount.yearFrac(valueDate, self._endDate)
        v = (1.0 + accFactor0 * self._fraRate) / (1.0 + accFactor1 * rate)
        return v

    def dump(self):
        print("")

########################################################################
