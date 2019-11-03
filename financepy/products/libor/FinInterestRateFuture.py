# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

# TODO: Add functionality around settlement
# TODO: Write test function

from math import exp

from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinGlobalVariables import gDaysInYear

###############################################################################


class FinInterestRateFuture(object):
    ''' Class for managing short term interest rate futures contracts. '''

    # Reference
    # https://www.cmegroup.com/education/files/eurodollar-futures-the-basics-file01.pdf

    def __init__(self,
                 lastTradingDate,
                 dayCountType,
                 contractSize):
        ''' Create an interest rate futures contract.'''

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Day Count Rule type " + str(dayCountType))

        # contract settles 2 days after last trading date
        self._lastTradingDate = lastTradingDate
        self._dayCountType = dayCountType
        self._contractSize = contractSize

        # CHECK
        self._lastSettlementDate = lastTradingDate.addDays(2)
        self._endOfInterestRatePeriod = lastTradingDate.addMonths(3)

###############################################################################

    def futuresRate(self,
                    settlementDate,
                    futuresPrice):
        ''' Calculate implied futures rate from the futures price.'''
        futuresRate = (100.0 - futuresPrice) / 100.0
        return futuresRate

###############################################################################

    def convexity(self,
                  settlementDate,
                  volatility,
                  a):
        ''' Calculation of the convexity adjustment between FRAs and interest
        rate futures using the Hull-White model as described in technical note. '''

        # Technical note
        # http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote1.pdf

        t1 = (self._lastTradingDate - settlementDate) / gDaysInYear
        t2 = (self._endOfInterestRatePeriod - settlementDate) / gDaysInYear

        # Hull White model for short rate dr = (theta(t)-ar) dt + sigma * dz
        # This reduces to Ho-Lee when a = 0 so to avoid divergences I provide
        # this numnerical limit
        if abs(a) > 1e-10:
            bt1t2 = (1.0 - exp(-a * (t2 - t1))) / a
            b0t1 = (1.0 - exp(-a * t1)) / a
            term = bt1t2 * (1.0 - exp(-2.0 * a * t1)) + 2.0 * a * (b0t1**2)
            c = bt1t2 * term * (volatility**2) / (t2 - t1) / 4.0 / a

            term = (t2 - t1) * 2.0 * a * t1 + 2.0 * a * t1
            c = (t2 - t1) * term * (volatility**2) / (t2 - t1) / 4.0 / a
        else:
            c = t1 * t2 * (volatility**2) / 2.0

        return c

##########################################################################

    def print(self):
        print("")

##########################################################################
