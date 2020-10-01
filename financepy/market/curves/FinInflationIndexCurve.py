##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import timesFromDates
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDate import daysInMonth
from ...finutils.FinGlobalVariables import gDaysInYear

###############################################################################

class FinInflationIndexCurve():
    ''' This is a curve calculated from a set of dates and CPI-like numbers. It
    should start at the issue date of the bond (or index). It also requires a
    lag in months. Here is a reference to the CPI curve used for TIPS.
    
    https://www.treasury.gov/about/organizational-structure/offices/Domestic-Finance/Documents/tips-presentation.pdf
    
    '''

###############################################################################

    def __init__(self,
                 indexDates: list,
                 indexValues: (list, np.ndarray),
                 lagInMonths: int = 3):

        checkArgumentTypes(self.__init__, locals())

        # Validate curve
        if len(indexDates) == 0:
            raise FinError("Dates has zero length")

        if len(indexDates) != len(indexValues):
            raise FinError("Dates and Values are not the same length")

        if lagInMonths < 0:
            raise FinError("Lag must be positive.")

        self._indexDates = np.array(indexDates)
        self._indexValues = np.array(indexValues)
        self._lagInMonths = lagInMonths
        self._baseDate = indexDates[0]

        self._indexTimes = timesFromDates(indexDates, self._baseDate)

        if testMonotonicity(self._indexTimes) is False:
            raise FinError("Times or dates are not sorted in increasing order")

###############################################################################

    def indexValue(self, dt: FinDate):
        ''' Calculate index value by interpolating the CPI curve '''

        lagMonthsAgoDt = dt.addMonths(-self._lagInMonths)
        
        cpiFirstDate = FinDate(1, lagMonthsAgoDt._m, lagMonthsAgoDt._y) 
        cpiSecondDate = cpiFirstDate.addMonths(1)
        
        cpiFirstTime = (cpiFirstDate - self._baseDate) / gDaysInYear
        cpiSecondTime = (cpiSecondDate - self._baseDate) / gDaysInYear
       
        cpiFirstValue = np.interp(cpiFirstTime,
                                  self._indexTimes,
                                  self._indexValues)

        cpiSecondValue = np.interp(cpiSecondTime,
                                   self._indexTimes,
                                   self._indexValues)
        
        d = dt._d
        m = dt._m
        y = dt._y       
        numDays = daysInMonth(m, y)
        v = cpiFirstValue + (d - 1) * (cpiSecondValue - cpiFirstValue) / numDays       
        return v

###############################################################################

    def indexRatio(self, dt: FinDate):
        ''' Calculate index value by interpolating the CPI curve '''

        vt = self.indexValue(dt)
        v0 = self.indexValue(self._baseDate)
        indexRatio = vt / v0
        return indexRatio

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("BASE DATE", self._baseDate)
        s += labelToString("INDEX LAG", self._lagInMonths)

        s += labelToString("DATES", "ZERO RATES")
        numPoints = len(self._indexValues)
        for i in range(0, numPoints):
            s += labelToString("%12s" % self._indexDates[i],
                               "%10.7f" % self._indexValues[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################

