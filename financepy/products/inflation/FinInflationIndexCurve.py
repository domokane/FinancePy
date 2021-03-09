##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import testMonotonicity
from ...utils.helpers import label_to_string
from ...utils.helpers import timesFromDates
from ...utils.helpers import check_argument_types
from ...utils.date import daysInMonth
from ...utils.global_vars import gDaysInYear

###############################################################################

class FinInflationIndexCurve():
    """ This is a curve calculated from a set of dates and CPI-like numbers. It
    should start at the issue date of the bond (or index). It also requires a
    lag in months. Here is a reference to the CPI curve used for TIPS.
    
    https://www.treasury.gov/about/organizational-structure/offices/Domestic-Finance/Documents/tips-presentation.pdf
    
    """

###############################################################################

    def __init__(self,
                 indexDates: list,
                 indexValues: (list, np.ndarray),
                 lagInMonths: int = 3):

        check_argument_types(self.__init__, locals())

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

    def indexValue(self, dt: Date):
        """ Calculate index value by interpolating the CPI curve """

        lagMonthsAgoDt = dt.addMonths(-self._lagInMonths)
        
        cpiFirstDate = Date(1, lagMonthsAgoDt._m, lagMonthsAgoDt._y)
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

    def indexRatio(self, dt: Date):
        """ Calculate index value by interpolating the CPI curve """

        vt = self.indexValue(dt)
        v0 = self.indexValue(self._baseDate)
        indexRatio = vt / v0
        return indexRatio

###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("BASE DATE", self._baseDate)
        s += label_to_string("INDEX LAG", self._lagInMonths)

        s += label_to_string("DATES", "ZERO RATES")
        num_points = len(self._indexValues)
        for i in range(0, num_points):
            s += label_to_string("%12s" % self._indexDates[i],
                               "%10.7f" % self._indexValues[i])

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################

