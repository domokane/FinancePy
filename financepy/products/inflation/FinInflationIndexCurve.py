##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.math import test_monotonicity
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates
from ...utils.helpers import check_argument_types
from ...utils.date import days_in_month
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
                 index_dates: list,
                 index_values: (list, np.ndarray),
                 lag_in_months: int = 3):

        check_argument_types(self.__init__, locals())

        # Validate curve
        if len(index_dates) == 0:
            raise FinError("Dates has zero length")

        if len(index_dates) != len(index_values):
            raise FinError("Dates and Values are not the same length")

        if lag_in_months < 0:
            raise FinError("Lag must be positive.")

        self._index_dates = np.array(index_dates)
        self._index_values = np.array(index_values)
        self._lag_in_months = lag_in_months
        self._base_date = index_dates[0]

        self._indexTimes = times_from_dates(index_dates, self._base_date)

        if test_monotonicity(self._indexTimes) is False:
            raise FinError("Times or dates are not sorted in increasing order")

###############################################################################

    def index_value(self, dt: Date):
        """ Calculate index value by interpolating the CPI curve """

        lagMonthsAgoDt = dt.add_months(-self._lag_in_months)

        cpiFirstDate = Date(1, lagMonthsAgoDt._m, lagMonthsAgoDt._y)
        cpiSecondDate = cpiFirstDate.add_months(1)

        cpiFirstTime = (cpiFirstDate - self._base_date) / gDaysInYear
        cpiSecondTime = (cpiSecondDate - self._base_date) / gDaysInYear

        cpiFirstValue = np.interp(cpiFirstTime,
                                  self._indexTimes,
                                  self._index_values)

        cpiSecondValue = np.interp(cpiSecondTime,
                                   self._indexTimes,
                                   self._index_values)

        d = dt._d
        m = dt._m
        y = dt._y
        numDays = days_in_month(m, y)
        v = cpiFirstValue + (d - 1) * (cpiSecondValue -
                                       cpiFirstValue) / numDays
        return v

###############################################################################

    def index_ratio(self, dt: Date):
        """ Calculate index value by interpolating the CPI curve """

        vt = self.index_value(dt)
        v0 = self.index_value(self._base_date)
        index_ratio = vt / v0
        return index_ratio

###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("BASE DATE", self._base_date)
        s += label_to_string("INDEX LAG", self._lag_in_months)

        s += label_to_string("DATES", "ZERO RATES")
        num_points = len(self._index_values)
        for i in range(0, num_points):
            s += label_to_string("%12s" % self._index_dates[i],
                                 "%10.7f" % self._index_values[i])

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
