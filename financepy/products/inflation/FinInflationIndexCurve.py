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
from ...utils.global_vars import g_days_in_year

###############################################################################


class FinInflationIndexCurve:
    """This is a curve calculated from a set of dates and CPI-like numbers. It
    should start at the issue date of the bond (or index). It also requires a
    lag in months. Here is a reference to the CPI curve used for TIPS.

    https://www.treasury.gov/about/organizational-structure/offices/Domestic-Finance/Documents/tips-presentation.pdf

    """

    ###############################################################################

    def __init__(
        self,
        index_dts: list,
        index_values: (list, np.ndarray),
        lag_in_months: int = 3,
    ):

        check_argument_types(self.__init__, locals())

        # Validate curve
        if len(index_dts) == 0:
            raise FinError("Dates has zero length")

        if len(index_dts) != len(index_values):
            raise FinError("Dates and Values are not the same length")

        if lag_in_months < 0:
            raise FinError("Lag must be positive.")

        self.index_dts = np.array(index_dts)
        self.index_values = np.array(index_values)
        self.lag_in_months = lag_in_months
        self.base_dt = index_dts[0]
        self.index_times = times_from_dates(index_dts, self.base_dt)

        if test_monotonicity(self.index_times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

    ###########################################################################

    def index_value(self, dt: Date):
        """Calculate index value by interpolating the CPI curve"""

        lagMonthsAgoDt = dt.add_months(-self.lag_in_months)

        cpi_first_dt = Date(1, lagMonthsAgoDt.m, lagMonthsAgoDt.y)
        cpi_second_dt = cpi_first_dt.add_months(1)

        cpi_first_time = (cpi_first_dt - self.base_dt) / g_days_in_year
        cpi_second_time = (cpi_second_dt - self.base_dt) / g_days_in_year

        cpi_first_value = np.interp(
            cpi_first_time, self.index_times, self.index_values
        )

        cpi_second_value = np.interp(
            cpi_second_time, self.index_times, self.index_values
        )

        d = dt.d
        m = dt.m
        y = dt.y
        num_days = days_in_month(m, y)
        v = (
            cpi_first_value
            + (d - 1) * (cpi_second_value - cpi_first_value) / num_days
        )
        return v

    ###########################################################################

    def index_ratio(self, dt: Date):
        """Calculate index value by interpolating the CPI curve"""

        vt = self.index_value(dt)
        v0 = self.index_value(self.base_dt)
        index_ratio = vt / v0
        return index_ratio

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("BASE DATE", self.base_dt)
        s += label_to_string("INDEX LAG", self.lag_in_months)

        s += label_to_string("DATES", "ZERO RATES")
        num_points = len(self.index_values)
        for i in range(0, num_points):
            s += label_to_string(
                "%12s" % self.index_dts[i], "%10.7f" % self.index_values[i]
            )

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
