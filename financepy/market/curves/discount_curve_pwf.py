# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import G_SMALL
from ...utils.math import test_monotonicity
from ...utils.frequency import FrequencyTypes
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates
from ...market.curves.discount_curve import DiscountCurve

########################################################################################


class DiscountCurvePWF(DiscountCurve):
    """Curve is made up of a series of zero rates sections with each having
    a piecewise flat zero rate. The default compounding assumption is
    continuous. The class inherits methods from DiscountCurve."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,
        zero_dts: list,
        zero_rates: Union[list, np.ndarray],
        freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
        dc_type: DayCountTypes = DayCountTypes.ACT_ACT_ISDA,
    ):
        """Creates a discount curve using a vector of times and zero rates
        that assumes that the zero rates are piecewise flat."""

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt

        if len(zero_dts) != len(zero_rates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zero_dts) == 0:
            raise FinError("Dates vector must have length > 0")

        self._zero_dts = zero_dts
        self._zero_rates = np.array(zero_rates)
        self.freq_type = freq_type
        self.dc_type = dc_type

        dc_times = times_from_dates(zero_dts, self.value_dt, self.dc_type)

        self._times = np.array(dc_times)

        if test_monotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

    ####################################################################################

    def _zero_rate(self, times: Union[float, np.ndarray, list]):
        """The piecewise flat zero rate is selected and returned."""

        if isinstance(times, float):
            times = np.array([times])

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, G_SMALL)

        zero_rates = []

        for t in times:
            l_index = 0
            found = 0

            num_times = len(self._times)
            for i in range(1, num_times):
                if self._times[i] > t:
                    l_index = i - 1
                    found = 1
                    break

            r0 = self._zero_rates[l_index]

            if found == 1:
                zero_rate = r0
            else:
                zero_rate = self._zero_rates[-1]

            zero_rates.append(zero_rate)

        return np.array(zero_rates)

    ####################################################################################

    def _fwd(self, times: Union[np.ndarray, list]):
        """Calculate the continuously compounded forward rate at the forward
        time provided. This is done by perturbing the time by a small amount
        and measuring the change in the log of the discount factor divided by
        the time increment dt."""

        dt = 1e-6
        times = np.maximum(times, dt)

        df1 = self.df_t(times - dt)
        df2 = self.df_t(times + dt)
        fwd = np.log(df1 / df2) / (2.0 * dt)

        # Prevent the discontinuous spikes by capping them at the sum of short
        # plus the long end zero rate !
        cap = self._zero_rates[0] + self._zero_rates[-1]
        fwd = np.minimum(fwd, cap)

        return fwd

    ####################################################################################

    def df(self, dates: Union[Date, list]):
        """Return discount factors given a single or vector of dates. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA."""

        # Get day count times to use with curve day count convention
        dc_times = times_from_dates(dates, self.value_dt, self.dc_type)

        zero_rates = self._zero_rate(dc_times)

        df = self._zero_to_df(
            self.value_dt, zero_rates, dc_times, self.freq_type, self.dc_type
        )

        return df

    ####################################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("DATE", "ZERO RATE")
        for i in range(0, len(self._zero_dts)):
            s += label_to_string(self._zero_dts[i], self._zero_rates[i])
        s += label_to_string("FREQUENCY", (self.freq_type))
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
