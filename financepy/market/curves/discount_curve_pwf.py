# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import G_SMALL
from ...utils.math import test_monotonicity
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates
from ...market.curves.discount_curve import DiscountCurve

###############################################################################


class DiscountCurvePWF(DiscountCurve):
    """Curve is made up of a series of zero rates sections with each having
    a piecewise flat zero rate. The default compounding assumption is
    continuous. The class inherits methods from DiscountCurve."""

    ###########################################################################

    def __init__(
        self,
        value_dt: Date,
        zero_dts: list[Date],
        zero_rates: Union[list, np.ndarray],
        input_freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Creates a discount curve using a vector of times and zero rates
        that assumes that the zero rates are piecewise flat."""

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt
        self._interp_type = None

        if len(zero_dts) != len(zero_rates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zero_dts) == 0:
            raise FinError("Dates vector must have length > 0")

        self._zero_dts = zero_dts
        self._df_dates = zero_dts

        if input_freq_type == FrequencyTypes.CONTINUOUS:
            self._cc_zero_rates = np.array(zero_rates)
        else:
            f = annual_frequency(input_freq_type)
            self._cc_zero_rates = f * np.log(1.0 + np.array(zero_rates) / f)

        self.input_freq_type = input_freq_type

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        dc_times = times_from_dates(zero_dts, self.value_dt, self.time_dc_type)

        self._times = np.array(dc_times)
        self._dfs = self.df_t(self._times)

        if test_monotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

    ###########################################################################

    def pwf_cc_zero_rate(self, times: Union[float, np.ndarray, list]):
        """The piecewise flat cc zero rate is selected and returned."""

        scalar_input = np.isscalar(times)
        times = np.atleast_1d(times).astype(float)

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, G_SMALL)

        cc_zero_rates = []

        for t in times:
            l_index = 0
            found = 0

            num_times = len(self._times)
            for i in range(1, num_times):
                if self._times[i] > t:
                    l_index = i - 1
                    found = 1
                    break

            r0 = self._cc_zero_rates[l_index]

            if found == 1:
                cc_zero_rate = r0
            else:
                cc_zero_rate = self._cc_zero_rates[-1]

            cc_zero_rates.append(cc_zero_rate)

        return np.array(cc_zero_rates)

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors for scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        cc_zero_rates = self.pwf_cc_zero_rate(times)
        dfs = np.exp(-cc_zero_rates * times)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        return DiscountCurvePWF(
            self.value_dt,
            self._zero_dts.copy(),
            self._cc_zero_rates + bump_size,
            input_freq_type=FrequencyTypes.CONTINUOUS,
            time_dc_type=self.time_dc_type,
        )

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("DATE", "CC ZERO RATES")
        for dt, rate in zip(self._zero_dts, self._cc_zero_rates):
            s += label_to_string(str(dt), f"{rate:12.8f}")

        s += label_to_string("INPUT FREQUENCY", self.input_freq_type)
        s += "\n"
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
