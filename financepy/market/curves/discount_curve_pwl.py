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

########################################################################################


class DiscountCurvePWL(DiscountCurve):
    """Curve is made up of a series of sections assumed to each have a
    piece-wise linear zero rate. The zero rate has a specified frequency
    which defaults to continuous. This curve inherits all of the extra methods
    from DiscountCurve."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,
        zero_dts: Union[Date, list],
        zero_rates: Union[list, np.ndarray],
        input_freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Curve is defined by a vector of increasing times and zero rates."""

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt
        self._interp_type = None

        if len(zero_dts) != len(zero_rates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zero_dts) < 2:
            raise FinError("Dates vector must have length at least 2")

        if input_freq_type == FrequencyTypes.CONTINUOUS:
            self._cc_zero_rates = np.array(zero_rates)
        else:
            f = annual_frequency(input_freq_type)
            self._cc_zero_rates = f * np.log(1.0 + np.array(zero_rates) / f)

        self._zero_dts = zero_dts
        self._df_dates = self._zero_dts
        self.input_freq_type = input_freq_type

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        dc_times = times_from_dates(zero_dts, self.value_dt, self.time_dc_type)
        self._times = np.array(dc_times)
        if test_monotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._dfs = self.df_t(self._times)

    ###########################################################################

    def pwl_cc_zero_rate(self, times: Union[list, np.ndarray]):
        """Calculate the piecewise linear cc zero rate. This is taken from the
        initial inputs. A simple linear interpolation scheme is used. If the
        user supplies a frequency type then a conversion is done."""

        scalar_input = np.isscalar(times)

        if np.isscalar(times):
            times = np.array([times])
        else:
            times = np.array(times)

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

            t0 = self._times[l_index]
            r0 = self._cc_zero_rates[l_index]
            t1 = self._times[l_index + 1]
            r1 = self._cc_zero_rates[l_index + 1]

            if found == 1:
                cc_zero_rate = ((t1 - t) * r0 + (t - t0) * r1) / (t1 - t0)
            else:
                cc_zero_rate = self._cc_zero_rates[-1]

            cc_zero_rates.append(cc_zero_rate)

        cc_zero_rates = np.asarray(cc_zero_rates, dtype=float)

        if scalar_input:
            return np.array(cc_zero_rates[0])
        else:
            return cc_zero_rates

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors for scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        cc_zero_rates = self.pwl_cc_zero_rate(times)
        dfs = np.exp(-cc_zero_rates * times)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        return DiscountCurvePWL(
            self.value_dt,
            self._zero_dts.copy(),
            self._cc_zero_rates + bump_size,
            input_freq_type=FrequencyTypes.CONTINUOUS,
            time_dc_type=self.time_dc_type,
        )

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("DATE", "CC ZERO RATE")
        for i in range(0, len(self._zero_dts)):
            s += label_to_string(self._zero_dts[i], self._cc_zero_rates[i])
        s += label_to_string("INPUT FREQUENCY", self.input_freq_type)
        s += "\n"

        # Then generic DiscountCurve info
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
