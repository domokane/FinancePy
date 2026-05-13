# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union
import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import times_from_dates
from ...market.curves.interpolator import InterpTypes

###############################################################################


class DiscountCurveFlat(DiscountCurve):
    """A very simple discount curve based on a single zero rate with its
    own specified compounding method. Hence, the curve is assumed to be flat.
    It is used for quick and dirty analysis and when limited information is
    available. It inherits several methods from DiscountCurve."""

    ###########################################################################

    def __init__(
        self,
        value_dt: Date,
        flat_rate: Union[float, np.ndarray],
        freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create a discount curve which is flat. This is very useful for
        quick testing and simply requires a curve date a rate and a compound
        frequency. As we have entered a rate, a corresponding day count
        convention must be used to specify how time periods are to be measured.
        As the curve is flat, no interpolation scheme is required.
        """

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt
        self.flat_rate = flat_rate
        self.freq_type = freq_type

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        # This is used by some inherited functions, so we choose the simplest
        self._interp_type = InterpTypes.FLAT_FWD_RATES

        # Set up an annual grid of times and discount factors for insight
        years = np.linspace(0.0, 10.0, 11)
        self._df_dates = self.value_dt.add_years(years)
        self._times = times_from_dates(
            self._df_dates, self.value_dt, self.time_dc_type
        )
        self._dfs = self.df_t(self._times)

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors from scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, 0.0)

        dfs = self._zero_to_df(self.flat_rate, times, self.freq_type)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        """Create a new FinDiscountCurveFlat object with the entire curve
        bumped up by the bumpsize. All other parameters are preserved."""

        disc_curve = DiscountCurveFlat(
            self.value_dt,
            self.flat_rate + bump_size,
            freq_type=self.freq_type,
            time_dc_type=self.time_dc_type,
        )
        return disc_curve

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("FLAT RATE", self.flat_rate)
        s += label_to_string("FREQUENCY TYPE", self.freq_type)

        # Then generic DiscountCurve info
        s += "\n"
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
