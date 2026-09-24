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

###############################################################################


class FlatDiscountCurve(DiscountCurve):
    """A flat discount curve defined by a single zero rate.

    The zero rate is expressed using the specified compounding frequency.
    The curve day-count convention determines how dates are converted to
    year fractions from the curve anchor date.

    The input rate is a zero rate, not a money-market or deposit rate.
    As the zero curve is flat, no interpolation scheme is required.
    """

    ###########################################################################

    def __init__(
        self,
        anchor_dt: Date,
        flat_zero_rate: float,
        freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
        curve_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ) -> None:
        """Create a discount curve which is flat. This is very useful for
        quick testing and simply requires a curve date a rate and a compound
        frequency. As we have entered a rate, a corresponding day count
        convention must be used to specify how time periods are to be measured.
        As the curve is flat, no interpolation scheme is required.
        """

        check_argument_types(self.__init__, locals())

        self.anchor_dt = anchor_dt
        self.flat_zero_rate = flat_zero_rate
        self.freq_type = freq_type

        if not isinstance(curve_dc_type, DayCountTypes):
            raise FinError("Invalid curve day count type.")

        self.curve_dc_type = curve_dc_type

        # This is used by some inherited functions, so we choose the simplest
        self._interp_type = None

        # Set up an annual grid of times and discount factors for insight
        years = np.linspace(0.0, 5.0, 6)
        self._df_dates = self.anchor_dt.add_years(years)
        self._times = times_from_dates(self.anchor_dt,
                                       self._df_dates,
                                       self.curve_dc_type)
        self._dfs = self.df_t(self._times)

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors from scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, 0.0)

        dfs = self._zero_to_df(self.flat_zero_rate,
                               times,
                               self.freq_type)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        """Create a new FlatDiscountCurve object with the entire curve
        bumped up by the bumpsize. All other parameters are preserved."""

        disc_curve = FlatDiscountCurve(
            self.anchor_dt,
            self.flat_zero_rate + bump_size,
            freq_type=self.freq_type,
            curve_dc_type=self.curve_dc_type,
        )
        return disc_curve

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("FLAT ZERO RATE", self.flat_zero_rate)
        s += label_to_string("FREQUENCY TYPE", self.freq_type)
        s += label_to_string("CURVE DC TYPE", self.curve_dc_type)

        # Then generic DiscountCurve info
        s += "\n"
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
