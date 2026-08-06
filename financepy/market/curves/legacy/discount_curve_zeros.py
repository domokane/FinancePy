# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np


from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.math import test_monotonicity
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...market.curves.interpolator import InterpTypes, Interpolator

import matplotlib.pyplot as plt

# TODO: Fix up __repr__ function

###############################################################################


class DiscountCurveZeros(DiscountCurve):
    """This is a curve calculated from a set of dates and zero rates. As we
    have rates as inputs, we need to specify the corresponding compounding
    frequency. Also to go from rates and dates to discount factors we need to
    compute the year fraction correctly and for this we require a day count
    convention. Finally, we need to interpolate the zero rate for the times
    between the zero rates given and for this we must specify an interpolation
    convention. The class inherits methods from DiscountCurve."""

    ###########################################################################

    def __init__(
        self,
        value_dt: Date,
        zero_dts: list,
        zero_rates: Union[list, np.ndarray],
        freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create the discount curve from a vector of dates and zero rates
        factors. The first date is the curve anchor. Then a vector of zero
        dates and then another same-length vector of rates. The rate is to the
        corresponding date. We must specify the compounding frequency of the
        zero rates and also a day count convention for calculating times which
        we must do to calculate discount factors. Finally we specify the
        interpolation scheme for off-grid dates."""

        check_argument_types(self.__init__, locals())

        # Validate curve
        if len(zero_dts) == 0:
            raise FinError("Dates has zero length")

        if len(zero_dts) != len(zero_rates):
            raise FinError("Dates and Rates are not the same length")

        if freq_type not in FrequencyTypes:
            raise FinError("Unknown Frequency type " + str(freq_type))

        self.value_dt = value_dt
        self.freq_type = freq_type

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        if zero_dts[0] <= value_dt:
            raise FinError("Zero rate dates must be strictly after the valuation date.")

        self.time_dc_type = time_dc_type
        zero_rates = np.asarray(zero_rates, dtype=float)

        zero_times = times_from_dates(value_dt, zero_dts, time_dc_type)
        zero_times = np.asarray(zero_times, dtype=float)

        if test_monotonicity(zero_times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        dfs = self._zero_to_df(zero_rates, zero_times, self.freq_type)

        self._times = np.concatenate(([0.0], zero_times))
        self._dfs = np.concatenate(([1.0], dfs))
        self._df_dates = [value_dt] + list(zero_dts)
        self._zero_dts = [value_dt] + list(zero_dts)
        self._zero_rates = np.concatenate(([np.nan], zero_rates))

        self._interp_type = interp_type
        self._interpolator = Interpolator(self._interp_type)
        self.fit(self._times, self._dfs)

    ###########################################################################

    def bump_parallel(self, bump_size: float):
        """Return a new curve with all zero rates bumped in parallel."""

        bumped_zero_rates = self._zero_rates + bump_size

        return ZeroRatesDiscountCurve(
            self.value_dt,
            self._zero_dts.copy(),
            bumped_zero_rates,
            freq_type=self.freq_type,
            interp_type=self._interp_type,
            time_dc_type=self.time_dc_type,
        )

    ###########################################################################

    def plot(
        self,
        title,
        times: np.ndarray = None,
        ymin: float = None,
        ymax: float = None,
        filename: str = None,
    ):
        """Display discount curve."""

        plt.rcParams.update(
            {
                "lines.linewidth": 3,
                "font.size": 14,
                "axes.labelsize": 14,
                "axes.titlesize": 16,
                "legend.fontsize": 14,
            }
        )

        plt.figure(figsize=(12, 6))
        plt.title(title)

        if times is None:
            tmax = np.max(self._times)
            times = np.linspace(0.0, tmax, int(tmax * 12.0))
        else:
            times = np.asarray(times, dtype=float)

        if np.any(times < 0.0):
            raise FinError("Plot times must be strictly positive.")

        times = np.maximum(times, 1e-8)

        # Annual zero rates
        zeros = self.zero_rate_t(times, self.freq_type)
        cc_fwds = self.fwd_rate_inst_t(times)

        # Bump the forwards by 1bp in case the curve is flat so they can be seen
        cc_fwds = cc_fwds + 1e-4

        plt.plot(times, zeros * 100, label="Zero Rates", color="blue")
        plt.plot(times, cc_fwds * 100, label="Inst Fwd Rates", color="orange")
        plt.plot(
            self._times,
            self._zero_rates * 100,
            "o",
            markersize=10,
            label="Input Zero Rates",
            color="blue",
        )

        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Rate (%)")
        plt.legend()

        if ymin is not None and ymax is not None:
            plt.ylim(ymin, ymax)

        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if filename is not None:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.02)

        plt.show()

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ZERO RATE FREQUENCY", self.freq_type)
        s += label_to_string("DATES", "ZERO RATES")

        for dt, rate in zip(self._zero_dts, self._zero_rates):
            s += label_to_string(str(dt), f"{rate:12.8f}")

        s += "\n"
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
