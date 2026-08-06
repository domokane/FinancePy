# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import G_SMALL
from ...utils.helpers import label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...utils.helpers import times_from_dates
from ...utils.day_count import DayCountTypes

###############################################################################


class DiscountCurvePoly(DiscountCurve):
    """Zero Rate Curve of a specified frequency parametrised using an
    arbitrary polynomial. The zero rate is assumed to be continuously
    compounded. The degree of the polynomial is determined by the number
    of coefficients supplied. We also need to specify a Day count
    convention for time calculations. The class inherits all of the
    methods from DiscountCurve."""

    ###########################################################################

    def __init__(
        self,
        value_dt: Date,
        coefficients: Union[list, np.ndarray],
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create zero rate curve parametrised using a cubic curve from
        coefficients and specifying a compounding frequency type and day count
        convention."""

        print("Warning: Deprecated. Use PolyDiscountCurve instead.")

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt
        self._coefficients = np.asarray(coefficients, dtype=float)

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type
        self._interp_type = None

        # Set up an annual grid of times and discount factors for insight
        years = np.linspace(0.0, 10.0, 11)
        self._df_dates = self.value_dt.add_years(years)
        self._times = times_from_dates(self.value_dt, self._df_dates, self.time_dc_type)
        self._dfs = self.df_t(self._times)

    ###########################################################################

    def poly_cc_zero_rate(self, times: Union[float, np.ndarray]):
        """Calculate cc zero rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use."""

        times, scalar_input = self._to_time_array(times)
        t = np.maximum(times, G_SMALL)

        zero_rate = 0.0
        for n, coeff in enumerate(self._coefficients):
            zero_rate += coeff * np.power(t, n)

        if scalar_input:
            return float(zero_rate[0])
        else:
            return zero_rate

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors for scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        zero_rates = self.poly_cc_zero_rate(times)
        dfs = np.exp(-zero_rates * times)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):

        bumped_coefficients = np.array(self._coefficients, copy=True)
        bumped_coefficients[0] += bump_size

        discount_curve = DiscountCurvePoly(
            self.value_dt,
            bumped_coefficients,
            time_dc_type=self.time_dc_type,
        )

        return discount_curve

    ###########################################################################

    def __repr__(self):
        """Display internal parameters of curve."""

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("POWER", "COEFFICIENT")
        for i, coeff in enumerate(self._coefficients):
            s += label_to_string(str(i), f"{coeff:12.8f}")

        # Then generic DiscountCurve info
        s += "\n"
        s += super().__repr__()
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
