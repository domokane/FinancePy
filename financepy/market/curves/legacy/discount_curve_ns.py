# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.global_vars import G_SMALL
from ...utils.error import FinError
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...utils.helpers import label_to_string
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates

########################################################################################


class DiscountCurveNS(DiscountCurve):
    """Implementation of Nelson-Siegel parametrisation of a discount curve.
    The internal rate is a continuously compounded rate but you can calculate
    alternative frequencies by providing a corresponding compounding frequency.
    A day count convention is needed to ensure that dates are converted to the
    correct time in years. The class inherits methods from DiscountCurve."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,
        beta_0: float,
        beta_1: float,
        beta_2: float,
        tau: float,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create a Nelson-Siegel discount curve. The parameters
        beta_0, beta_1, beta_2 and tau define continuously compounded
        zero rates. Other inherited methods may convert these zero
        rates to alternative compounding conventions.
        """

        check_argument_types(self.__init__, locals())
        self._interp_type = None

        if tau <= 0:
            raise FinError("Tau must be positive")

        self.value_dt = value_dt
        self._beta_0 = beta_0
        self._beta_1 = beta_1
        self._beta_2 = beta_2
        self._tau = tau
        self.time_dc_type = time_dc_type

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        # Set up an annual grid of times and discount factors for insight
        years = np.linspace(0.0, 10.0, 11)
        self._df_dates = self.value_dt.add_years(years)
        self._times = times_from_dates( self.value_dt, self._df_dates,self.time_dc_type)
        self._dfs = self.df_t(self._times)

    ###########################################################################

    def ns_cc_zero_rate(self, times: Union[float, np.ndarray]):
        """Zero rate for Nelson-Siegel curve parametrisation. This means that
        the t vector must use the curve day count."""

        t = np.maximum(times, G_SMALL)

        theta = t / self._tau
        e = np.exp(-theta)
        zero_rate = self._beta_0
        zero_rate += self._beta_1 * (1.0 - e) / theta
        zero_rate += self._beta_2 * ((1.0 - e) / theta - e)
        return zero_rate

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors for scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        zero_rates = self.ns_cc_zero_rate(times)
        dfs = np.exp(-zero_rates * times)

        if scalar_input:
            return float(dfs[0])
        else:
            return np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):

        discount_curve = DiscountCurveNS(
            self.value_dt,
            self._beta_0 + bump_size,
            self._beta_1,
            self._beta_2,
            self._tau,
            time_dc_type=self.time_dc_type,
        )

        return discount_curve

    ####################################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("PARAMETER", "VALUE")
        s += label_to_string("BETA_0", self._beta_0)
        s += label_to_string("BETA_1", self._beta_1)
        s += label_to_string("BETA_2", self._beta_2)
        s += label_to_string("TAU", self._tau)

        # Then generic DiscountCurve info
        s += "\n"
        s += super().__repr__()
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
