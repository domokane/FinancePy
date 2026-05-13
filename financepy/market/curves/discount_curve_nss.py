# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.global_vars import G_SMALL
from ...utils.helpers import label_to_string
from ...utils.error import FinError
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates

########################################################################################


class DiscountCurveNSS(DiscountCurve):
    """Implementation of Nelson-Siegel-Svensson parametrisation of the
    zero rate curve. The zero rate is assumed to be continuously compounded.
    This can be changed when calling for zero rates. A day count convention is
    needed to ensure that dates are converted to the correct time in years. The
    class inherits methods from DiscountCurve."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,
        beta_0: float,
        beta_1: float,
        beta_2: float,
        beta_3: float,
        tau_1: float,
        tau_2: float,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create a FinDiscountCurveNSS object by passing in curve valuation
        date plus the 4 different beta values and the 2 tau values. The zero
        rates produced by this parametrisation have an implicit compounding
        convention that defaults to continuous but can be overriden."""

        check_argument_types(self.__init__, locals())

        if tau_1 <= 0:
            raise FinError("Tau1 must be positive")

        if tau_2 <= 0:
            raise FinError("Tau2 must be positive")

        self.value_dt = value_dt
        self._beta_0 = beta_0
        self._beta_1 = beta_1
        self._beta_2 = beta_2
        self._beta_3 = beta_3
        self._tau_1 = tau_1
        self._tau_2 = tau_2
        self._interp_type = None

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        # Set up an annual grid of times and discount factors for insight
        years = np.linspace(0.0, 10.0, 11)
        self._df_dates = self.value_dt.add_years(years)
        self._times = times_from_dates(
            self._df_dates, self.value_dt, self.time_dc_type
        )
        self._dfs = self.df_t(self._times)

    ####################################################################################

    def nss_cc_zero_rate(self, times: Union[float, np.ndarray]):
        """Calculation of zero rates given a single time or a numpy vector of
        times. This function can return a single zero rate or a vector of zero
        rates. The compounding frequency must be provided."""

        t = np.maximum(times, G_SMALL)

        theta_1 = t / self._tau_1
        theta_2 = t / self._tau_2
        e_1 = np.exp(-theta_1)
        e_2 = np.exp(-theta_2)
        zero_rate = self._beta_0
        zero_rate += self._beta_1 * (1.0 - e_1) / theta_1
        zero_rate += self._beta_2 * ((1.0 - e_1) / theta_1 - e_1)
        zero_rate += self._beta_3 * ((1.0 - e_2) / theta_2 - e_2)
        return zero_rate

    ###########################################################################

    def df_t(self, t: Union[float, list, np.ndarray]):
        """Return discount factors for scalar or vector times."""

        times, scalar_input = self._to_time_array(t)
        times = np.maximum(times, G_SMALL)

        zero_rates = self.nss_cc_zero_rate(times)
        dfs = np.exp(-zero_rates * times)

        return float(dfs[0]) if scalar_input else np.asarray(dfs, dtype=float)

    ###########################################################################

    def bump_parallel(self, bump_size: float):

        discount_curve = DiscountCurveNSS(
            self.value_dt,
            self._beta_0 + bump_size,
            self._beta_1,
            self._beta_2,
            self._beta_3,
            self._tau_1,
            self._tau_2,
            time_dc_type=self.time_dc_type,
        )

        return discount_curve

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("PARAMETER", "VALUE")
        s += label_to_string("BETA0", f"{self._beta_0:12.8f}")
        s += label_to_string("BETA1", f"{self._beta_1:12.8f}")
        s += label_to_string("BETA2", f"{self._beta_2:12.8f}")
        s += label_to_string("BETA3", f"{self._beta_3:12.8f}")
        s += label_to_string("TAU1", f"{self._tau_1:12.8f}")
        s += label_to_string("TAU2", f"{self._tau_2:12.8f}")

        # Then generic DiscountCurve info
        s += "\n"
        s += super().__repr__()
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
