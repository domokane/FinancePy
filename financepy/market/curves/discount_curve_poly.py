##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import gSmall
from ...utils.helpers import label_to_string
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...utils.frequency import FrequencyTypes
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates


###############################################################################


class DiscountCurvePoly(DiscountCurve):
    """ Zero Rate Curve of a specified frequency parametrised using a cubic
    polynomial. The zero rate is assumed to be continuously compounded but
    this can be amended by providing a frequency when extracting zero rates.
    We also need to specify a Day count convention for time calculations.
    The class inherits all of the methods from FinDiscountCurve. """

    def __init__(self,
                 valuation_date: Date,
                 coefficients: (list, np.ndarray),
                 freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_ACT_ISDA):
        """ Create zero rate curve parametrised using a cubic curve from
        coefficients and specifying a compounding frequency type and day count
        convention. """

        check_argument_types(self.__init__, locals())

        self._valuation_date = valuation_date
        self._coefficients = coefficients
        self._power = len(coefficients) - 1
        self._freq_type = freq_type
        self._day_count_type = day_count_type

    ###############################################################################

    def zero_rate(self,
                  dts: (list, Date),
                  freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
                  day_count_type: DayCountTypes = DayCountTypes.ACT_360):
        """ Calculation of zero rates with specified frequency according to
        polynomial parametrisation. This method overrides FinDiscountCurve.
        The parametrisation is not strictly in terms of continuously compounded
        zero rates, this function allows other compounding and day counts.
        This function returns a single or vector of zero rates given a vector
        of dates so must use Numpy functions. The default frequency is a
        continuously compounded rate and ACT ACT day counting. """

        if isinstance(freq_type, FrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if isinstance(day_count_type, DayCountTypes) is False:
            raise FinError("Invalid Day Count type.")

        # Get day count times to use with curve day count convention
        dc_times = times_from_dates(
            dts, self._valuation_date, self._day_count_type)

        # We now get the discount factors using these times
        zero_rates = self._zero_rate(dc_times)

        # Now get the discount factors using curve conventions
        dfs = self._zero_to_df(self._valuation_date,
                               zero_rates,
                               dc_times,
                               self._freq_type,
                               self._day_count_type)

        # Convert these to zero rates in the required frequency and day count
        zero_rates = self._df_to_zero(dfs, dts, freq_type, day_count_type)
        return zero_rates

    ###############################################################################

    def _zero_rate(self,
                   times: (float, np.ndarray)):
        """ Calculate the zero rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use. The compounding frequency defaults to that specified in the
        constructor of the curve object. Which may be annual to continuous. """

        t = np.maximum(times, gSmall)

        zero_rate = 0.0
        for n in range(0, len(self._coefficients)):
            zero_rate += self._coefficients[n] * np.power(t, n)

        return zero_rate

    ###############################################################################

    def df(self,
           dates: (list, Date)):
        """ Calculate the fwd rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use. The compounding frequency defaults to that specified in the
        constructor of the curve object. """

        # Get day count times to use with curve day count convention
        dc_times = times_from_dates(dates,
                                    self._valuation_date,
                                    self._day_count_type)

        # We now get the discount factors using these times
        zero_rates = self._zero_rate(dc_times)

        # Now get the discount factors using curve conventions
        dfs = self._zero_to_df(self._valuation_date,
                               zero_rates,
                               dc_times,
                               self._freq_type,
                               self._day_count_type)

        return dfs

    ###############################################################################

    def __repr__(self):
        """ Display internal parameters of curve. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("POWER", "COEFFICIENT")
        for i in range(0, len(self._coefficients)):
            s += label_to_string(str(i), self._coefficients[i])
        s += label_to_string("FREQUENCY", self._freq_type)

        return s

    ###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
