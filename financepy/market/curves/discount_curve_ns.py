##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import gSmall
from ...utils.error import FinError
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import check_argument_types
from ...utils.helpers import label_to_string
from ...utils.day_count import DayCountTypes
from ...utils.helpers import times_from_dates


###############################################################################


class DiscountCurveNS(DiscountCurve):
    """ Implementation of Nelson-Siegel parametrisation of a discount curve.
    The internal rate is a continuously compounded rate but you can calculate
    alternative frequencies by providing a corresponding compounding frequency.
    A day count convention is needed to ensure that dates are converted to the
    correct time in years. The class inherits methods from FinDiscountCurve."""

    def __init__(self,
                 valuation_date: Date,
                 beta0: float,
                 beta1: float,
                 beta2: float,
                 tau: float,
                 freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_ACT_ISDA):
        """ Creation of a FinDiscountCurveNS object. Parameters are provided
        individually for beta0, beta1, beta2 and tau. The zero rates produced
        by this parametrisation have an implicit compounding convention that
        defaults to continuous but which can be overridden. """

        check_argument_types(self.__init__, locals())

        if tau <= 0:
            raise FinError("Tau must be positive")

        self._valuation_date = valuation_date
        self._beta0 = beta0
        self._beta1 = beta1
        self._beta2 = beta2
        self._tau = tau
        self._freq_type = freq_type
        self._day_count_type = day_count_type

    ###############################################################################

    def zero_rate(self,
                  dates: (list, Date),
                  freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
                  day_count_type: DayCountTypes = DayCountTypes.ACT_360):
        """ Calculation of zero rates with specified frequency according to
        NS parametrisation. This method overrides that in FinDiscountCurve.
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

        # Convert these to zero rates in the required frequency and day count
        zero_rates = self._df_to_zero(dfs,
                                      dates,
                                      freq_type,
                                      day_count_type)

        return zero_rates

    ###############################################################################

    def _zero_rate(self,
                   times: (float, np.ndarray)):
        """ Zero rate for Nelson-Siegel curve parametrisation. This means that
        the t vector must use the curve day count."""

        t = np.maximum(times, gSmall)

        theta = t / self._tau
        e = np.exp(-theta)
        zero_rate = self._beta0
        zero_rate += self._beta1 * (1.0 - e) / theta
        zero_rate += self._beta2 * ((1.0 - e) / theta - e)
        return zero_rate

    ###############################################################################

    def df(self,
           dates: (Date, list)):
        """ Return discount factors given a single or vector of dates. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA. """

        # Get day count times to use with curve day count convention
        dc_times = times_from_dates(dates,
                                    self._valuation_date,
                                    self._day_count_type)

        zero_rates = self._zero_rate(dc_times)

        df = self._zero_to_df(self._valuation_date,
                              zero_rates,
                              dc_times,
                              self._freq_type,
                              self._day_count_type)

        return df

    ###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("PARAMETER", "VALUE")
        s += label_to_string("BETA0", self._beta0)
        s += label_to_string("BETA1", self._beta1)
        s += label_to_string("BETA2", self._beta2)
        s += label_to_string("TAU", self._tau)
        s += label_to_string("FREQUENCY", (self._freq_type))
        s += label_to_string("DAY_COUNT", (self._day_count_type))
        return s

    ###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
