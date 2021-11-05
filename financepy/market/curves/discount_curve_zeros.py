##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

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
from .interpolator import InterpTypes, Interpolator


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################

class DiscountCurveZeros(DiscountCurve):
    """ This is a curve calculated from a set of dates and zero rates. As we
    have rates as inputs, we need to specify the corresponding compounding
    frequency. Also to go from rates and dates to discount factors we need to
    compute the year fraction correctly and for this we require a day count
    convention. Finally, we need to interpolate the zero rate for the times
    between the zero rates given and for this we must specify an interpolation
    convention. The class inherits methods from FinDiscountCurve. """

###############################################################################

    def __init__(self,
                 valuation_date: Date,
                 zero_dates: list,
                 zero_rates: (list, np.ndarray),
                 freq_type: FrequencyTypes = FrequencyTypes.ANNUAL,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_ACT_ISDA,
                 interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES):
        """ Create the discount curve from a vector of dates and zero rates
        factors. The first date is the curve anchor. Then a vector of zero
        dates and then another same-length vector of rates. The rate is to the
        corresponding date. We must specify the compounding frequency of the
        zero rates and also a day count convention for calculating times which
        we must do to calculate discount factors. Finally we specify the
        interpolation scheme for off-grid dates."""

        check_argument_types(self.__init__, locals())

        # Validate curve
        if len(zero_dates) == 0:
            raise FinError("Dates has zero length")

        if len(zero_dates) != len(zero_rates):
            raise FinError("Dates and Rates are not the same length")

        if freq_type not in FrequencyTypes:
            raise FinError("Unknown Frequency type " + str(freq_type))

        if day_count_type not in DayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(day_count_type))

        self._valuation_date = valuation_date
        self._freq_type = freq_type
        self._day_count_type = day_count_type
        self._interp_type = interp_type

        self._zero_rates = np.array(zero_rates)
        self._zero_dates = zero_dates

        self._times = times_from_dates(
            zero_dates, valuation_date, day_count_type)

        if test_monotonicity(self._times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        dfs = self._zero_to_df(self._valuation_date,
                               self._zero_rates,
                               self._times,
                               self._freq_type,
                               self._day_count_type)

        self._dfs = np.array(dfs)
        self._interpolator = Interpolator(self._interp_type)
        self._interpolator.fit(self._times, self._dfs)

# ###############################################################################

#     def bump(self, bump_size):
#         """ Calculate the continuous forward rate at the forward date. """

#         times = self._times.copy()
#         discount_factors = self._discount_factors.copy()

#         n = len(self._times)
#         for i in range(0, n):
#             t = times[i]
#             discount_factors[i] = discount_factors[i] * np.exp(-bump_size*t)

#         discCurve = FinDiscountCurve(self._valuation_date, times,
#                                      discount_factors,
#                                      self._interp_type)

#         return discCurve

###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUATION DATE", self._valuation_date)
        s += label_to_string("FREQUENCY TYPE", (self._freq_type))
        s += label_to_string("DAY COUNT TYPE", (self._day_count_type))
        s += label_to_string("INTERP TYPE", (self._interp_type))

        s += label_to_string("DATES", "ZERO RATES")
        num_points = len(self._times)
        for i in range(0, num_points):
            s += label_to_string("%12s" % self._zero_dates[i],
                                 "%10.7f" % self._zero_rates[i])

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
