##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.Frequency import FinFrequencyTypes
from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.DayCount import FinDayCountTypes
from ...utils.Math import testMonotonicity
from ...utils.FinHelperFunctions import labelToString
from ...utils.FinHelperFunctions import timesFromDates
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...utils.FinHelperFunctions import checkArgumentTypes
from .FinInterpolator import FinInterpTypes, FinInterpolator


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################

class FinDiscountCurveZeros(FinDiscountCurve):
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
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 freq_type: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 day_count_type: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA,
                 interpType: FinInterpTypes = FinInterpTypes.FLAT_FWD_RATES):
        """ Create the discount curve from a vector of dates and zero rates
        factors. The first date is the curve anchor. Then a vector of zero
        dates and then another same-length vector of rates. The rate is to the
        corresponding date. We must specify the compounding frequency of the
        zero rates and also a day count convention for calculating times which
        we must do to calculate discount factors. Finally we specify the
        interpolation scheme for off-grid dates."""

        checkArgumentTypes(self.__init__, locals())

        # Validate curve
        if len(zeroDates) == 0:
            raise FinError("Dates has zero length")

        if len(zeroDates) != len(zeroRates):
            raise FinError("Dates and Rates are not the same length")

        if freq_type not in FinFrequencyTypes:
            raise FinError("Unknown Frequency type " + str(freq_type))

        if day_count_type not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(day_count_type))

        self._valuation_date = valuation_date
        self._freq_type = freq_type
        self._day_count_type = day_count_type
        self._interpType = interpType

        self._zeroRates = np.array(zeroRates)
        self._zeroDates = zeroDates

        self._times = timesFromDates(zeroDates, valuation_date, day_count_type)

        if testMonotonicity(self._times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        dfs = self._zeroToDf(self._valuation_date,
                             self._zeroRates,
                             self._times,
                             self._freq_type,
                             self._day_count_type)

        self._dfs = np.array(dfs)
        self._interpolator = FinInterpolator(self._interpType)
        self._interpolator.fit(self._times, self._dfs)

# ###############################################################################

#     def bump(self, bumpSize):
#         """ Calculate the continuous forward rate at the forward date. """

#         times = self._times.copy()
#         discountFactors = self._discountFactors.copy()

#         n = len(self._times)
#         for i in range(0, n):
#             t = times[i]
#             discountFactors[i] = discountFactors[i] * np.exp(-bumpSize*t)

#         discCurve = FinDiscountCurve(self._valuation_date, times,
#                                      discountFactors,
#                                      self._interpType)

#         return discCurve

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUATION DATE", self._valuation_date)
        s += labelToString("FREQUENCY TYPE", (self._freq_type))
        s += labelToString("DAY COUNT TYPE", (self._day_count_type))
        s += labelToString("INTERP TYPE", (self._interpType))

        s += labelToString("DATES", "ZERO RATES")
        num_points = len(self._times)
        for i in range(0, num_points):
            s += labelToString("%12s" % self._zeroDates[i],
                               "%10.7f" % self._zeroRates[i])

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################

