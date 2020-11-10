##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import timesFromDates
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from .FinInterpolator import FinInterpTypes, FinInterpolator


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################

class FinDiscountCurveZeros(FinDiscountCurve):
    ''' This is a curve calculated from a set of dates and zero rates. As we
    have rates as inputs, we need to specify the corresponding compounding
    frequency. Also to go from rates and dates to discount factors we need to
    compute the year fraction correctly and for this we require a day count
    convention. Finally, we need to interpolate the zero rate for the times
    between the zero rates given and for this we must specify an interpolation
    convention. The class inherits methods from FinDiscountCurve. '''

###############################################################################

    def __init__(self,
                 valuationDate: FinDate,
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 freqType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA,
                 interpType: FinInterpTypes = FinInterpTypes.FLAT_FWD_RATES):
        ''' Create the discount curve from a vector of dates and zero rates
        factors. The first date is the curve anchor. Then a vector of zero
        dates and then another same-length vector of rates. The rate is to the
        corresponding date. We must specify the compounding frequency of the
        zero rates and also a day count convention for calculating times which
        we must do to calculate discount factors. Finally we specify the
        interpolation scheme for off-grid dates.'''

        checkArgumentTypes(self.__init__, locals())

        # Validate curve
        if len(zeroDates) == 0:
            raise FinError("Dates has zero length")

        if len(zeroDates) != len(zeroRates):
            raise FinError("Dates and Rates are not the same length")

        if freqType not in FinFrequencyTypes:
            raise FinError("Unknown Frequency type " + str(freqType))

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(dayCountType))

        self._valuationDate = valuationDate
        self._freqType = freqType
        self._dayCountType = dayCountType
        self._interpType = interpType

        self._zeroRates = np.array(zeroRates)
        self._zeroDates = zeroDates

        self._times = timesFromDates(zeroDates, valuationDate, dayCountType)

        if testMonotonicity(self._times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        dfs = self._zeroToDf(self._valuationDate,
                             self._zeroRates,
                             self._times,
                             self._freqType,
                             self._dayCountType)

        self._dfs = np.array(dfs)
        self._interpolator = FinInterpolator(self._interpType)
        self._interpolator.fit(self._times, self._dfs)

# ###############################################################################

#     def bump(self, bumpSize):
#         ''' Calculate the continuous forward rate at the forward date. '''

#         times = self._times.copy()
#         discountFactors = self._discountFactors.copy()

#         n = len(self._times)
#         for i in range(0, n):
#             t = times[i]
#             discountFactors[i] = discountFactors[i] * np.exp(-bumpSize*t)

#         discCurve = FinDiscountCurve(self._valuationDate, times,
#                                      discountFactors,
#                                      self._interpType)

#         return discCurve

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUATION DATE", self._valuationDate)
        s += labelToString("FREQUENCY TYPE", (self._freqType))
        s += labelToString("DAY COUNT TYPE", (self._dayCountType))
        s += labelToString("INTERP TYPE", (self._interpType))

        s += labelToString("DATES", "ZERO RATES")
        numPoints = len(self._times)
        for i in range(0, numPoints):
            s += labelToString("%12s" % self._zeroDates[i],
                               "%10.7f" % self._zeroRates[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################

