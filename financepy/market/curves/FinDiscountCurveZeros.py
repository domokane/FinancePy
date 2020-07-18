##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from .FinInterpolate import FinInterpMethods
from ...finutils.FinHelperFunctions import labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################

def timesFromDatesFlat(dt, valuationDate, dayCountType):

    dayCount = FinDayCount(dayCountType)

    if isinstance(dt, FinDate):
        numDates = 1
        times = [None]
        times[0] = dayCount.yearFrac(valuationDate, dt)
    elif isinstance(dt, list) and isinstance(dt[0], FinDate):
        numDates = len(dt)
        times = []
        for i in range(0, numDates):
            t = dayCount.yearFrac(valuationDate, dt[i])
            times.append(t)
    else:
        raise FinError("Discount factor must take dates.")

    times = np.array(times)
    return times

###############################################################################


class FinDiscountCurveZeros(FinDiscountCurve):
    ''' This is a curve calculated from a set of dates and zero rates. As we
    have rates as inputs, we need to specify the corresponding compounding
    frequency. Also to go from rates and dates to discount factors we need to
    compute the year fraction correctly and for this we require a day count
    convention. Finally, we need to interpolate the zero rate for the times
    between the zero rates given and for this we must specify an interpolation
    convention. '''

###############################################################################

    def __init__(self,
                 valuationDate,
                 dates,
                 zeroRates,
                 frequencyType=FinFrequencyTypes.ANNUAL,
                 dayCountType=FinDayCountTypes.ACT_ACT_ISDA,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Create the discount curve from a vector of dates and zero rates
        factors. The first date is the curve anchor. Then a vector of zero
        dates and then another same-length vector of rates. The rate is to the
        corresponding date. We must specify the compounding frequency of the
        zero rates and also a day count convention for calculating times which
        we must do to calculate discount factors. Finally we specify the
        interpolation scheme. '''

        # Validate curve
        if len(dates) == 0:
            raise FinError("Dates has zero length")

        if len(dates) != len(zeroRates):
            raise FinError("Dates and Rates are not the same length")

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Unknown Frequency type " + str(frequencyType))

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(dayCountType))

        self._valuationDate = valuationDate
        self._frequencyType = frequencyType
        self._dayCountType = dayCountType
        self._interpMethod = interpMethod

        self._times = timesFromDatesFlat(dates, valuationDate, dayCountType)
        self._zeroRates = zeroRates

        if testMonotonicity(self._times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        self._buildCurvePoints()

###############################################################################

    def _buildCurvePoints(self):
        ''' Hidden function to extract discount factors from zero rates. '''

        f = FinFrequency(self._frequencyType)
        values = []

        # We just calculate the discount factors using times as provided
        numTimes = len(self._times)

        for i in range(0, numTimes):
            t = self._times[i]
            r = self._zeroRates[i]

            if self._frequencyType == FinFrequencyTypes.CONTINUOUS:
                df = np.exp(-r*t)
            else:
                df = 1.0 / np.power(1.0 + r/f, f * t)

            values.append(df)

        self._discountFactors = np.array(values)

###############################################################################

    def bump(self, bumpSize):
        ''' Calculate the continuous forward rate at the forward date. '''

        times = self._times.copy()
        discountFactors = self._discountFactors.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            discountFactors[i] = discountFactors[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._curveDate,
                                     times,
                                     discountFactors,
                                     self._interpMethod)

        return discCurve

###############################################################################

    def __repr__(self):

        numPoints = len(self._times)
        s = labelToString("TIMES", "DISCOUNT FACTORS")
        for i in range(0, numPoints):
            s += labelToString(self._times[i], self._values[i])

        return s

###############################################################################
