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
from ...finutils.FinGlobalVariables import gDaysInYear


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################


class FinDiscountCurveZeros(FinDiscountCurve):
    ''' This is a curve calculated from a set of times and zero rates.
    '''

###############################################################################

    def __init__(self, curveDate,
                 timesOrDates,
                 zeroRates,
                 frequencyType=FinFrequencyTypes.ANNUAL,
                 dayCountType=FinDayCountTypes.ACT_ACT_ISDA,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Create the discount curve from a vector of times and discount
        factors. First date is the curve anchor and first rates should be zero
        as it starts and ends on that date and so has no impact. '''

        # Validate curve
        if len(timesOrDates) == 0:
            raise FinError("Times has zero length")

        if len(timesOrDates) != len(zeroRates):
            raise FinError("Times and Values are not the same")

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Unknown Frequency type " + str(frequencyType))

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(dayCountType))

        freq = FinFrequency(frequencyType)

        times = []
        values = []

        if isinstance(timesOrDates[0], float):

            # We just calculate the discount factors using times as provided
            for i in range(0, len(timesOrDates)):
                t = timesOrDates[i]
                if t < 0.0:
                    raise FinError("Times must be > 0.")

                r = zeroRates[i]

                if freq == -1:
                    df = np.exp(-r*t)
                else:
                    df = 1.0 / np.power(1.0 + r/freq, freq * t)

                times.append(t)
                values.append(df)

        elif isinstance(timesOrDates[0], FinDate):

            # Now extract discount factors which depend on the zero rates which
            # have been quoted witha frequency and day count convention

            dc = FinDayCount(dayCountType)
            for i in range(0, len(timesOrDates)):
                t = (timesOrDates[i] - curveDate) / gDaysInYear
                if t < 0.0:
                    raise FinError("Times must be > 0.")

                alpha = dc.yearFrac(curveDate, timesOrDates[i])
                r = zeroRates[i]

                if freq == -1:
                    df = np.exp(-r*t)
                else:
                    df = 1.0 / np.power(1.0 + r/freq, freq * alpha)

                times.append(t)
                values.append(df)
        else:
            raise FinError("Input timeOrDates must be list of times or dates")

        times = np.array(times)

        if testMonotonicity(times) is False:
            raise FinError("Times or dates are not sorted in increasing order")

        self._curveDate = curveDate
        self._times = times
        self._values = np.array(values)
        self._dayCountType = dayCountType
        self._frequencyType = frequencyType
        self._interpMethod = interpMethod

###############################################################################

    def bump(self, bumpSize):
        ''' Calculate the continuous forward rate at the forward date. '''

        times = self._times.copy()
        values = self._values.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            values[i] = values[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._curveDate,
                                     times,
                                     values,
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
