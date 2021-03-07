##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import gSmall
from ...utils.math import testMonotonicity
from ...utils.frequency import FrequencyTypes
from ...utils.helpers import labelToString
from ...utils.helpers import check_argument_types
from ...utils.day_count import DayCountTypes
from ...utils.helpers import timesFromDates
from ...market.discount.curve import DiscountCurve

###############################################################################


class DiscountCurvePWL(DiscountCurve):
    """ Curve is made up of a series of sections assumed to each have a
    piece-wise linear zero rate. The zero rate has a specified frequency
    which defaults to continuous. This curve inherits all of the extra methods
    from FinDiscountCurve. """

    def __init__(self,
                 valuation_date: Date,
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 freq_type: FrequencyTypes = FrequencyTypes.CONTINUOUS,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_ACT_ISDA):
        """ Curve is defined by a vector of increasing times and zero rates."""

        check_argument_types(self.__init__, locals())

        self._valuation_date = valuation_date

        if len(zeroDates) != len(zeroRates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zeroDates) == 0:
            raise FinError("Dates vector must have length > 0")

        self._zeroRates = np.array(zeroRates)
        self._zeroDates = zeroDates
        self._freq_type = freq_type
        self._day_count_type = day_count_type

        dcTimes = timesFromDates(zeroDates,
                                 self._valuation_date,
                                 self._day_count_type)

        self._times = np.array(dcTimes)

        if testMonotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

###############################################################################

    def _zeroRate(self,
                  times: (list, np.ndarray)):
        """ Calculate the piecewise linear zero rate. This is taken from the
        initial inputs. A simple linear interpolation scheme is used. If the
        user supplies a frequency type then a conversion is done. """

        if isinstance(times, float):
            times = np.array([times])

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, 1e-6)

        zeroRates = []

        for t in times:
            l_index = 0
            found = 0

            numTimes = len(self._times)
            for i in range(1, numTimes):
                if self._times[i] > t:
                    l_index = i - 1
                    found = 1
                    break

            t0 = self._times[l_index]
            r0 = self._zeroRates[l_index]
            t1 = self._times[l_index+1]
            r1 = self._zeroRates[l_index+1]

            if found == 1:
                zeroRate = ((t1 - t) * r0 + (t - t0) * r1)/(t1 - t0)
            else:
                zeroRate = self._zeroRates[-1]

            zeroRates.append(zeroRate)

        return np.array(zeroRates)

###############################################################################

    def df(self,
           dates: (Date, list)):
        """ Return discount factors given a single or vector of dates. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA. """

        # Get day count times to use with curve day count convention
        dcTimes = timesFromDates(dates,
                                 self._valuation_date,
                                 self._day_count_type)

        zeroRates = self._zeroRate(dcTimes)

        df = self._zeroToDf(self._valuation_date,
                            zeroRates,
                            dcTimes,
                            self._freq_type,
                            self._day_count_type)

        return df

###############################################################################

    # def _df(self,
    #         t: (float, np.ndarray)):
    #     """ Returns the discount factor at time t taking into account the
    #     piecewise flat zero rate curve and the compunding frequency. """

    #     r = self._zeroRate(t, self._freq_type)
    #     df = zeroToDf(r, t, self._freq_type)
    #     return df

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("DATE", "ZERO RATE")
        for i in range(0, len(self._zeroDates)):
            s += labelToString(self._zeroDates[i], self._zeroRates[i])
        s += labelToString("FREQUENCY", (self._freq_type))
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
