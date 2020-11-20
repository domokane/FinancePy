##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinHelperFunctions import timesFromDates
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################


class FinDiscountCurvePWF(FinDiscountCurve):
    ''' Curve is made up of a series of zero rates sections with each having
    a piecewise flat zero rate. The default compounding assumption is
    continuous. The class inherits methods from FinDiscountCurve. '''

    def __init__(self,
                 valuationDate: FinDate,
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 freqType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA):
        ''' Creates a discount curve using a vector of times and zero rates
        that assumes that the zero rates are piecewise flat. '''

        checkArgumentTypes(self.__init__, locals())

        self._valuationDate = valuationDate

        if len(zeroDates) != len(zeroRates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zeroDates) == 0:
            raise FinError("Dates vector must have length > 0")

        self._zeroDates = zeroDates
        self._zeroRates = np.array(zeroRates)
        self._freqType = freqType
        self._dayCountType = dayCountType

        dcTimes = timesFromDates(zeroDates,
                                 self._valuationDate,
                                 self._dayCountType)

        self._times = np.array(dcTimes)

        if testMonotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

###############################################################################

    def _zeroRate(self,
                  times: (float, np.ndarray, list)):
        ''' The piecewise flat zero rate is selected and returned. '''

        if isinstance(times, float):
            times = np.array([times])

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, gSmall)

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

            r0 = self._zeroRates[l_index]

            if found == 1:
                zeroRate = r0
            else:
                zeroRate = self._zeroRates[-1]

            zeroRates.append(zeroRate)

        return np.array(zeroRates)

###############################################################################

    def _fwd(self,
             times: (np.ndarray, list)):
        ''' Calculate the continuously compounded forward rate at the forward
        time provided. This is done by perturbing the time by a small amount
        and measuring the change in the log of the discount factor divided by
        the time increment dt.'''

        dt = 1e-6
        times = np.maximum(times, dt)

        df1 = self._df(times-dt)
        df2 = self._df(times+dt)
        fwd = np.log(df1/df2)/(2.0*dt)

        # Prevent the discontinuous spikes by capping them at the sum of short
        # plus the long end zero rate !
        cap = (self._zeroRates[0] + self._zeroRates[-1])
        fwd = np.minimum(fwd, cap)

        return fwd

###############################################################################

    def df(self,
           dates: (FinDate, list)):
        ''' Return discount factors given a single or vector of dates. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA. '''

        # Get day count times to use with curve day count convention
        dcTimes = timesFromDates(dates,
                                 self._valuationDate,
                                 self._dayCountType)

        zeroRates = self._zeroRate(dcTimes)

        df = self._zeroToDf(self._valuationDate,
                            zeroRates,
                            dcTimes,
                            self._freqType,
                            self._dayCountType)

        return df

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("DATE", "ZERO RATE")
        for i in range(0, len(self._zeroDates)):
            s += labelToString(self._zeroDates[i], self._zeroRates[i])
        s += labelToString("FREQUENCY", (self._freqType))
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
