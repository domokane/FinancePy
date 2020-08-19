##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinMath import testMonotonicity
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinFrequency import zeroToDf
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinHelperFunctions import labelToString

###############################################################################


class FinDiscountCurvePWF(FinDiscountCurve):
    ''' Curve is made up of a series of zero rates sections with each having
    a piecewise flat zero rate. The default compounding assumption is
    continuous. The class inherits methods from FinDiscountCurve. '''

    def __init__(self,
                 valuationDate: FinDate,
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 frequencyType: FinFrequencyTypes=FinFrequencyTypes.CONTINUOUS):
        ''' Creates a discount curve using a vector of times and zero rates
        that assumes that the zero rates are piecewise flat. '''

        checkArgumentTypes(self.__init__, locals())

        self._valuationDate = valuationDate

        if len(zeroDates) != len(zeroRates):
            raise FinError("Dates and rates vectors must have same length")

        if len(zeroDates) == 0:
            raise FinError("Dates vector must have length > 0")

        self._times = []
        for dt in zeroDates:
            t = (dt - self._valuationDate) / gDaysInYear
            self._times.append(t)

#        if self._times[0] > 0.0:
#            self._times.insert(0, 0.0)
#            r0 = zeroRates[0]
#            zeroRates = zeroRates.tolist()
#            zeroRates.insert(0, r0)

        self._times = np.array(self._times)
        self._zeroRates = np.array(zeroRates)
        self._zeroDates = zeroDates

        if testMonotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._frequencyType = frequencyType

###############################################################################

    # def zeroRate(self,
    #              dt: (list, FinDate),
    #              frequencyType=FinFrequencyTypes.CONTINUOUS):
    #     ''' Calculate the zero rate to maturity date. The compounding frequency
    #     of the rate defaults to continuous which is useful for supplying a rate
    #     to theoretical models such as Black-Scholes that require a continuously
    #     compounded zero rate as input. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     zeroRates = self._zeroRate(times, frequencyType)

    #     if isinstance(dt, FinDate):
    #         print(dt)
    #         print(zeroRates)
    #         return zeroRates[0]
    #     else:
    #         return np.array(zeroRates)

###############################################################################

    def _zeroRate(self,
                  times: (float, np.ndarray, list),
                  frequencyType=FinFrequencyTypes.CONTINUOUS):
        ''' The piecewise flat zero rate is selected and returned. '''

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

            r0 = self._zeroRates[l_index]

            if found == 1:
                zeroRate = r0
            else:
                zeroRate = self._zeroRates[-1]

            zeroRates.append(zeroRate)

        return np.array(zeroRates)

###############################################################################

    # def fwd(self,
    #         dt: FinDate):
    #     ''' Returns the continuously compounded forward rate at time t. As the
    #     curve is piecewise flat in zero rate space, its shape in forward rates
    #     requires a numerical calculation of the fwd rate. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     fwds = self._fwd(times)

    #     if isinstance(dt, FinDate):
    #         return fwds
    #     else:
    #         return np.array(fwds)

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

    # def df(self,
    #        dt):
    #     ''' Function to calculate a discount factor from a date or a
    #     vector of dates. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     dfs = self._df(times)

    #     if isinstance(dt, FinDate):
    #         return dfs[0]
    #     else:
    #         return np.array(dfs)

###############################################################################

    def _df(self,
            t: (float, np.ndarray)):
        ''' Given a time or vector of times this returns the corresponding
        vector of discount factors. '''

        r = self._zeroRate(t)
        df = zeroToDf(r, t, self._frequencyType)
        return df

###############################################################################

    def __repr__(self):
        s = type(self).__name__ + "\n"
        s += labelToString("DATE", "ZERO RATE")
        for i in range(0, len(self._zeroDates)):
            s += labelToString(self._zeroDates[i], self._zeroRates[i])
        s += labelToString("FREQUENCY", (self._frequencyType))
        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
