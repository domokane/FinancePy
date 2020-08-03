##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

###############################################################################
# TODO
# Inherit from FinCurve and add df method
# Put in a convention for the rate
# Use Frequency object
###############################################################################

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinFrequency import zeroToDf, dfToZero
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString


###############################################################################


class FinDiscountCurvePWL(FinDiscountCurve):
    ''' Curve is made up of a series of sections assumed to each have a
    piece-wise linear zero rate. The zero rate has a specified frequency
    which defaults to continuous. This curve inherits all of the extra methods
    from FinDiscountCurve. '''

    def __init__(self,
                 valuationDate: FinDate,
                 zeroDates: list,
                 zeroRates: (list, np.ndarray),
                 frequencyType: FinFrequencyTypes= FinFrequencyTypes.CONTINUOUS):
        ''' Curve is defined by a vector of increasing times and zero rates.'''

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

        print(zeroDates)
        print(zeroRates)

        self._times = np.array(self._times)
        self._zeroRates = np.array(zeroRates)
        self._zeroDates = zeroDates

        if testMonotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._frequencyType = frequencyType

###############################################################################

    # def zeroRate(self,
    #              dt: FinDate,
    #              frequencyType=FinFrequencyTypes.CONTINUOUS):
    #     ''' Calculate the zero rate to maturity date. The compounding frequency
    #     of the rate defaults to continuous which is useful for supplying a rate
    #     to theoretical models such as Black-Scholes that require a continuously
    #     compounded zero rate as input. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     zeroRates = self._zeroRate(times, frequencyType)

    #     if isinstance(dt, FinDate):
    #         return zeroRates[0]
    #     else:
    #         return np.array(zeroRates)

###############################################################################

    def _zeroRate(self,
                  times: (list, np.ndarray),
                  frequencyType: FinFrequencyTypes):
        ''' Calculate the piecewise linear zero rate. This is taken from the
        initial inputs. A simple linear interpolation scheme is used. If the
        user supplies a frequency type then a conversion is done. '''

        if isinstance(times, float):
            times = np.array([times])

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

    # def parRate(self,
    #             dt: FinDate,
    #             frequencyType=FinFrequencyTypes.ANNUAL):
    #     ''' Calculate the par rate to maturity date. This is the rate paid by a
    #     bond that has a price of par today. For a par swap rate, use the swap
    #     rate function that takes in the swap details. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     parRates = self._parRate(times, frequencyType)

    #     if isinstance(dt, FinDate):
    #         return parRates[0]
    #     else:
    #         return np.array(parRates)

##########################################################################

    # def _parRate(self,
    #              times: list,
    #              frequencyType=FinFrequencyTypes.ANNUAL):
    #     ''' Calculate the zero rate to maturity date but with times as inputs.
    #     This function is used internally and should be discouraged for external
    #     use. The compounding frequency defaults to continuous. '''

    #     if frequencyType == FinFrequencyTypes.SIMPLE:
    #         raise FinError("Cannot calculate par rate with simple yield freq.")
    #     elif frequencyType == FinFrequencyTypes.CONTINUOUS:
    #         raise FinError("Cannot calculate par rate with continuous freq.")

    #     f = FinFrequency(frequencyType)
    #     dt = 1.0 / f

    #     if np.any(times < 0.0):
    #         raise FinError("Unable to calculate par rate as one t < 0.0")

    #     parRateList = []

    #     for t in times:

    #         flowTimes = pv01Times(t, f)
    #         pv01 = 0.0
    #         for tFlow in flowTimes:
    #             pv01 += self._df(tFlow)

    #         pv01 = pv01 * dt
    #         dft = self._df(t)
    #         parRate = (1.0 - dft) / pv01
    #         parRateList.append(parRate)

    #     return np.array(parRateList)

###############################################################################

    # def _fwd(self,
    #          times: (np.ndarray, float)):
    #     ''' Returns the continuously compounded forward rate at time t. As the
    #     curve is piecewise linear in zero rate space, its shape in forward 
    #     rates requires a numerical calculation of the fwd rate. '''

    #     dt = 1e-6
    #     times = np.maximum(times, dt)

    #     df1 = self._df(times-dt)
    #     df2 = self._df(times+dt)
    #     fwd = np.log(df1/df2)/(2.0*dt)
    #     return fwd

###############################################################################

    # def df(self,
    #        dt: (list, FinDate)):
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
        ''' Returns the discount factor at time t taking into account the
        piecewise flat zero rate curve and the compunding frequency. '''

        r = self._zeroRate(t, self._frequencyType)
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
