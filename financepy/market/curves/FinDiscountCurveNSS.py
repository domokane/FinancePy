##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinError import FinError
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes

###############################################################################


class FinDiscountCurveNSS(FinDiscountCurve):
    ''' Implementation of Nelson-Siegel-Svensson parametrisation of the
    zero rate curve. The zero rate is assumed to be continuously compounded.
    This can be changed when calling for zero rates. The class inherits lots
    of methods from FinDiscountCurve. '''

    def __init__(self,
                 valuationDate: FinDate,
                 beta0: float,
                 beta1: float,
                 beta2: float,
                 beta3: float,
                 tau1: float,
                 tau2: float,
                 frequencyType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS):
        ''' Create a FinDiscountCurveNSS object by passing in curve valuation
        date plus the 4 different beta values and the 2 tau values. The zero
        rates produced by this parametrisation have an implicit compounding
        convention that defaults to continuous but can be overriden. '''

        checkArgumentTypes(self.__init__, locals())

        if tau1 <= 0:
            raise FinError("Tau1 must be positive")

        if tau2 <= 0:
            raise FinError("Tau2 must be positive")

        self._valuationDate = valuationDate
        self._beta0 = beta0
        self._beta1 = beta1
        self._beta2 = beta2
        self._beta3 = beta3
        self._tau1 = tau1
        self._tau2 = tau2
        self._frequencyType = frequencyType

###############################################################################

    # def zeroRate(self,
    #              dt: (list, FinDate),
    #              frequencyType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS):
    #     ''' Calculation of zero rates with specified frequency. This
    #     function can return a vector of zero rates given a vector of
    #     times so must use Numpy functions. Default frequency is a
    #     continuously compounded rate. '''

    #     times = timesFromDates(dt, self._valuationDate)
    #     zeroRates = self._zeroRate(times, frequencyType)

    #     if frequencyType != self._frequencyType:
    #         dfs = zeroToDf(zeroRates, self._frequencyType)
    #         zeroRates = dfToZero(dfs, frequencyType)

    #     if isinstance(dt, FinDate):
    #         return zeroRates[0]
    #     else:
    #         return np.array(zeroRates)

###############################################################################

    def _zeroRate(self,
                  times: (float, np.ndarray),
                  frequencyType: FinFrequencyTypes):
        ''' Calculation of zero rates given a single time or a numpy vector of
        times. This function can return a single zero rate or a vector of zero
        rates. The compounding frequency must be provided. '''

        if isinstance(times, float):
            times = np.array([times])

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, 1e-6)

        theta1 = times / self._tau1
        theta2 = times / self._tau2
        e1 = np.exp(-theta1)
        e2 = np.exp(-theta2)
        zeroRate = self._beta0
        zeroRate += self._beta1 * (1.0 - e1) / theta1
        zeroRate += self._beta2 * ((1.0 - e1) / theta1 - e1)
        zeroRate += self._beta3 * ((1.0 - e2) / theta2 - e2)
        return zeroRate

##########################################################################

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
        ''' Discount factor for Nelson-Siegel-Svensson curve
        parametrisation. '''
        r = self._zeroRate(t, self._frequencyType)
        return np.exp(-r * t)

###############################################################################

    # def fwd(self,
    #         dt: FinDate):
    #     ''' Calculate the continuously compounded forward rate at the forward
    #     FinDate provided. This is done by perturbing the time by a small amount
    #     and measuring the change in the log of the discount factor divided by
    #     the time increment dt.'''

    #     times = timesFromDates(dt, self._valuationDate)
    #     fwds = self._fwd(times)

    #     if isinstance(dt, FinDate):
    #         return fwds
    #     else:
    #         return np.array(fwds)

##########################################################################

    # def _fwd(self,
    #          times: (np.ndarray, list)):
    #     ''' Calculate the continuously compounded forward rate at the forward
    #     time provided. This is done by perturbing the time by a small amount
    #     and measuring the change in the log of the discount factor divided by
    #     the time increment dt.'''

    #     dt = 1e-6
    #     times = np.maximum(times, dt)
    #     df1 = self._df(times-dt)
    #     df2 = self._df(times+dt)
    #     fwd = np.log(df1/df2)/(2.0*dt)
    #     return fwd

###############################################################################

    def __repr__(self):
        s = type(self).__name__ + "\n"
        s += labelToString("PARAMETER", "VALUE")
        s += labelToString("BETA0", self._beta0)
        s += labelToString("BETA1", self._beta1)
        s += labelToString("BETA2", self._beta2)
        s += labelToString("BETA3", self._beta3)
        s += labelToString("TAU1", self._tau1)
        s += labelToString("TAU2", self._tau2)
        s += labelToString("FREQUENCY", (self._frequencyType))
        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
