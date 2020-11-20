##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinError import FinError
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinHelperFunctions import timesFromDates

###############################################################################


class FinDiscountCurveNSS(FinDiscountCurve):
    ''' Implementation of Nelson-Siegel-Svensson parametrisation of the
    zero rate curve. The zero rate is assumed to be continuously compounded.
    This can be changed when calling for zero rates. A day count convention is
    needed to ensure that dates are converted to the correct time in years. The
    class inherits methods from FinDiscountCurve.'''

    def __init__(self,
                 valuationDate: FinDate,
                 beta0: float,
                 beta1: float,
                 beta2: float,
                 beta3: float,
                 tau1: float,
                 tau2: float,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA):
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
        self._freqType = freqType
        self._dayCountType = dayCountType

###############################################################################

    def zeroRate(self,
                 dates: (list, FinDate),
                 freqType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360):
        ''' Calculation of zero rates with specified frequency according to
        NSS parametrisation. This method overrides that in FinDiscountCurve.
        The NSS parametrisation is no strictly terms of continuously compounded
        zero rates, this function allows other compounding and day counts.
        This function returns a single or vector of zero rates given a vector
        of dates so must use Numpy functions. The default frequency is a
        continuously compounded rate and ACT ACT day counting. '''

        if isinstance(freqType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if isinstance(dayCountType, FinDayCountTypes) is False:
            raise FinError("Invalid Day Count type.")

        # Get day count times to use with curve day count convention
        dcTimes = timesFromDates(dates,
                                 self._valuationDate,
                                 self._dayCountType)

        # We now get the discount factors using these times
        zeroRates = self._zeroRate(dcTimes)

        # Now get the discount factors using curve conventions
        dfs = self._zeroToDf(self._valuationDate,
                             zeroRates,
                             dcTimes,
                             self._freqType,
                             self._dayCountType)

        # Convert these to zero rates in the required frequency and day count
        zeroRates = self._dfToZero(dfs,
                                   dates,
                                   freqType,
                                   dayCountType)

        if isinstance(dates, FinDate):
            return zeroRates[0]
        else:
            return np.array(zeroRates)

        return zeroRates

###############################################################################

    def _zeroRate(self,
                  times: (float, np.ndarray)):
        ''' Calculation of zero rates given a single time or a numpy vector of
        times. This function can return a single zero rate or a vector of zero
        rates. The compounding frequency must be provided. '''

        t = np.maximum(times, gSmall)

        theta1 = t / self._tau1
        theta2 = t / self._tau2
        e1 = np.exp(-theta1)
        e2 = np.exp(-theta2)
        zeroRate = self._beta0
        zeroRate += self._beta1 * (1.0 - e1) / theta1
        zeroRate += self._beta2 * ((1.0 - e1) / theta1 - e1)
        zeroRate += self._beta3 * ((1.0 - e2) / theta2 - e2)
        return zeroRate

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

        if isinstance(dates, FinDate):
            return df[0]
        else:
            return df

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("PARAMETER", "VALUE")
        s += labelToString("BETA0", self._beta0)
        s += labelToString("BETA1", self._beta1)
        s += labelToString("BETA2", self._beta2)
        s += labelToString("BETA3", self._beta3)
        s += labelToString("TAU1", self._tau1)
        s += labelToString("TAU2", self._tau2)
        s += labelToString("FREQUENCY", (self._freqType))
        s += labelToString("DAY_COUNT", (self._dayCountType))
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
