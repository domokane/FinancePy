##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinFrequency import zeroToDf

###############################################################################


class FinDiscountCurvePoly(FinDiscountCurve):
    ''' Zero Rate Curve of a specified frequency parametrised using a cubic
    polynomial. The zero rate is assumed to be continuously compounded but
    this can be amended by providing a frequency when extracting zero rates. 
    The class inherits all of the methods from FinDiscountCurve. '''

    def __init__(self,
                 valuationDate: FinDate,
                 coefficients: (list, np.ndarray),
                 frequencyType: FinFrequencyTypes=FinFrequencyTypes.CONTINUOUS):
        ''' Create zero rate curve parametrised using a cubic curve from
        coefficients and specifying a compounding frequency type. '''

        checkArgumentTypes(self.__init__, locals())

        self._valuationDate = valuationDate
        self._coefficients = coefficients
        self._power = len(coefficients) - 1
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
# TODO : NEED TO CONVERT ZERO RATE TO REQUESTED FREQUENCY

    def _zeroRate(self,
                  times: (float, np.ndarray),
                  frequencyType: FinFrequencyTypes = None):
        ''' Calculate the zero rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use. The compounding frequency defaults to that specified in the
        constructor of the curve object. '''

        if frequencyType is None:
            frequencyType = self._frequencyType

        if isinstance(times, float):
            times = np.array([times])

        if np.any(times < 0.0):
            raise FinError("All times must be positive")

        times = np.maximum(times, 1e-6)

        zerosList = []

        for t in times:
            zeroRate = 0.0
            for n in range(0, len(self._coefficients)):
                zeroRate += self._coefficients[n] * (t**n)
            zerosList.append(zeroRate)
        return np.array(zerosList)

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
        ''' Given a time or list of times this returns the corresponding
        vector of discount factors. '''

        r = self._zeroRate(t)
        df = zeroToDf(r, t, self._frequencyType)
        return df

###############################################################################

    # def _fwd(self,
    #          t: float):
    #     ''' Returns the continuously compounded forward rate at time t. '''

    #     dt = 1e-5
    #     df1 = self._df(t)
    #     df2 = self._df(t+dt)
    #     fwdRate = -(df2 - df1)/dt
    #     return fwdRate

###############################################################################

    def __repr__(self):
        ''' Display internal parameters of curve. '''
        s = type(self).__name__ + "\n"
        s += labelToString("POWER", "COEFFICIENT")
        for i in range(0, len(self._coefficients)):
            s += labelToString(str(i), self._coefficients[i])
        s += labelToString("FREQUENCY", (self._frequencyType))

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
