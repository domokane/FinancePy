###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import numpy as np

###############################################################################

from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinFrequency import zeroToDf

###############################################################################
# TODO: Do I need to add a day count to ensure rate and times are linked in
#       the correct way URGENT
###############################################################################


class FinDiscountCurveFlat(FinDiscountCurve):
    ''' A very simple discount curve based on a single zero rate with its
    own specified compounding method. Hence the curve is assumed to be flat.
    It is used for quick and dirty analysis and when limited information is
    available. It inherits several methods from FinDiscountCurve. '''

###############################################################################

    def __init__(self,
                 valuationDate: FinDate,
                 flatRate: float,
                 frequencyType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA):
        ''' Create a discount curve which is flat. This is very useful for
        quick testing and simply requires a curve date and a rate and also a
        frequency. As we have entered a rate, a corresponding day count
        convention must be used to specify how time periods are to be measured.
        As the curve is flat, no interpolation scheme is required.
        '''
        checkArgumentTypes(self.__init__, locals())

        self._valuationDate = valuationDate
        self._flatRate = flatRate
        self._frequencyType = frequencyType
        self._dayCountType = dayCountType

        # Set up a grid of times and discount factors for functions that
        # need a simple curve representation in order to use NUMBA
        self._times = np.linspace(0.0, 10.0, 121)
        self._dfValues = self._df(self._times)

###############################################################################

    def bump(self,
             bumpSize: float):
        ''' Creates a new FinDiscountCurveFlat object with the entire curve
        bumped up by the bumpsize. All other parameters are preserved.'''

        rBumped = self._flatRate + bumpSize
        discCurve = FinDiscountCurveFlat(self._valuationDate,
                                         rBumped,
                                         frequencyType=self._frequencyType,
                                         dayCountType=self._dayCountType)
        return discCurve

###############################################################################

    # def df(self,
    #        dt: FinDate):
    #     ''' Function to calculate a discount factor from a date or a
    #     vector of dates. Times are calculated according to a specified
    #     convention. '''

    #     times = timesFromDates(dt, self._valuationDate, self._dayCountType)
    #     dfs = self._df(times)

    #     if isinstance(dt, FinDate):
    #         return dfs[0]
    #     else:
    #         return dfs

###############################################################################

    def _df(self,
            t: (np.ndarray, float)):
        ''' Return the discount factor given a single or vector of times. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA. '''

        r = self._flatRate
        df = zeroToDf(r, t, self._frequencyType)
        return df

###############################################################################

    def __repr__(self):
        s = type(self).__name__ + "\n"
        s += labelToString("FLAT RATE", (self._flatRate))
        s += labelToString("FREQUENCY", (self._frequencyType))
        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
