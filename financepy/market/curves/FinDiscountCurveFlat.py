###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import numpy as np

###############################################################################

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCountTypes, FinDayCount
from ...finutils.FinFrequency import FinFrequencyTypes, FinFrequency
from ...finutils.FinHelperFunctions import inputTime
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################
# TODO: Do I need to add a day count to ensure rate and times are linked in
#       the correct way URGENT
###############################################################################


class FinDiscountCurveFlat(FinDiscountCurve):
    ''' A trivially simple discount curve based on a single zero rate with its
    own specified compounding method. Hence the curve is assumed to be flat.
    It is used for quick and dirty analysis and when limited information is
    available. '''

###############################################################################

    def __init__(self,
                 curveDate: FinDate,
                 flatRate: float,
                 rateFrequencyType: FinFrequencyTypes=FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA):
        ''' Create a discount curve which is flat. This is very useful for
        quick testing and simply requires a curve date and a rate and also a
        frequency. As we have entered a rate, a corresponding day count
        convention must be used to specify how time periods are to be measured.
        '''

        checkArgumentTypes(self.__init__, locals())

        self._curveDate = curveDate
        self._flatRate = flatRate
        self._rateFrequencyType = rateFrequencyType
        self._rateFrequency = FinFrequency(rateFrequencyType)
        self._dayCountType = dayCountType
        self._dayCount = FinDayCount(dayCountType)

        # Some methods require these members but don't do this if rates is
        # an array of doubles else it will fail due to vectorisations of
        # different lengths

        self._times = None
        self._values = None

        if isinstance(self._flatRate, np.ndarray) is False:
            self._times = np.linspace(0.0, 10.0, 40)
            self._values = self.df(self._times)

###############################################################################

    # def zeroRate(self,
    #              dt: FinDate,
    #              zeroRateFrequencyType):
    #     ''' Return the zero rate which is simply the curve rate. '''

    #     t = inputTime(dt, self)

    #     if isinstance(dt, FinDate):
    #         t = self._dayCount.yearFrac(self._curveDate, dt)

    #     if zeroRateFrequencyType == self._rateFrequencyType:
    #         return self._rate

    #     if self._rateFrequencyType == FinFrequencyTypes.CONTINUOUS:
    #         df = np.exp(-self._rate * t)
    #     else:
    #         df = (1.0 + self._rate/self._cmpdFreq)**(-t*self._cmpdFreq)

    #     if zeroRateFrequencyType == FinFrequencyTypes.CONTINUOUS:
    #         r = -np.log(df) / t
    #     else:
    #         f = FinFrequency(zeroRateFrequencyType)
    #         r = (df**(-1.0/t) - 1) * f
    #     return r

###############################################################################

    def bump(self, bumpSize):
        ''' Calculate the continuous forward rate at the forward date. '''
        rBumped = self._flatRate + bumpSize
        discCurve = FinDiscountCurveFlat(self._curveDate, rBumped,
                                         compoundingFreq=self._cmpdFreq)
        return discCurve

###############################################################################

    def df(self, dt):
        ''' Return the discount factor based on the compounding approach. '''

        if isinstance(dt, FinDate):
            t = self._dayCount.yearFrac(self._curveDate, dt)
        else:
            t = dt

        r = self._flatRate
        if self._rateFrequencyType == FinFrequencyTypes.CONTINUOUS:
            df = np.exp(-r * t)
        else:
            df = (1 + r/self._rateFrequency)**(-t*self._rateFrequency)

        return df

###############################################################################

    def __repr__(self):

        s = labelToString("Flat rate:%9.5f" % (self._flatRate))
        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
