##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from math import exp

from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinHelperFunctions import inputTime
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinError import FinError
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################


class FinDiscountCurveNS(FinDiscountCurve):
    ''' Implementation of Nelson-Siegel parametrisation of a discount curve.
    The internal rate is a continuously compounded rate but you can calculate
    alternative frequencies by providing a corresponding compounding frequency.
    '''

    def __init__(self,
                 curveDate,
                 params):
        ''' Creation of a Nelson-Siegel curve. Parameters are provided as a
        list or vector of 4 values for beta1, beta2, beta3 and tau. '''

        self._curveDate = curveDate

        if len(params) != 4:
            raise FinError("NS requires a vector/list of 4 parameters.")

        if params[3] <= 0:
            raise FinError("Tau must be positive")

        self._beta1 = params[0]
        self._beta2 = params[1]
        self._beta3 = params[2]
        self._tau = params[3]

###############################################################################

    def zeroRate(self, dt, frequencyType=FinFrequencyTypes.CONTINUOUS):
        ''' Calculation of zero rates with specified frequency. This
        function can return a vector of zero rates given a vector of
        times so must use Numpy functions. Default frequency is a
        continuously compounded rate. '''

        t = inputTime(dt, self)
        f = FinFrequency(frequencyType)
        t = t + gSmall  # To avoid overflows if t=0.0

        theta = t / self._tau
        e = np.exp(-theta)
        zeroRate = self._beta1 + self._beta2 * (1.0 - e) / theta
        zeroRate += self._beta3 * ((1.0 - e) / theta - e)
        df = np.exp(-zeroRate * t)

        if f == 0:  # Simple interest
            r = (1.0/df-1.0)/t
        if f == -1:  # Continuous
            r = zeroRate
        else:
            r = (df**(-1.0/t)-1.0) * f

        return r

###############################################################################

    def fwd(self, dt):
        ''' Calculation of continuously compounded forward rates. This
        function can return a vector of instantaneous forward rates given a
        vector of times. '''
        t = inputTime(dt, self)
        theta = t / self._tau
        e = np.exp(-theta)
        fwdRate = self.beta1 + self.beta2 * e + self.beta3 * theta * e
        return fwdRate

###############################################################################

    def df(self, dt):
        ''' Discount factor for Nelson-Siegel curve parametrisation. '''
        t = inputTime(dt, self)
        r = self.zeroRate(t)
        return exp(-r * t)

###############################################################################
