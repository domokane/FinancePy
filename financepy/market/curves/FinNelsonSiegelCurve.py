# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np
from math import exp

from ...finutils.FinGlobalVariables import gSmall
from ...market.curves.FinCurve import FinCurve

##########################################################################
# TODO: Inherit from FinCurve and add df method
# TODO: Add fitting optimiser to take in bonds and do a best fit
##########################################################################


##########################################################################
##########################################################################

class FinNelsonSiegelCurve(FinCurve):
    ''' Implementation of Nelson-Siegel parametrisation of the zero rate curve '''

    def __init__(self, beta1, beta2, beta3, tau):

        if tau <= 0:
            raise ValueError("Tau must be positive")

        self._beta1 = beta1
        self._beta2 = beta2
        self._beta3 = beta3
        self._tau = tau

##########################################################################

    def zero(self, t):
        ''' Calculation of zero rates. This function can return a vector
        of zero rates given a vector of times. '''

        if np.any(t < 0.0):
            raise ValueError("All times must be positive")

        t = t + gSmall  # To avoid overflows when t=0.0
        theta = t / self._tau
        expTerm = np.exp(-theta)
        zeroRate = self._beta1
        zeroRate += self._beta2 * (1.0 - expTerm) / theta
        zeroRate += self._beta3 * ((1.0 - expTerm) / theta - expTerm)
        return zeroRate

##########################################################################

    def fwd(self, t):
        ''' Calculation of forward rates. This function can return a vector
        of forward rates given a vector of times. '''

        theta = t / self._tau
        expTerm = np.exp(-theta)
        fwdRate = self.beta1
        fwdRate += self.beta2 * expTerm
        fwdRate += self.beta3 * theta * expTerm
        return fwdRate

##########################################################################

    def df(self, t):
        ''' Discount factor for Nelson-Siegel curve parametrisation. '''
        r = self.zero(t)
        return exp(-r * t)

#############################################################################
