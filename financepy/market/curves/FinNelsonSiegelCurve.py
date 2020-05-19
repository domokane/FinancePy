# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np
from math import exp

from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinHelperFunctions import inputTime, inputFrequency
from ...finutils.FinHelperFunctions import labelToString

##########################################################################


class FinNelsonSiegelCurve():
    ''' Implementation of Nelson-Siegel parametrisation of a rate curve. 
    The default is a continuously compounded rate but you can override 
    this by providing a corresponding compounding frequency. '''

    def __init__(self, curveDate, params, cmpdFreq=-1):
        ''' Creation of a Nelson-Siegel curve. Parameters are provided as a 
        list or vector of 4 values for beta1, beta2, beta3 and tau. '''

        self._curveDate = curveDate

        if len(params) != 4:
            raise ValueError("NS requires a vector/list of 4 parameters.")

        if params[3] <= 0:
            raise ValueError("Tau must be positive")

        self._beta1 = params[0]
        self._beta2 = params[1]
        self._beta3 = params[2]
        self._tau = params[3]

##########################################################################

    def zeroRate(self, dt, compoundingFreq=-1):
        ''' Calculation of zero rates with specified frequency. This 
        function can return a vector of zero rates given a vector of 
        times so must use Numpy functions. '''

        t = inputTime(dt, self)
        f = inputFrequency(compoundingFreq)
        t = t + gSmall  # To avoid overflows when t=0.0

        theta = t / self._tau
        e = np.exp(-theta)
        zeroRate = self._beta1 + self._beta2 * (1.0 - e) / theta
        zeroRate += self._beta3 * ((1.0 - e) / theta - e)
        df = np.exp(-zeroRate * t)

        if f == 0:  # Simple interest
            r = (1.0/df-1.0)/t
        if f == -1:  # Continuous
            r = -np.log(df)/t
        else:
            r = (df**(-1.0/t)-1.0) * f

        return r

##########################################################################

    def fwd(self, dt):
        ''' Calculation of forward rates. This function can return a vector
        of instantaneous forward rates given a vector of times. '''
        t = inputTime(dt, self)
        theta = t / self._tau
        e = np.exp(-theta)
        fwdRate = self.beta1 + self.beta2 * e + self.beta3 * theta * e
        return fwdRate

##########################################################################

    def df(self, dt):
        ''' Discount factor for Nelson-Siegel curve parametrisation. '''
        t = inputTime(dt, self)
        r = self.zero(t)
        return exp(-r * t)

#############################################################################


class FinNelsonSiegelSvenssonCurve():
    ''' Implementation of Nelson-Siegel-Svensson parametrisation of the
    zero rate curve '''

    def __init__(self, beta1, beta2, beta3, beta4, tau1, tau2):

        if tau1 <= 0:
            raise ValueError("Tau1 must be positive")

        if tau2 <= 0:
            raise ValueError("Tau2 must be positive")

        self._beta1 = beta1
        self._beta2 = beta2
        self._beta3 = beta3
        self._beta4 = beta4
        self._tau1 = tau1
        self._tau2 = tau2

##########################################################################

    def zero(self, t):
        ''' Calculation of zero rates. This function can return a vector
        of zero rates given a vector of times. '''

        if np.any(t < 0.0):
            raise ValueError("All times must be positive")

        t = t + gSmall  # To avoid overflows when t=0.0
        theta1 = t / self._tau1
        theta2 = t / self._tau2
        expTerm1 = np.exp(-theta1)
        expTerm2 = np.exp(-theta2)
        zeroRate = self._beta1
        zeroRate += self._beta2 * (1.0 - expTerm1) / theta1
        zeroRate += self._beta3 * ((1.0 - expTerm1) / theta1 - expTerm1)
        zeroRate += self._beta4 * ((1.0 - expTerm2) / theta2 - expTerm2)
        return zeroRate

##########################################################################

    def fwd(self, t):
        ''' Calculation of forward rates. This function uses Numpy so can return
        a vector of forward rates given a Numpy array vector of times. '''

        theta1 = t / self._tau1
        theta2 = t / self._tau2
        expTerm1 = np.exp(-theta1)
        expTerm2 = np.exp(-theta2)
        fwdRate = self._beta1
        fwdRate += self._beta2 * expTerm1
        fwdRate += self._beta3 * theta1 * expTerm1
        fwdRate += self._beta4 * theta2 * expTerm2
        return fwdRate

##########################################################################

    def df(self, t):
        ''' Discount factor for Nelson-Siegel-Svensson curve 
        parametrisation. '''
        r = self.zero(t)
        return exp(-r * t)

#############################################################################
