# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:10:12 2019

@author: Dominic
"""

from math import log, exp
#from numba import njit, float64, int64

##########################################################################
# dr = theta(t) dt + sigma * dW
##########################################################################

# TO DO - DECIDE WHETHER TO OO MODEL


class FinHullWhiteTree():

    def __init__(self,
                 discountCurve,
                 a,
                 sigma):
        self._discountCurve = discountCurve
        self._a = a
        self._sigma = sigma


class FinHullWhiteRateModel():

    def __init__(self, discountCurve,
                 a,
                 sigma):
        ''' Create Hull White Model from time zero discount curve. '''
        self._discountCurve = discountCurve
        self._a = a
        self._sigma = sigma

##########################################################################

    def P(self,
          r1,  # short rate at time t1
          t1,  # forward start time t1
          t2):  # forward maturity t2
        ''' From discount curve at time 0 get price of 1 dollar at time t2 as of t1 '''

        tau = t2 - t1
        b = (1.0 - exp(-self._a * tau)) / self._a

        dt = 1e-10
        pt1 = self._discountCurve.df(t1)
        pt1p = self._discountCurve.df(t1 + dt)
        pt2 = self._discountCurve.df(t1)
        f0t1 = -log(pt1p / pt1) / dt

        sigma2 = self._sigma**2
        q1 = exp(-self._a * t1)
        q2 = exp(-self._a * t2)

        lnA = log(pt2 / pt1) + b * f0t1
        lnA = lnA - sigma2 * (q2 - q1) * (q2 - q1) * \
            (exp(2.0 * self._a * t1) - 1.0) / 4 / (self._a**3)
        df = exp(lnA - r1 * b)
        return df

##########################################################################
