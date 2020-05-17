# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:10:12 2019

@author: Dominic
"""

from math import log, exp

##########################################################################
# dr = theta(t) dt + sigma * dW
##########################################################################

# TO DO - DECIDE WHETHER TO OO MODEL


class FinModelRatesHL():

    def __init__(self, discountCurve, sigma):
        self._discountCurve = discountCurve
        self._sigma = sigma

##########################################################################

    def P(self,
          r1,  # short rate at time t1
          t1,  # foward start time t1
          t2):  # forward maturity t2

        tau = t2 - t1
        dt = 1e-10
        pt1 = self._discountCurve.df(t1)
        pt1p = self._discountCurve.df(t1 + dt)
        pt2 = self._discountCurve.df(t2)
        f0t1 = -log(pt1p / pt1) / dt

        sigma2 = self._sigma**2
        lnA = log(pt2 / pt1) + tau * f0t1 - 0.5 * sigma2 * (tau**2)
        df = exp(lnA - r1 * tau)
        return df

##########################################################################
