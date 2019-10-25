# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np

################################################################################
# TODO
# Generalise to polynomial and not just cubic or rename as cubic
# Implement a spline
################################################################################

from ...market.curves.FinCurve import FinCurve

################################################################################

class FinPolynomialCurve(FinCurve):
    ''' Curve with zero rate parametrised as a cubic polynomial. '''

    def __init__(self,a=0,b=0,c=0,d=0):
        ''' Create cubic curve from coefficients. '''

        self.c0 = a
        self.c1 = b
        self.c2 = c
        self.c3 = d

    def zero(self,t):
        zeroRate = self.c0 + self.c1*t + self.c2*t*t + self.c3*t*t*t
        return zeroRate

    def fwd(self,t):
        fwdRate = self.c1 + 2*self.c2*t + 3*self.c3*t*t
        return fwdRate

    def df(self,t):
        r = self.zero(t)
        return np.exp(-r*t)

################################################################################

