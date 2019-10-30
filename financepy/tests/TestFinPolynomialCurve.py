# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import sys
sys.path.append("..//..")

################################################################################
# TODO
# Inherit from FinCurve and add df method
# Put in a convention for the rate
# Use Frequency object
################################################################################


import numpy as np
#import matplotlib.pyplot as plt

from financepy.market.curves.FinPolynomialCurve import FinPolynomialCurve
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinPolynomialCurve(): 
    
    times = np.linspace(0.0,10.0,5)

    curve1 = FinPolynomialCurve(1.,0.,0.,0.)
    factor1loading = curve1.zero(times)
    curve2 = FinPolynomialCurve(1.,1.,0.,0.)
    factor2loading = curve2.zero(times)
    curve3 = FinPolynomialCurve(1.,1.,1.,0)
    factor3loading = curve3.zero(times)
    curve4 = FinPolynomialCurve(1.,1.,1.,1.)
    factor4loading = curve4.zero(times)

    testCases.header("FACTOR LOADING")
    testCases.print(factor1loading)
    testCases.print(factor2loading)
    testCases.print(factor3loading)
    testCases.print(factor4loading)

#    plt.figure(figsize = (6,4))
#    plt.plot(times,scaleVector(factor1loading,1),label='c0=1');
#    plt.plot(times,scaleVector(factor2loading,1),label='c0=c1=1');
#    plt.plot(times,scaleVector(factor3loading,1),label='c0=c1=c2=1');
#    plt.plot(times,scaleVector(factor4loading,1),label='c0=c1=c2=c3=1');
#    plt.ylim((0,1000))
#
#    plt.title('Factor Loadings in Nelson-Siegel Model');
#    plt.xlabel('Time (years)');
#    plt.ylabel('Loading');
#    plt.legend(loc='best')


test_FinPolynomialCurve()
testCases.compareTestCases()