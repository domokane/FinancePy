# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.FinPolynomialCurve import FinPolynomialCurve
from financepy.finutils.FinDate import FinDate
import numpy as np
import sys
sys.path.append("..//..")

##########################################################################
# TODO
# Inherit from FinCurve and add df method
# Put in a convention for the rate
# Use Frequency object
##########################################################################

showPlots = True

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinPolynomialCurve():

    times = np.linspace(0.00, 10.0, 20)
    curveDate = FinDate(2019, 2, 2)
    coeffs = [0.0004, -0.0001, 0.00000010]
    curve1 = FinPolynomialCurve(curveDate, coeffs)
    zeros = curve1.zeroRate(times)
    fwds = curve1.fwd(times)

    if showPlots:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        plt.plot(times, zeros, label="Zeros")
        plt.plot(times, fwds, label="Forwards")
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate')
        plt.legend(loc='best')


test_FinPolynomialCurve()
testCases.compareTestCases()
