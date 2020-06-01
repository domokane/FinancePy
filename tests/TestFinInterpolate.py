# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.market.curves.FinInterpolate import interpolate, FinInterpMethods
import numpy as np
import math
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinInterpolate():

    import time

    xValues = np.array([0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0])
    a = -0.1
    b = 0.002

    yValues = []
    for x in xValues:
        y = math.exp(a * x + b * x * x)
        yValues.append(y)

    yValues = np.array(yValues)

    xInterpolateValues = np.linspace(0.0, 10.0, 20)

    testCases.header("METHOD", "X", "Y_INTERPOLATED")

    for method in FinInterpMethods:

        yInterpValues = []
        start = time.time()

        for x in xInterpolateValues:
            y_int = interpolate(x, xValues, yValues, method.value)
            testCases.print(method, x, y_int)
            yInterpValues.append(y_int)

        end = time.time()

        import matplotlib.pyplot as plt
        plt.figure(figsize=(12, 10))
        plt.plot(xValues, yValues, color='r', marker='o')
        plt.plot(xInterpolateValues, yInterpValues, color='b', label=str(method))
        plt.legend()

    xp = np.array([0.2, 0.4, 0.45, 0.6, 0.82, 0.93, 0.99])
    yp = np.array([0.4, 0.9, 0.32, 0.2, 0.22, 0.10, 0.28])
    n = 10000

    testCases.header("LABEL", "TIME")

    start = time.time()
    for i in range(0, n):
        interpolate(0.8, xp, yp, method.value)
    end = time.time()
    testCases.print("10000 Interpolations", end - start)


test_FinInterpolate()
testCases.compareTestCases()
