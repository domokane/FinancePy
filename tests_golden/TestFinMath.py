###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from financepy.utils.math import normcdf_integrate
from financepy.utils.math import N
from financepy.utils.math import normcdf_slow
from financepy.utils.math import norminvcdf
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinMath():

    xValues = np.linspace(-6.0, 6.0, 13)

    start = time.time()

    testCases.header("FUNCTION", "X", "Y")
    for x in xValues:
        y = N(x)
        testCases.print("NORMCDF1", x, y)

    end = time.time()
    duration = end - start
    testCases.header("LABEL", "TIME")
    testCases.print("Fast N(x) takes ", duration)

    ##########################################################################

    testCases.header("FUNCTION", "X", "Y")

    start = time.time()
    for x in xValues:
        y = normcdf_slow(x)
        testCases.print("NORMCDF2", x, y)

    end = time.time()
    duration = end - start
    testCases.header("LABEL", "TIME")
    testCases.print("Slow N(x) takes ", duration)

    ##########################################################################

    testCases.header("FUNCTION", "X", "Y")

    start = time.time()
    for x in xValues:
        y = normcdf_integrate(x)
        testCases.print("NORMCDF INTEGRATE", x, y)

    end = time.time()
    duration = end - start

    testCases.header("LABEL", "TIME")
    testCases.print("Trapezium N(x) takes ", duration)

    ##########################################################################

    xValues = np.linspace(-6.0, 6.0, 20)

    testCases.header("X", "Y1", "Y2", "Y3", "DIFF1", "DIFF2")

    for x in xValues:
        y1 = N(x)
        y2 = normcdf_slow(x)
        y3 = normcdf_integrate(x)
        diff1 = y3 - y1
        diff2 = y3 - y2
        testCases.print(x, y1, y2, y3, diff1, diff2)

    ##########################################################################

    xValues = np.linspace(-6.0, 6.0, 20)

    testCases.header("X", "Y1", "Y2", "INV_Y1", "INV_Y2", "DIFF1", "DIFF2")

    for x_in in xValues:
        y1 = N(x_in)
        y2 = normcdf_slow(x_in)
        x_out1 = norminvcdf(y1)
        x_out2 = norminvcdf(y2)
        diff1 = x_out1 - x_in
        diff2 = x_out2 - x_in
        testCases.print(x, y1, y2, x_out1, x_out2, diff1, diff2)


##########################################################################


test_FinMath()
testCases.compareTestCases()
