# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import numpy as np

import add_fp_to_path

from financepy.utils.math import normcdf_integrate
from financepy.utils.math import normcdf
from financepy.utils.math import normcdf_slow
from financepy.utils.math import norminvcdf
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_math():

    x_values = np.linspace(-6.0, 6.0, 13)

    start = time.time()

    test_cases.header("FUNCTION", "X", "Y")
    for x in x_values:
        y = normcdf(x)
        test_cases.print("NORMCDF1", x, y)

    end = time.time()
    duration = end - start
    test_cases.header("LABEL", "TIME")
    test_cases.print("Fast normcdf(x) takes ", duration)

    test_cases.header("FUNCTION", "X", "Y")

    start = time.time()
    for x in x_values:
        y = normcdf_slow(x)
        test_cases.print("NORMCDF2", x, y)

    end = time.time()
    duration = end - start
    test_cases.header("LABEL", "TIME")
    test_cases.print("Slow normcdf(x) takes ", duration)

    test_cases.header("FUNCTION", "X", "Y")

    start = time.time()
    for x in x_values:
        y = normcdf_integrate(x)
        test_cases.print("NORMCDF INTEGRATE", x, y)

    end = time.time()
    duration = end - start

    test_cases.header("LABEL", "TIME")
    test_cases.print("Trapezium normcdf(x) takes ", duration)

    x_values = np.linspace(-6.0, 6.0, 20)

    test_cases.header("X", "Y1", "Y2", "Y3", "DIFF1", "DIFF2")

    for x in x_values:
        y1 = normcdf(x)
        y2 = normcdf_slow(x)
        y3 = normcdf_integrate(x)
        diff1 = y3 - y1
        diff2 = y3 - y2
        test_cases.print(x, y1, y2, y3, diff1, diff2)

    x_values = np.linspace(-6.0, 6.0, 20)

    test_cases.header("X", "Y1", "Y2", "INV_Y1", "INV_Y2", "DIFF1", "DIFF2")

    for x_in in x_values:
        y1 = normcdf(x_in)
        y2 = normcdf_slow(x_in)
        x_out1 = norminvcdf(y1)
        x_out2 = norminvcdf(y2)
        diff1 = x_out1 - x_in
        diff2 = x_out2 - x_in
        test_cases.print(x, y1, y2, x_out1, x_out2, diff1, diff2)


########################################################################################

test_fin_math()
test_cases.compare_test_cases()
