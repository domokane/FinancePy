# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import matplotlib.pyplot as plt
import math
import numpy as np
from financepy.market.curves.interpolator import Interpolator, InterpTypes
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)

plot_graphs = False



########################################################################################


def test_FinInterpolate():

    import time

    x_values = np.array([0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0])
    a = -0.1
    b = 0.002

    y_values = []
    for x in x_values:
        y = math.exp(a * x + b * x * x)
        y_values.append(y)

    y_values = np.array(y_values)

    x_interpolate_values = np.linspace(0.0, 10.0, 20)

    test_cases.header("METHOD", "X", "Y_INTERPOLATED")

    for interp_type in InterpTypes:

        y_interp_values = []
        start = time.time()

        interpolator = Interpolator(interp_type)
        interpolator.fit(x_values, y_values)

        for x in x_interpolate_values:
            y_int = interpolator.interpolate(x)
            test_cases.print(interp_type, x, y_int)
            y_interp_values.append(y_int)

        end = time.time()

        if plot_graphs:
            plt.figure(figsize=(12, 10))
            plt.plot(x_values, y_values, color="r", marker="o")
            plt.plot(
                x_interpolate_values,
                y_interp_values,
                color="b",
                label=str(interp_type),
            )
            plt.legend()

    xp = np.array([0.2, 0.4, 0.45, 0.6, 0.82, 0.93, 0.99])
    yp = np.array([0.4, 0.9, 0.32, 0.2, 0.22, 0.10, 0.28])
    n = 10000

    test_cases.header("LABEL", "TIME")
    interpolator = Interpolator(interp_type)
    interpolator.fit(xp, yp)

    start = time.time()
    for i in range(0, n):
        interpolator.interpolate(0.8)
    end = time.time()
    test_cases.print("10000 Interpolations", end - start)


########################################################################################


test_FinInterpolate()
test_cases.compare_test_cases()

########################################################################################

