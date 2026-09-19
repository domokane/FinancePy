# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import time
import numpy as np

import add_fp_to_path

from financepy.utils.math import normcdf_integrate
from financepy.utils.math import normcdf
from financepy.utils.math import normcdf_slow
from financepy.utils.math import norminvcdf


########################################################################################


def test_fin_math():

    x_values = np.linspace(-6.0, 6.0, 13)

    start = time.time()

    print("FUNCTION", "X", "Y")
    for x in x_values:
        y = normcdf(x)
        print("NORMCDF1", x, y)

    end = time.time()
    duration = end - start
    print("LABEL", "TIME")
    print("Fast normcdf(x) takes ", duration)

    print("FUNCTION", "X", "Y")

    start = time.time()
    for x in x_values:
        y = normcdf_slow(x)
        print("NORMCDF2", x, y)

    end = time.time()
    duration = end - start
    print("LABEL", "TIME")
    print("Slow normcdf(x) takes ", duration)

    print("FUNCTION", "X", "Y")

    start = time.time()
    for x in x_values:
        y = normcdf_integrate(x)
        print("NORMCDF INTEGRATE", x, y)

    end = time.time()
    duration = end - start

    print("LABEL", "TIME")
    print("Trapezium normcdf(x) takes ", duration)

    x_values = np.linspace(-6.0, 6.0, 20)

    print("X", "Y1", "Y2", "Y3", "DIFF1", "DIFF2")

    for x in x_values:
        y1 = normcdf(x)
        y2 = normcdf_slow(x)
        y3 = normcdf_integrate(x)
        diff1 = y3 - y1
        diff2 = y3 - y2
        print(x, y1, y2, y3, diff1, diff2)

    x_values = np.linspace(-6.0, 6.0, 20)

    print("X", "Y1", "Y2", "INV_Y1", "INV_Y2", "DIFF1", "DIFF2")

    for x_in in x_values:
        y1 = normcdf(x_in)
        y2 = normcdf_slow(x_in)
        x_out1 = norminvcdf(y1)
        x_out2 = norminvcdf(y2)
        diff1 = x_out1 - x_in
        diff2 = x_out2 - x_in
        print(x, y1, y2, x_out1, x_out2, diff1, diff2)


########################################################################################

test_fin_math()
