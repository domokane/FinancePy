
# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
from numba import jit
import time

import add_fp_to_path

from financepy.models.sobol import get_uniform_sobol, get_gaussian_sobol


########################################################################################


def test_fin_sobol():

    num_points = 1000
    dimensions = 3

    points = get_uniform_sobol(num_points, dimensions)

    for d in range(dimensions):
        av = 0.0
        var = 0.0

        for point in points[:, d]:
            av += point
            var += point**2

        av /= num_points
        var /= num_points

        av_error = abs(av - (1 / 2))
        var_error = abs(var - (1 / 3))
        assert av_error < 0.002
        assert var_error < 0.002

    num_repeats = 100
    num_dimensions = 10

    print("LABEL", "TIME")
    start = time.time()
    for _ in range(num_repeats):
        get_uniform_sobol(1000, num_dimensions)
    end = time.time()
    print("Average time taken", (end - start) / num_repeats)

    start = time.time()
    for _ in range(num_repeats):
        get_gaussian_sobol(1000, num_dimensions)
    end = time.time()
    print("Average time taken", (end - start) / num_repeats)


@jit(cache=True, nopython=True)

########################################################################################


def test_fin_sobol_cache():

    return get_uniform_sobol(2, 2)


########################################################################################

test_fin_sobol()
test_fin_sobol_cache()
