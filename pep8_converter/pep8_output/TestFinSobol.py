from numba import jit
import time
from financepy.models.sobol import get_uniform_sobol, get_gaussian_sobol
from FinTestCases import FinTestCases, global_test_case_mode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test__fin_sobol():

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

    test_cases.header("LABEL", "TIME")
    start = time.time()
    for _ in range(num_repeats):
        get_uniform_sobol(1000, num_dimensions)
    end = time.time()
    test_cases.print("Average time taken", (end - start) / num_repeats)

    start = time.time()
    for _ in range(num_repeats):
        get_gaussian_sobol(1000, num_dimensions)
    end = time.time()
    test_cases.print("Average time taken", (end - start) / num_repeats)




@jit(cache=True, nopython=True)

################################################################################


def test__fin_sobol_cache():

    return get_uniform_sobol(2, 2)




################################################################################

test__fin_sobol()
test__fin_sobol_cache()
test_cases.compare_test_cases()
