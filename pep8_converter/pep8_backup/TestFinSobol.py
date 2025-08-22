from numba import jit
import time
from financepy.models.sobol import get_uniform_sobol, get_gaussian_sobol
from FinTestCases import FinTestCases, global_test_case_mode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_FinSobol():

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

        avError = abs(av - (1 / 2))
        varError = abs(var - (1 / 3))
        assert avError < 0.002
        assert varError < 0.002

    numRepeats = 100
    numDimensions = 10

    test_cases.header("LABEL", "TIME")
    start = time.time()
    for _ in range(numRepeats):
        get_uniform_sobol(1000, numDimensions)
    end = time.time()
    test_cases.print("Average time taken", (end - start) / numRepeats)

    start = time.time()
    for _ in range(numRepeats):
        get_gaussian_sobol(1000, numDimensions)
    end = time.time()
    test_cases.print("Average time taken", (end - start) / numRepeats)


########################################################################################


@jit(cache=True, nopython=True)
def test_FinSobolCache():
    return get_uniform_sobol(2, 2)


########################################################################################


test_FinSobol()
test_FinSobolCache()
test_cases.compare_test_cases()
