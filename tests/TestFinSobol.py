
import time
from numba import jit

import sys
sys.path.append("..")

from financepy.models.FinSobol import getUniformSobol, getGaussianSobol

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinSobol():

    numPoints = 1000
    dimensions = 3

    points = getUniformSobol(numPoints, dimensions)

    for d in range(dimensions):
        av = 0.0
        var = 0.0

        for point in points[:, d]:
            av += point
            var += point ** 2

        av /= numPoints
        var /= numPoints

        avError = abs(av - (1/2))
        varError = abs(var - (1/3))
        assert(avError < 0.002)
        assert(varError < 0.002)

    numRepeats = 100
    numDimensions = 10

    testCases.header("LABEL", "TIME")
    start = time.time()
    for _ in range(numRepeats):
        getUniformSobol(1000, numDimensions)
    end = time.time()
    testCases.print("Average time taken", (end - start) / numRepeats)

    start = time.time()
    for _ in range(numRepeats):
        getGaussianSobol(1000, numDimensions)
    end = time.time()
    testCases.print("Average time taken", (end - start) / numRepeats)


###############################################################################

@jit(cache=True, nopython=True)
def test_FinSobolCache():
    return getUniformSobol(2, 2)

###############################################################################


test_FinSobol()
test_FinSobolCache()
testCases.compareTestCases()
