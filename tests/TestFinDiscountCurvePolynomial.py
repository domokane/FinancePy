###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.FinDiscountCurvePoly import FinDiscountCurvePoly
from financepy.finutils.FinDate import FinDate
import numpy as np
import sys
sys.path.append("..//..")

##########################################################################
# TODO
# Inherit from FinDiscountCurve and add df method
# Put in a convention for the rate
# Use Frequency object
##########################################################################

showPlots = True

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinDiscountCurvePolynomial():

    times = np.linspace(0.00, 10.0, 21)
    curveDate = FinDate(2019, 2, 2)
    dates = curveDate.addYears(times)
    coeffs = [0.0004, -0.0001, 0.00000010]
    curve1 = FinDiscountCurvePoly(curveDate, coeffs)
    zeros = curve1.zeroRate(dates)
    fwds = curve1.fwd(dates)

    if showPlots:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        plt.plot(times, zeros, label="Zeros")
        plt.plot(times, fwds, label="Forwards")
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate')
        plt.legend(loc='best')


test_FinDiscountCurvePolynomial()
testCases.compareTestCases()
