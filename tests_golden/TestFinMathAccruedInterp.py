##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import matplotlib.pyplot as plt
import numpy as np
from financepy.utils.math import accrued_interpolator
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

pltGraph = False

##########################################################################


def test_FinMathAccdInterpolator():

    coupon_times = [0.0, 4.000087613144162,
                    4.495649459810208, 5.002162949496498]
    coupon_flows = [0.0, 0.0, 0.03461111111111111, 0.035194444444444445]

    tree_times = [0., 0.12498563, 0.24997125, 0.37495688, 0.4999425, 0.62492813,
                  0.74991376, 0.87489938, 0.99988501, 1.12487063, 1.24985626, 1.37484189,
                  1.49982751, 1.62481314, 1.74979876, 1.87478439, 1.99977002, 2.12475564,
                  2.24974127, 2.37472689, 2.49971252, 2.62469815, 2.74968377, 2.8746694,
                  2.99965502, 3.12464065, 3.24962628, 3.3746119, 3.49959753, 3.62458315,
                  3.74956878, 3.87455441, 3.99954003, 4.12452566, 4.24951128, 4.37449691,
                  4.49948253, 4.62446816, 4.74945379, 4.87443941, 4.99942504, 5.12441066]

    coupon_times = np.array(coupon_times)
    coupon_flows = np.array(coupon_flows)

    values = []

    for t in tree_times:
        v = accrued_interpolator(t, coupon_times, coupon_flows)
        values.append(v)

    testCases.header("VALUE")
    testCases.print(values)

    if pltGraph:
        plt.plot(tree_times, values)

##########################################################################


test_FinMathAccdInterpolator()
testCases.compareTestCases()
