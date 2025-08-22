# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np
from financepy.utils.math import accrued_interpolator
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)

plt_graph = False

########################################################################################


def test_fin_math_accd_interpolator():

    cpn_times = [0.0, 4.000087613144162, 4.495649459810208, 5.002162949496498]
    cpn_flows = [0.0, 0.0, 0.03461111111111111, 0.035194444444444445]

    tree_times = [
        0.0,
        0.12498563,
        0.24997125,
        0.37495688,
        0.4999425,
        0.62492813,
        0.74991376,
        0.87489938,
        0.99988501,
        1.12487063,
        1.24985626,
        1.37484189,
        1.49982751,
        1.62481314,
        1.74979876,
        1.87478439,
        1.99977002,
        2.12475564,
        2.24974127,
        2.37472689,
        2.49971252,
        2.62469815,
        2.74968377,
        2.8746694,
        2.99965502,
        3.12464065,
        3.24962628,
        3.3746119,
        3.49959753,
        3.62458315,
        3.74956878,
        3.87455441,
        3.99954003,
        4.12452566,
        4.24951128,
        4.37449691,
        4.49948253,
        4.62446816,
        4.74945379,
        4.87443941,
        4.99942504,
        5.12441066,
    ]

    cpn_times = np.array(cpn_times)
    cpn_flows = np.array(cpn_flows)

    values = []

    for t in tree_times:
        v = accrued_interpolator(t, cpn_times, cpn_flows)
        values.append(v)

    test_cases.header("VALUE")
    test_cases.print(values)

    if plt_graph:
        plt.plot(tree_times, values)


########################################################################################

test_fin_math_accd_interpolator()
test_cases.compare_test_cases()
