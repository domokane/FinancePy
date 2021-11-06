###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.math import normcdf_slow, N, normcdf_integrate
from financepy.utils.math import accrued_interpolator
import numpy as np


x = -3.0


def test_normcdf1():
    result = N(x)
    assert round(result, 5) == 0.00135


def test_normcdf2():
    result = normcdf_slow(x)
    assert round(result, 5) == 0.00135


def test_normcdf_integrate():
    result = normcdf_integrate(x)
    assert round(result, 5) == 0.00135


def test_accrued_interpolator():
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

    assert [round(x * 1000, 6) for x in values] == \
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 8.691022, 17.420288, 26.149555, 0.266336,
         8.950803, 17.63527, 26.319737, 35.004204, 0.0]
