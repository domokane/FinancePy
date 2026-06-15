# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat


results = [
    [
        0.9958,
        0.9837,
        0.9714,
        0.9592,
        0.9472,
        0.9356,
        0.9239,
        0.9123,
        0.9009,
        0.8900,
        0.8788,
        0.8678,
        0.8569,
        0.8466,
        0.8360,
        0.8255,
        0.8151,
        0.8053,
        0.7952,
        0.7852,
    ],
    [
        0.9959,
        0.9841,
        0.9721,
        0.9602,
        0.9484,
        0.9371,
        0.9256,
        0.9144,
        0.9033,
        0.8926,
        0.8817,
        0.8709,
        0.8603,
        0.8501,
        0.8397,
        0.8294,
        0.8193,
        0.8096,
        0.7997,
        0.7899,
    ],
    [
        0.9958,
        0.9839,
        0.9717,
        0.9597,
        0.9478,
        0.9364,
        0.9248,
        0.9134,
        0.9022,
        0.8914,
        0.8803,
        0.8694,
        0.8587,
        0.8484,
        0.8379,
        0.8275,
        0.8173,
        0.8075,
        0.7975,
        0.7877,
    ],
    [
        0.9958,
        0.9838,
        0.9716,
        0.9595,
        0.9475,
        0.9360,
        0.9244,
        0.9129,
        0.9016,
        0.8907,
        0.8796,
        0.8687,
        0.8579,
        0.8475,
        0.8370,
        0.8266,
        0.8163,
        0.8065,
        0.7964,
        0.7865,
    ],
    [
        0.9958,
        0.9837,
        0.9714,
        0.9593,
        0.9473,
        0.9358,
        0.9241,
        0.9126,
        0.9012,
        0.8903,
        0.8792,
        0.8682,
        0.8573,
        0.8470,
        0.8364,
        0.8259,
        0.8156,
        0.8057,
        0.7957,
        0.7857,
    ],
]

########################################################################################


def test_fin_flat_curve():

    curve_dt = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_dt.add_months(months)
    compounding = FrequencyTypes.CONTINUOUS

    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)

    print(dfs)
    print(results[0])
    assert np.allclose(dfs, results[0], atol=1e-3)

    compounding = FrequencyTypes.ANNUAL
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert np.allclose(dfs, results[1], atol=1e-3)

    compounding = FrequencyTypes.SEMI_ANNUAL
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert np.allclose(dfs, results[2], atol=1e-3)

    compounding = FrequencyTypes.QUARTERLY
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert np.allclose(dfs, results[3], atol=1e-3)

    compounding = FrequencyTypes.MONTHLY
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert np.allclose(dfs, results[4], atol=1e-3)

if __name__ == "__main__":
    test_fin_flat_curve()
