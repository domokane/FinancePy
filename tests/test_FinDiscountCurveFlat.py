###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat


results = [
    [0.9958, 0.9837, 0.9714, 0.9592, 0.9472, 0.9356,
     0.9239, 0.9124, 0.9010,  0.8901, 0.8789, 0.8679,
     0.8571, 0.8467, 0.8361, 0.8256, 0.8153, 0.8054,
     0.7953, 0.7853],
    [0.9959, 0.9841, 0.9721, 0.9602, 0.9485, 0.9371,
     0.9257, 0.9144, 0.9033, 0.8926, 0.8817, 0.8709,
     0.8603, 0.8501, 0.8397, 0.8294, 0.8193, 0.8096,
     0.7997, 0.7899],
    [0.9958, 0.9839, 0.9717, 0.9597, 0.9478, 0.9364,
     0.9248, 0.9134, 0.9022, 0.8914, 0.8803, 0.8694,
     0.8587, 0.8484, 0.8379, 0.8275, 0.8173, 0.8075,
     0.7975, 0.7877],
    [0.9958, 0.9838, 0.9716, 0.9595, 0.9475, 0.9360,
     0.9244, 0.9129, 0.9016, 0.8907, 0.8796, 0.8687,
     0.8579, 0.8475, 0.8370, 0.8266, 0.8163, 0.8065,
     0.7964, 0.7865],
    [0.9958, 0.9837, 0.9714, 0.9593, 0.9473, 0.9358,
     0.9241, 0.9126, 0.9012, 0.8903, 0.8792, 0.8682,
     0.8573, 0.8470, 0.8364, 0.8259, 0.8156, 0.8057,
     0.7957, 0.7857]
]


def test_FinFlatCurve():
    curve_date = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_date.add_months(months)
    compounding = FrequencyTypes.CONTINUOUS

    flat_curve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert [round(x, 4) for x in dfs] == results[0]

    compounding = FrequencyTypes.ANNUAL
    flat_curve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert [round(x, 4) for x in dfs] == results[1]

    compounding = FrequencyTypes.SEMI_ANNUAL
    flat_curve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert [round(x, 4) for x in dfs] == results[2]

    compounding = FrequencyTypes.QUARTERLY
    flat_curve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert [round(x, 4) for x in dfs] == results[3]

    compounding = FrequencyTypes.MONTHLY
    flat_curve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flat_curve.df(dates)
    assert [round(x, 4) for x in dfs] == results[4]
