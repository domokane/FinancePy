###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
import numpy as np


tValues = np.array([0.0, 3.0, 5.0, 10.0])
rValues = np.array([0.04, 0.07, 0.08, 0.09])
df_values = np.exp(-tValues*rValues)
tInterpValues = np.linspace(0.0, 12.0, 49)

curve_date = Date(1, 1, 2019)

tDates = curve_date.add_years(tValues)
tInterpDates = curve_date.add_years(tInterpValues)


def test_FinInterpolatedForwards():
    interp_type = InterpTypes.FLAT_FWD_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.949, 0.932, 0.916, 0.901, 0.885, 0.869, 0.855,
         0.840, 0.825, 0.811, 0.792, 0.773, 0.755, 0.737, 0.720, 0.703, 0.687,
         0.670, 0.654, 0.638, 0.622, 0.607, 0.592, 0.577, 0.563, 0.549, 0.536,
         0.523, 0.510, 0.497, 0.485, 0.473, 0.461, 0.450, 0.439, 0.428, 0.417,
         0.407, 0.397, 0.387, 0.378, 0.368, 0.359, 0.350, 0.342, 0.333]

    interp_type = InterpTypes.LINEAR_FWD_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.949, 0.932, 0.916, 0.901, 0.885, 0.869, 0.855,
         0.840, 0.825, 0.811, 0.796, 0.781, 0.764, 0.747, 0.729, 0.710, 0.691,
         0.671, 0.655, 0.639, 0.624, 0.609, 0.595, 0.580, 0.566, 0.552, 0.539,
         0.526, 0.513, 0.500, 0.488, 0.475, 0.463, 0.451, 0.440, 0.429, 0.418,
         0.407, 0.397, 0.387, 0.378, 0.368, 0.359, 0.350, 0.342, 0.333]

    interp_type = InterpTypes.LINEAR_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.949, 0.932, 0.916, 0.901, 0.885, 0.869, 0.855,
         0.840, 0.825, 0.811, 0.794, 0.776, 0.759, 0.741, 0.724, 0.706, 0.688,
         0.671, 0.656, 0.641, 0.626, 0.612, 0.598, 0.584, 0.570, 0.556, 0.542,
         0.529, 0.516, 0.503, 0.490, 0.478, 0.465, 0.453, 0.441, 0.430, 0.418,
         0.407, 0.398, 0.389, 0.380, 0.372, 0.364, 0.356, 0.348, 0.340]

    interp_type = InterpTypes.FINCUBIC_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.950, 0.934, 0.918, 0.903, 0.888, 0.873, 0.858,
         0.843, 0.827, 0.811, 0.795, 0.778, 0.760, 0.742, 0.725, 0.707, 0.688,
         0.671, 0.653, 0.636, 0.620, 0.603, 0.588, 0.573, 0.558, 0.544, 0.530,
         0.517, 0.504, 0.491, 0.480, 0.468, 0.457, 0.446, 0.436, 0.426, 0.416,
         0.407, 0.398, 0.389, 0.381, 0.372, 0.365, 0.357, 0.349, 0.342]

    interp_type = InterpTypes.NATCUBIC_LOG_DISCOUNT
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.985, 0.969, 0.954, 0.939, 0.924, 0.908, 0.893, 0.877, 0.861,
         0.845, 0.828, 0.811, 0.794, 0.776, 0.758, 0.740, 0.723, 0.705, 0.688,
         0.671, 0.654, 0.638, 0.622, 0.607, 0.592, 0.577, 0.563, 0.549, 0.536,
         0.523, 0.510, 0.497, 0.485, 0.473, 0.461, 0.450, 0.439, 0.428, 0.417,
         0.407, 0.397, 0.387, 0.378, 0.368, 0.359, 0.350, 0.342, 0.333]

    interp_type = InterpTypes.NATCUBIC_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.950, 0.934, 0.918, 0.903, 0.888, 0.873, 0.858,
         0.843, 0.827, 0.811, 0.795, 0.778, 0.760, 0.742, 0.724, 0.707, 0.688,
         0.671, 0.653, 0.636, 0.620, 0.604, 0.589, 0.574, 0.559, 0.545, 0.531,
         0.518, 0.505, 0.493, 0.481, 0.470, 0.458, 0.447, 0.437, 0.427, 0.417,
         0.407, 0.397, 0.388, 0.379, 0.370, 0.361, 0.352, 0.343, 0.334]

    interp_type = InterpTypes.PCHIP_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.983, 0.966, 0.949, 0.932, 0.916, 0.901, 0.885, 0.869, 0.855,
         0.840, 0.825, 0.811, 0.796, 0.780, 0.762, 0.743, 0.725, 0.706, 0.687,
         0.670, 0.655, 0.639, 0.623, 0.608, 0.593, 0.578, 0.564, 0.549, 0.535,
         0.522, 0.508, 0.495, 0.483, 0.471, 0.459, 0.447, 0.437, 0.426, 0.416,
         0.407, 0.398, 0.390, 0.382, 0.374, 0.368, 0.362, 0.356, 0.351]

    interp_type = InterpTypes.PCHIP_LOG_DISCOUNT
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 3) for x in dfInterpValues] == \
        [1.0, 0.986, 0.972, 0.957, 0.941, 0.925, 0.910, 0.893, 0.877, 0.861,
         0.844, 0.827, 0.811, 0.794, 0.777, 0.758, 0.740, 0.722, 0.705, 0.687,
         0.670, 0.654, 0.639, 0.623, 0.608, 0.594, 0.579, 0.565, 0.551, 0.538,
         0.525, 0.512, 0.499, 0.487, 0.475, 0.463, 0.451, 0.440, 0.429, 0.418,
         0.407, 0.397, 0.387, 0.376, 0.367, 0.357, 0.348, 0.339, 0.330]
