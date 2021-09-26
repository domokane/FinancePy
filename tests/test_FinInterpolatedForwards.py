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
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9829, 0.9659, 0.949, 0.9325, 0.9163, 0.9005, 0.8848, 
         0.8693, 0.8544, 0.8397, 0.825,  0.8106, 0.7918, 0.7733, 0.755, 
         0.7371, 0.7201, 0.7032, 0.6866, 0.6703, 0.6538, 0.6378, 0.6219, 
         0.6064, 0.5917, 0.5771, 0.5628, 0.5488, 0.5354, 0.5223, 0.5093, 
         0.4966, 0.4845, 0.4726, 0.4609, 0.4494, 0.4383, 0.4276, 0.4169, 
         0.4066, 0.3967, 0.3869, 0.3773, 0.3679, 0.359, 0.3501, 0.3414, 0.3329]

    interp_type = InterpTypes.LINEAR_FWD_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9829, 0.9659, 0.949, 0.9325, 0.9163, 0.9005, 0.8848, 
         0.8693, 0.8544, 0.8397, 0.825, 0.8106, 0.7961, 0.7805, 0.7639, 
         0.7464, 0.7286, 0.7099, 0.6904, 0.6703, 0.6546, 0.6392, 0.6238, 
         0.6088, 0.5944, 0.5801, 0.5659, 0.552, 0.5387, 0.5255, 0.5124, 
         0.4995, 0.4872, 0.4751, 0.463, 0.4512, 0.4397, 0.4285, 0.4174, 
         0.4066, 0.3967, 0.3869, 0.3773, 0.3679, 0.359, 0.3501, 0.3414, 0.3329]

    interp_type = InterpTypes.LINEAR_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9829, 0.9659, 0.949, 0.9325, 0.9163, 0.9005, 0.8848, 
         0.8693, 0.8544, 0.8397, 0.825, 0.8106, 0.7935, 0.7762, 0.7585, 
         0.7408, 0.7235, 0.7059, 0.6881, 0.6703, 0.6554, 0.6406, 0.6259, 
         0.6113, 0.5972, 0.5832, 0.5692, 0.5554, 0.5421, 0.5288, 0.5156, 
         0.5026, 0.4901, 0.4776, 0.4652, 0.453, 0.4412, 0.4295, 0.4179, 
         0.4066, 0.3977, 0.3888, 0.3801, 0.3716, 0.3635, 0.3554, 0.3474, 0.3396]

    interp_type = InterpTypes.FINCUBIC_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.983, 0.9663, 0.9499, 0.9338, 0.9183, 0.903, 0.8878, 
         0.8725, 0.8576, 0.8424, 0.8267, 0.8106, 0.7944, 0.7775, 0.7599, 
         0.7421, 0.7244, 0.7064, 0.6882, 0.6703, 0.6529, 0.636, 0.6193, 
         0.603, 0.5876, 0.5724, 0.5576, 0.5433, 0.5297, 0.5165, 0.5035, 
         0.4911, 0.4793, 0.4678, 0.4566, 0.4459, 0.4356, 0.4256, 0.4159, 
         0.4066, 0.3977, 0.389, 0.3805, 0.3722, 0.3643, 0.3566, 0.3489, 0.3414]

    interp_type = InterpTypes.NATCUBIC_LOG_DISCOUNT
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9847, 0.9694, 0.9541, 0.9387, 0.9235, 0.9082, 0.8925, 
         0.8766, 0.8608, 0.8445, 0.8278, 0.8106, 0.7934, 0.7758, 0.7579, 
         0.7399, 0.7224, 0.7049, 0.6874, 0.6703, 0.6538, 0.6377, 0.6219, 
         0.6064, 0.5916, 0.5771, 0.5627, 0.5487, 0.5354, 0.5222, 0.5092, 
         0.4966, 0.4845, 0.4726, 0.4608, 0.4494, 0.4383, 0.4276, 0.4169, 
         0.4066, 0.3967, 0.3869, 0.3773, 0.3679, 0.359, 0.3502, 0.3415, 0.333]

    interp_type = InterpTypes.NATCUBIC_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.983, 0.9663, 0.9499, 0.9338, 0.9183, 0.9031, 0.8878, 
         0.8726, 0.8576, 0.8424, 0.8267, 0.8106, 0.7943, 0.7774, 0.7599, 
         0.742, 0.7243, 0.7063, 0.6882, 0.6703, 0.653, 0.6362, 0.6196, 
         0.6035, 0.5882, 0.5733, 0.5586, 0.5444, 0.531, 0.5179, 0.5051, 
         0.4927, 0.4809, 0.4694, 0.4582, 0.4472, 0.4367, 0.4265, 0.4164, 
         0.4066, 0.3971, 0.3878, 0.3785, 0.3693, 0.3604, 0.3515, 0.3425, 0.3336]

    interp_type = InterpTypes.PCHIP_ZERO_RATES
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9829, 0.9659, 0.949, 0.9325, 0.9163, 0.9005, 0.8848, 
         0.8693, 0.8544, 0.8397, 0.825, 0.8106, 0.7959, 0.7795, 0.7617, 
         0.7431, 0.7244, 0.7056, 0.6873, 0.6703, 0.6545, 0.6388, 0.6231, 
         0.6077, 0.5928, 0.5779, 0.5632, 0.5488, 0.535, 0.5213, 0.5079, 
         0.4949, 0.4825, 0.4704, 0.4585, 0.4472, 0.4363, 0.426, 0.416, 
         0.4066, 0.3978, 0.3894, 0.3815, 0.3742, 0.3675, 0.3613, 0.3557, 0.3506]

    interp_type = InterpTypes.PCHIP_LOG_DISCOUNT
    discount_curve = DiscountCurve(
        curve_date, tDates, df_values, interp_type)
    dfInterpValues = discount_curve.df(tInterpDates)
    assert [round(x, 4) for x in dfInterpValues] == \
        [1.0, 0.9862, 0.9717, 0.9566, 0.9411, 0.9254, 0.9095, 0.8932, 
         0.8767, 0.8604, 0.8439, 0.8272, 0.8106, 0.7939, 0.7764, 0.7583, 
         0.74, 0.7222, 0.7045, 0.6871, 0.6703, 0.6543, 0.6386, 0.6231, 
         0.6079, 0.5934, 0.579, 0.5647, 0.5508, 0.5375, 0.5243, 0.5113, 
         0.4985, 0.4863, 0.4743, 0.4623, 0.4506, 0.4393, 0.4282, 0.4173, 
         0.4066, 0.3963, 0.3862, 0.3761, 0.3663, 0.3569, 0.3476, 0.3385, 0.3295]
