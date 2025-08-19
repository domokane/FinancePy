###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.interpolator import Interpolator, InterpTypes
import numpy as np
import math


xValues = np.array([0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0])
A = -0.1
B = 0.002

yValues = []
for x in xValues:
    y = math.exp(A * x + B * x * x)
    yValues.append(y)

yValues = np.array(yValues)

xInterpolateValues = np.linspace(0.0, 10.0, 20)


def test_FinInterpolate_Runs():
    for interp_type in InterpTypes:

        yInterpValues = []

        interpolator = Interpolator(interp_type)
        interpolator.fit(xValues, yValues)

        for x in xInterpolateValues:
            y_int = interpolator.interpolate(x)
            yInterpValues.append(y_int)


def test_FinInterpolate_Recovers_Inputs():
    for interp_type in InterpTypes:

        yInterpValues = []

        interpolator = Interpolator(interp_type)
        interpolator.fit(xValues, yValues)

        for x in xValues:
            y_int = interpolator.interpolate(x)
            yInterpValues.append(y_int)

        yInterpValues = np.array(yInterpValues)
        assert np.linalg.norm(yValues - yInterpValues) / np.linalg.norm(yValues) <= 1e-6


def test_FLAT_FWD_RATES():
    interp_type = InterpTypes.FLAT_FWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 0
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.0
    assert round(y_int, 4) == 1.0

    index = 5
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.6316
    assert round(y_int, 4) == 0.7797

    index = 10
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.2632
    assert round(y_int, 4) == 0.6260


def test_LINEAR_FWD_RATES():
    interp_type = InterpTypes.LINEAR_FWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 3
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.5789
    assert round(y_int, 4) == 0.8581

    index = 15
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 7.8947
    assert round(y_int, 4) == 0.5119

    index = 19
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 10.0
    assert round(y_int, 4) == 0.4493


def test_LINEAR_ZERO_RATES():
    interp_type = InterpTypes.LINEAR_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 8
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.2105
    assert round(y_int, 4) == 0.6800

    index = 13
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.8421
    assert round(y_int, 4) == 0.5540

    index = 18
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 9.4737
    assert round(y_int, 4) == 0.4640


def test_FINCUBIC_ZERO_RATES():
    interp_type = InterpTypes.FINCUBIC_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 1
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.5263
    assert round(y_int, 4) == 0.9493

    index = 6
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 3.1579
    assert round(y_int, 4) == 0.7439

    index = 11
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.7895
    assert round(y_int, 4) == 0.6007


def test_NATCUBIC_LOG_DISCOUNT():
    interp_type = InterpTypes.NATCUBIC_LOG_DISCOUNT

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 4
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.1053
    assert round(y_int, 4) == 0.8174

    index = 9
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.7368
    assert round(y_int, 4) == 0.6512

    index = 14
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 7.3684
    assert round(y_int, 4) == 0.5355


def test_NATCUBIC_ZERO_RATES():
    interp_type = InterpTypes.NATCUBIC_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 2
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.0526
    assert round(y_int, 4) == 0.9021

    index = 7
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 3.6842
    assert round(y_int, 4) == 0.7109

    index = 12
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.3158
    assert round(y_int, 4) == 0.5759


def test_PCHIP_ZERO_RATES():
    interp_type = InterpTypes.PCHIP_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 0
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.0
    assert round(y_int, 4) == 1.0

    index = 5
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.6316
    assert round(y_int, 4) == 0.7793

    index = 10
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.2632
    assert round(y_int, 4) == 0.6244


def test_PCHIP_LOG_DISCOUNT():
    interp_type = InterpTypes.PCHIP_LOG_DISCOUNT

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 3
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.5789
    assert round(y_int, 4) == 0.8582

    index = 8
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.2105
    assert round(y_int, 4) == 0.6796

    index = 13
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.8421
    assert round(y_int, 4) == 0.5551


def test_LINEAR_ONFWD_RATES_empty_fit():
    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([], [])
    assert round(interpolator.interpolate(1.0), 4) == 1.0


def test_LINEAR_ONFWD_RATES_single_value_at_origin():
    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.0], [1.0])
    assert round(interpolator.interpolate(1.0), 4) == 1.0


def test_LINEAR_ONFWD_RATES_single_value_not_at_origin():
    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.1], [0.9])
    assert round(interpolator.interpolate(0.0), 4) == 1.0
    assert round(interpolator.interpolate(0.05), 4) == 0.9487
    assert round(interpolator.interpolate(0.1), 4) == 0.9
    assert round(interpolator.interpolate(1.0), 4) == 0.3487


def test_LINEAR_ONFWD_RATES_two_values_including_origin():
    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.0, 0.1], [1.0, 0.9])
    assert round(interpolator.interpolate(0.0), 4) == 1.0
    assert round(interpolator.interpolate(0.05), 4) == 0.9487
    assert round(interpolator.interpolate(0.1), 4) == 0.9
    assert round(interpolator.interpolate(1.0), 4) == 0.3487


def test_LINEAR_ONFWD_RATES():
    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(xValues, yValues)

    index = 3
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.8583
    print(y_int)

    index = 8
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.6802

    index = 13
    x = xInterpolateValues[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.5537


if __name__ == "__main__":
    # test_LINEAR_ONFWD_RATES_empty_fit()
    # test_LINEAR_ONFWD_RATES_single_value_at_origin()
    # test_LINEAR_ONFWD_RATES_single_value_not_at_origin()
    test_LINEAR_ONFWD_RATES_two_values_including_origin()
    # test_LINEAR_ONFWD_RATES()
