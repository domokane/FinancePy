# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import math
import numpy as np
from financepy.market.curves.interpolator import Interpolator, InterpTypes


x_values = np.array([0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0])
a = -0.1
b = 0.002

y_values = []
for x in x_values:
    y = math.exp(a * x + b * x * x)
    y_values.append(y)

y_values = np.array(y_values)

x_interpolate_values = np.linspace(0.0, 10.0, 20)

########################################################################################


def test_fin_interpolate__runs():

    for interp_type in InterpTypes:

        y_interp_values = []

        interpolator = Interpolator(interp_type)
        interpolator.fit(x_values, y_values)

        for x in x_interpolate_values:
            y_int = interpolator.interpolate(x)
            y_interp_values.append(y_int)


########################################################################################


def test_fin_interpolate__recovers__inputs():

    for interp_type in InterpTypes:

        y_interp_values = []

        interpolator = Interpolator(interp_type)
        interpolator.fit(x_values, y_values)

        for x in x_values:
            y_int = interpolator.interpolate(x)
            y_interp_values.append(y_int)

        y_interp_values = np.array(y_interp_values)
        assert (
            np.linalg.norm(y_values - y_interp_values) / np.linalg.norm(y_values)
            <= 1e-6
        )


########################################################################################


def test_flat_fwd_rates():

    interp_type = InterpTypes.FLAT_FWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 0
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.0
    assert round(y_int, 4) == 1.0

    index = 5
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.6316
    assert round(y_int, 4) == 0.7797

    index = 10
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.2632
    assert round(y_int, 4) == 0.6260


########################################################################################


def test_linear_fwd_rates():

    interp_type = InterpTypes.LINEAR_FWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 3
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.5789
    assert round(y_int, 4) == 0.8581

    index = 15
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 7.8947
    assert round(y_int, 4) == 0.5119

    index = 19
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 10.0
    assert round(y_int, 4) == 0.4493


########################################################################################


def test_linear_zero_rates():

    interp_type = InterpTypes.LINEAR_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 8
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.2105
    assert round(y_int, 4) == 0.6800

    index = 13
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.8421
    assert round(y_int, 4) == 0.5540

    index = 18
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 9.4737
    assert round(y_int, 4) == 0.4640


########################################################################################


def test_fincubic_zero_rates():

    interp_type = InterpTypes.FINCUBIC_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 1
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.5263
    assert round(y_int, 4) == 0.9493

    index = 6
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 3.1579
    assert round(y_int, 4) == 0.7439

    index = 11
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.7895
    assert round(y_int, 4) == 0.6007


########################################################################################


def test_natcubic_log_discount():

    interp_type = InterpTypes.NATCUBIC_LOG_DISCOUNT

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 4
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.1053
    assert round(y_int, 4) == 0.8174

    index = 9
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.7368
    assert round(y_int, 4) == 0.6512

    index = 14
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 7.3684
    assert round(y_int, 4) == 0.5355


########################################################################################


def test_natcubic_zero_rates():

    interp_type = InterpTypes.NATCUBIC_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 2
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.0526
    assert round(y_int, 4) == 0.9021

    index = 7
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 3.6842
    assert round(y_int, 4) == 0.7109

    index = 12
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.3158
    assert round(y_int, 4) == 0.5759


########################################################################################


def test_pchip_zero_rates():

    interp_type = InterpTypes.PCHIP_ZERO_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 0
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 0.0
    assert round(y_int, 4) == 1.0

    index = 5
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 2.6316
    assert round(y_int, 4) == 0.7793

    index = 10
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 5.2632
    assert round(y_int, 4) == 0.6244


########################################################################################


def test_pchip_log_discount():

    interp_type = InterpTypes.PCHIP_LOG_DISCOUNT

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 3
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 1.5789
    assert round(y_int, 4) == 0.8582

    index = 8
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 4.2105
    assert round(y_int, 4) == 0.6796

    index = 13
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(x, 4) == 6.8421
    assert round(y_int, 4) == 0.5551


########################################################################################


def test_linear_onfwd_rates_empty_fit():

    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([], [])
    assert round(interpolator.interpolate(1.0), 4) == 1.0


########################################################################################


def test_linear_onfwd_rates_single_value_at_origin():

    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.0], [1.0])
    assert round(interpolator.interpolate(1.0), 4) == 1.0


########################################################################################


def test_linear_onfwd_rates_single_value_not_at_origin():

    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.1], [0.9])
    assert round(interpolator.interpolate(0.0), 4) == 1.0
    assert round(interpolator.interpolate(0.05), 4) == 0.9487
    assert round(interpolator.interpolate(0.1), 4) == 0.9
    assert round(interpolator.interpolate(1.0), 4) == 0.3487


########################################################################################


def test_linear_onfwd_rates_two_values_including_origin():

    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit([0.0, 0.1], [1.0, 0.9])
    assert round(interpolator.interpolate(0.0), 4) == 1.0
    assert round(interpolator.interpolate(0.05), 4) == 0.9487
    assert round(interpolator.interpolate(0.1), 4) == 0.9
    assert round(interpolator.interpolate(1.0), 4) == 0.3487


########################################################################################


def test_linear_onfwd_rates():

    interp_type = InterpTypes.LINEAR_ONFWD_RATES

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    index = 3
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.8583
    print(y_int)

    index = 8
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.6802

    index = 13
    x = x_interpolate_values[index]
    y_int = interpolator.interpolate(x)
    assert round(y_int, 4) == 0.5537


########################################################################################

########################################################################################

if __name__ == "__main__":
    # test_linear_onfwd_rates_empty_fit()
    # test_linear_onfwd_rates_single_value_at_origin()
    # test_linear_onfwd_rates_single_value_not_at_origin()
    test_linear_onfwd_rates_two_values_including_origin()
    # test_linear_onfwd_rates()
