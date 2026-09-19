# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest
import numpy as np

from financepy.models.black_shifted import BlackShifted
from financepy.utils.global_types import OptionTypes


def test_black_shifted_stores_inputs():

    model = BlackShifted(
        volatility=0.20,
        shift=0.01,
        implementation=7,
    )

    assert model.volatility == 0.20
    assert model.shift == 0.01
    assert model.implementation == 7


def test_black_shifted_call_value():

    model = BlackShifted(0.20, 0.01)

    value = model.value(
        forward_rate=0.03,
        strike_rate=0.025,
        time_to_expiry=2.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(value, 0.006816152949408656, atol=1e-12)


def test_black_shifted_put_value():

    model = BlackShifted(0.20, 0.01)

    value = model.value(
        forward_rate=0.03,
        strike_rate=0.025,
        time_to_expiry=2.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_PUT,
    )

    assert np.isclose(value, 0.002066152949408657, atol=1e-12)


def test_black_shifted_reduces_to_black_when_shift_is_zero():

    shifted = BlackShifted(0.20, 0.00)
    unshifted = BlackShifted(0.20, 0.00)

    f = 0.03
    k = 0.025
    t = 2.0
    df = 0.95

    v1 = shifted.value(f, k, t, df, OptionTypes.EUROPEAN_CALL)
    v2 = unshifted.value(f, k, t, df, OptionTypes.EUROPEAN_CALL)

    assert np.isclose(v1, v2, atol=1e-14)


def test_black_shifted_handles_negative_rates_with_positive_shift():

    model = BlackShifted(0.20, 0.02)

    value = model.value(
        forward_rate=-0.005,
        strike_rate=-0.010,
        time_to_expiry=1.0,
        df=0.99,
        call_or_put=OptionTypes.EUROPEAN_CALL,
    )

    assert value > 0.0


def test_black_shifted_zero_expiry_call():

    model = BlackShifted(0.20, 0.01)

    value = model.value(
        forward_rate=0.03,
        strike_rate=0.025,
        time_to_expiry=0.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(value, 0.95 * 0.005, atol=1e-14)


def test_black_shifted_zero_expiry_put():

    model = BlackShifted(0.20, 0.01)

    value = model.value(
        forward_rate=0.02,
        strike_rate=0.025,
        time_to_expiry=0.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_PUT,
    )

    assert np.isclose(value, 0.95 * 0.005, atol=1e-14)


def test_black_shifted_zero_volatility_call():

    model = BlackShifted(0.0, 0.01)

    value = model.value(
        forward_rate=0.03,
        strike_rate=0.025,
        time_to_expiry=2.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(value, 0.95 * 0.005, atol=1e-14)


def test_black_shifted_zero_volatility_put():

    model = BlackShifted(0.0, 0.01)

    value = model.value(
        forward_rate=0.02,
        strike_rate=0.025,
        time_to_expiry=2.0,
        df=0.95,
        call_or_put=OptionTypes.EUROPEAN_PUT,
    )

    assert np.isclose(value, 0.95 * 0.005, atol=1e-14)


def test_black_shifted_rejects_negative_expiry():

    model = BlackShifted(0.20, 0.01)

    with pytest.raises(ValueError):
        model.value(
            forward_rate=0.03,
            strike_rate=0.025,
            time_to_expiry=-1.0,
            df=0.95,
            call_or_put=OptionTypes.EUROPEAN_CALL,
        )


def test_black_shifted_rejects_negative_discount_factor():

    model = BlackShifted(0.20, 0.01)

    with pytest.raises(ValueError):
        model.value(
            forward_rate=0.03,
            strike_rate=0.025,
            time_to_expiry=1.0,
            df=-0.95,
            call_or_put=OptionTypes.EUROPEAN_CALL,
        )


def test_black_shifted_rejects_invalid_shifted_forward():

    model = BlackShifted(0.20, 0.01)

    with pytest.raises(ValueError):
        model.value(
            forward_rate=-0.02,
            strike_rate=0.025,
            time_to_expiry=1.0,
            df=0.95,
            call_or_put=OptionTypes.EUROPEAN_CALL,
        )


def test_black_shifted_rejects_invalid_shifted_strike():

    model = BlackShifted(0.20, 0.01)

    with pytest.raises(ValueError):
        model.value(
            forward_rate=0.03,
            strike_rate=-0.02,
            time_to_expiry=1.0,
            df=0.95,
            call_or_put=OptionTypes.EUROPEAN_CALL,
        )


def test_black_shifted_rejects_non_european_option_type():

    model = BlackShifted(0.20, 0.01)

    with pytest.raises(ValueError):
        model.value(
            forward_rate=0.03,
            strike_rate=0.025,
            time_to_expiry=1.0,
            df=0.95,
            call_or_put=OptionTypes.AMERICAN_CALL,
        )


def test_black_shifted_repr():

    model = BlackShifted(0.20, 0.01, implementation=3)

    s = repr(model)

    assert "BlackShifted" in s
    assert "VOLATILITY" in s
    assert "SHIFT" in s
    assert "IMPLEMENTATION" in s
