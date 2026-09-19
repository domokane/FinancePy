# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.models.black import (
    Black,
    BlackTypes,
    black_value,
    black_delta,
    black_gamma,
    black_vega,
    black_theta,
    implied_volatility,
)
from financepy.utils.global_types import OptionTypes
from financepy.utils.error import FinError

FWD = 100.0
STRIKE = 95.0
TEXP = 2.0
RATE = 0.04
DF = np.exp(-RATE * TEXP)
VOL = 0.20


def test_black_constructor_stores_inputs():

    model = Black(VOL, BlackTypes.ANALYTICAL, num_steps=100)

    assert model.volatility == VOL
    assert model.implementation_type is BlackTypes.ANALYTICAL
    assert model.num_steps == 100


def test_black_rejects_negative_volatility():

    with pytest.raises(FinError):
        Black(-0.01)


def test_black_call_value():

    model = Black(VOL)

    value = model.value(
        FWD,
        STRIKE,
        TEXP,
        DF,
        OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(value, 12.59471597859074, atol=1e-12)


def test_black_put_value():

    model = Black(VOL)

    value = model.value(
        FWD,
        STRIKE,
        TEXP,
        DF,
        OptionTypes.EUROPEAN_PUT,
    )

    assert np.isclose(value, 7.979134246657562, atol=1e-12)


def test_black_put_call_parity():

    model = Black(VOL)

    call = model.value(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)
    put = model.value(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_PUT)

    assert np.isclose(call - put, DF * (FWD - STRIKE), atol=1e-12)


def test_black_zero_expiry_call_intrinsic():

    model = Black(VOL)

    value = model.value(
        FWD,
        STRIKE,
        0.0,
        1.0,
        OptionTypes.EUROPEAN_CALL,
    )

    assert value == max(FWD - STRIKE, 0.0)


def test_black_zero_expiry_put_intrinsic():

    model = Black(VOL)

    value = model.value(
        FWD,
        STRIKE,
        0.0,
        1.0,
        OptionTypes.EUROPEAN_PUT,
    )

    assert value == max(STRIKE - FWD, 0.0)


def test_black_zero_vol_call_intrinsic_discounted():

    model = Black(0.0)

    value = model.value(
        FWD,
        STRIKE,
        TEXP,
        DF,
        OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(value, DF * max(FWD - STRIKE, 0.0), atol=1e-14)


def test_black_zero_vol_put_intrinsic_discounted():

    model = Black(0.0)

    value = model.value(
        FWD,
        STRIKE,
        TEXP,
        DF,
        OptionTypes.EUROPEAN_PUT,
    )

    assert np.isclose(value, DF * max(STRIKE - FWD, 0.0), atol=1e-14)


def test_black_delta_matches_finite_difference():

    model = Black(VOL)
    h = 1.0e-4

    analytic = model.delta(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    up = model.value(FWD + h, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)
    down = model.value(FWD - h, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    fd = (up - down) / (2.0 * h)

    assert np.isclose(analytic, fd, atol=1e-7)


def test_black_gamma_matches_finite_difference():

    model = Black(VOL)
    h = 1.0e-3

    analytic = model.gamma(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    up = model.value(FWD + h, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)
    mid = model.value(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)
    down = model.value(FWD - h, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    fd = (up - 2.0 * mid + down) / (h * h)

    assert np.isclose(analytic, fd, atol=5e-6)


def test_black_vega_matches_finite_difference():

    model = Black(VOL)
    h = 1.0e-5

    analytic = model.vega(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    up = Black(VOL + h).value(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)
    down = Black(VOL - h).value(FWD, STRIKE, TEXP, DF, OptionTypes.EUROPEAN_CALL)

    fd = (up - down) / (2.0 * h)

    assert np.isclose(analytic, fd, atol=1e-6)


def test_black_function_and_class_value_agree():

    model = Black(VOL)

    class_value = model.value(
        FWD,
        STRIKE,
        TEXP,
        DF,
        OptionTypes.EUROPEAN_CALL,
    )

    function_value = black_value(
        FWD,
        TEXP,
        STRIKE,
        RATE,
        VOL,
        OptionTypes.EUROPEAN_CALL,
    )

    assert np.isclose(class_value, function_value, atol=1e-14)


def test_black_implied_volatility_recovers_input_vol():

    price = black_value(
        FWD,
        TEXP,
        STRIKE,
        RATE,
        VOL,
        OptionTypes.EUROPEAN_CALL,
    )

    implied_vol = implied_volatility(
        FWD,
        TEXP,
        RATE,
        STRIKE,
        price,
        OptionTypes.EUROPEAN_CALL,
        debug_print=False,
    )

    assert np.isclose(implied_vol, VOL, atol=1e-6)


@pytest.mark.parametrize(
    "bad_args",
    [
        (-1.0, STRIKE, TEXP, DF),
        (FWD, -1.0, TEXP, DF),
        (FWD, STRIKE, -1.0, DF),
        (FWD, STRIKE, TEXP, 0.0),
        (FWD, STRIKE, TEXP, -0.1),
        (FWD, STRIKE, TEXP, 1.1),
    ],
)
def test_black_rejects_bad_inputs(bad_args):

    model = Black(VOL)

    with pytest.raises(FinError):
        model.value(*bad_args, OptionTypes.EUROPEAN_CALL)


def test_black_rejects_unsupported_option_type_for_analytical_model():

    model = Black(VOL, BlackTypes.ANALYTICAL)

    with pytest.raises(FinError):
        model.value(
            FWD,
            STRIKE,
            TEXP,
            DF,
            OptionTypes.AMERICAN_CALL,
        )


def test_black_repr():

    model = Black(VOL, BlackTypes.ANALYTICAL, num_steps=100)

    s = repr(model)

    assert "Black" in s
    assert "VOLATILITY" in s
    assert "IMPLEMENTATION" in s
    assert "NUMSTEPS" in s
