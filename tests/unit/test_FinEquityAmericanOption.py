# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_american_analytic import BlackScholesAlgorithmTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_american_option import (
    EquityAmericanOption,
)


value_dt = Date(1, 1, 2016)
expiry_dt = Date(1, 1, 2017)
stock_price = 50.0
interest_rate = 0.06
dividend_yield = 0.04
volatility = 0.40
strike_price = 50.0
num_steps = 100

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

model = BlackScholes(volatility)


@pytest.mark.parametrize("spot", [50, 50.0, np.int64(50), np.float64(50), np.array(50.0)])
def test_scalar_spot_representations_return_one_price(spot):
    """Python/NumPy scalar representations have the same contract value."""
    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.AMERICAN_PUT, 3.0)
    value = option.value(value_dt, spot, discount_curve, dividend_curve, model)
    assert np.ndim(value) == 0
    assert round(value / 3.0, 3) == 7.255


########################################################################################


def test_european_put():

    put_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT
    )

    value = put_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value, 3) == 7.081


########################################################################################


def test_american_put():

    put_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.AMERICAN_PUT
    )

    value = put_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value, 3) == 7.255


########################################################################################


def test_european_call():

    call_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
    )

    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value, 3) == 8.035


########################################################################################


def test_american_call():

    call_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.AMERICAN_CALL
    )

    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value, 3) == 8.056
