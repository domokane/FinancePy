# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest
from pytest import approx

from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.models.black_scholes import BlackScholesTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_types import OptionTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date


value_dt = Date(8, 5, 2015)
expiry_dt = Date(15, 1, 2016)

strike_price = 130.0
stock_price = 127.62
volatility = 0.20
interest_rate = 0.001
dividend_yield = 0.0163

opt_type = OptionTypes.AMERICAN_CALL
eu_option_type = OptionTypes.EUROPEAN_CALL

am_option = EquityAmericanOption(expiry_dt, strike_price, opt_type)

ameu_option = EquityAmericanOption(expiry_dt, strike_price, eu_option_type)

eu_option = EquityVanillaOption(expiry_dt, strike_price, eu_option_type)

discount_curve = FlatDiscountCurve(
    value_dt, interest_rate, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
)

dividend_curve = FlatDiscountCurve(
    value_dt, dividend_yield, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
)

num_steps_per_year = 400

model_tree = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps_per_year)

########################################################################################


def test_black_scholes():

    v = am_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_tree
    )
    assert round(v, 4) == 6.8398

    model_approx = BlackScholes(volatility, BlackScholesTypes.BARONE_ADESI)

    v = am_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_approx
    )

    assert round(v, 4) == 6.8277

    v = ameu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_tree
    )

    assert round(v, 4) == 6.7512

    v = eu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model_tree
    )

    assert round(v, 4) == 6.7493


########################################################################################


def test_barone_adesi_zero_and_near_zero_rates():
    """The BAW rate ratio has a finite limit as the risk-free rate tends to zero."""

    baw_model = BlackScholes(0.20, BlackScholesTypes.BARONE_ADESI)
    tree_model = BlackScholes(0.20, BlackScholesTypes.CRR_TREE, 2000)
    analytical_model = BlackScholes(0.20, BlackScholesTypes.ANALYTICAL)

    call_at_zero = baw_model.value(
        100.0, 1.0, 100.0, 0.0, 0.02, OptionTypes.AMERICAN_CALL
    )
    call_near_zero = baw_model.value(
        100.0, 1.0, 100.0, 1.0e-12, 0.02, OptionTypes.AMERICAN_CALL
    )
    call_tree = tree_model.value(
        100.0, 1.0, 100.0, 0.0, 0.02, OptionTypes.AMERICAN_CALL
    )

    assert call_at_zero == approx(call_near_zero, abs=1.0e-6)
    assert call_at_zero == approx(call_tree, abs=2.0e-2)

    put_at_zero = baw_model.value(
        100.0, 1.0, 100.0, 0.0, 0.0, OptionTypes.AMERICAN_PUT
    )
    put_near_zero = baw_model.value(
        100.0, 1.0, 100.0, 1.0e-12, 0.0, OptionTypes.AMERICAN_PUT
    )
    european_put = analytical_model.value(
        100.0, 1.0, 100.0, 0.0, 0.0, OptionTypes.EUROPEAN_PUT
    )

    assert put_at_zero == approx(european_put, abs=1.0e-12)
    assert put_near_zero == approx(put_at_zero, abs=1.0e-9)


########################################################################################


def test_bjerksund_stensland():

    # Valuation of American call option as in Bjerksund and Sensland's paper published in 1993.
    # See Table 2 in https://www.sciencedirect.com/science/article/abs/pii/095652219390009H

    # value_dt and exipry_dt are set so that time to maturity becomes 0.25

    value_dt = Date(8, 5, 2015)
    expiry_dt = Date(7, 8, 2015, hh=6)
    interest_rate = 0.08
    volatility = 0.40
    borrow_rate = 0.04
    strike_price = 100
    stock_prices = [80.0, 90.0, 100.0, 110.0, 120.0]

    # model setting
    discount_curve = FlatDiscountCurve(
        value_dt,
        interest_rate,
        FrequencyTypes.CONTINUOUS,
        DayCountTypes.ACT_365F,
    )

    borrow_curve = FlatDiscountCurve(
        value_dt,
        borrow_rate,
        FrequencyTypes.CONTINUOUS,
        DayCountTypes.ACT_365F,
    )

    model = BlackScholes(volatility, BlackScholesTypes.BJERKSUND_STENSLAND)

    # american call case
    am_call_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.AMERICAN_CALL
    )
    values = []

    for stock_price in stock_prices:

        value = am_call_option.value(
            value_dt, stock_price, discount_curve, borrow_curve, model
        )

        values.append(round(value, 2))

    assert values == [1.29, 3.82, 8.35, 14.80, 22.71]

    # american put case
    am_put_option = EquityAmericanOption(
        expiry_dt, strike_price, OptionTypes.AMERICAN_PUT
    )

    values = []

    for stock_price in stock_prices:

        value = am_put_option.value(
            value_dt, stock_price, discount_curve, borrow_curve, model
        )

        values.append(round(value, 2))

    assert values == [20.53, 12.91, 7.42, 3.93, 1.93]


########################################################################################


def test_black_scholes_fd():
    """
    Assert finite difference model matches tree model to at least 1 dp
    """
    params = {"num_samples": 200, "theta": 0.5}
    model = BlackScholes(
        volatility, bs_type=BlackScholesTypes.FINITE_DIFFERENCE, params=params
    )

    v = am_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert v == approx(6.8391, 1e-1)

    v = ameu_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert v == approx(6.7510, 1e-1)

    v = eu_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert v == approx(6.7493, 1e-1)


########################################################################################


def test_bjerksund_stensland_zero_dividend_equals_european():
    """With no dividend yield an American call is never exercised early, so the
    Bjerksund-Stensland value must equal the European value instead of failing on
    a zero division in the trigger price."""
    from financepy.models.black_scholes_analytic import (
        bjerksund_stensland_value,
        european_value,
    )

    s, t, k, r, v = 100.0, 1.0, 105.0, 0.05, 0.30
    call = OptionTypes.AMERICAN_CALL.value
    put = OptionTypes.AMERICAN_PUT.value
    euro_call = european_value(s, t, k, r, 0.0, v, OptionTypes.EUROPEAN_CALL.value)

    assert bjerksund_stensland_value(s, t, k, r, 0.0, v, call) == pytest.approx(
        euro_call
    )
    # Negative dividend yield (cost of carry above the riskless rate)
    assert bjerksund_stensland_value(s, t, k, r, -0.01, v, call) == pytest.approx(
        european_value(s, t, k, r, -0.01, v, OptionTypes.EUROPEAN_CALL.value)
    )
    # The put-call transformation maps a put with a zero rate to the same case
    euro_put = european_value(s, t, k, 0.0, 0.02, v, OptionTypes.EUROPEAN_PUT.value)
    assert bjerksund_stensland_value(s, t, k, 0.0, 0.02, v, put) == pytest.approx(
        euro_put
    )


def test_bjerksund_stensland_never_below_european():
    """The approximation is a lower bound on the American value and must not fall
    below the European value as the dividend yield approaches zero."""
    from financepy.models.black_scholes_analytic import (
        bjerksund_stensland_value,
        european_value,
    )

    s, t, k, r, v = 100.0, 1.0, 105.0, 0.05, 0.30
    for q in (1.0e-12, 1.0e-6, 1.0e-3, 0.01, 0.03, 0.05, 0.08):
        american = bjerksund_stensland_value(
            s, t, k, r, q, v, OptionTypes.AMERICAN_CALL.value
        )
        european = european_value(s, t, k, r, q, v, OptionTypes.EUROPEAN_CALL.value)
        assert american >= european - 1.0e-12
