# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest
from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_asian_option import AsianOptionValuationTypes
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.utils.global_types import OptionTypes

value_dt = Date(1, 1, 2014)
start_averaging_dt = Date(1, 6, 2014)
expiry_dt = Date(1, 1, 2015)
stock_price = 100.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
num_observations = 120  # daily as we have a half year
accrued_avg = None
k = 100
seed = 1976
num_paths = 5000

model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

asian_option = EquityAsianOption(
    start_averaging_dt,
    expiry_dt,
    k,
    OptionTypes.EUROPEAN_CALL,
    num_observations,
)

########################################################################################


def test_geometric():

    value_geometric = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.KEMNA_VORST,
        accrued_avg,
    )

    assert round(value_geometric, 4) == 6.9723


########################################################################################


def test_turnbull_wakeman():

    value_turnbull_wakeman = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_avg,
    )

    assert round(value_turnbull_wakeman, 4) == 7.0895


########################################################################################


def test_curran():

    value_curran = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.CURRAN,
        accrued_avg,
    )

    assert round(value_curran, 4) == 7.0857


########################################################################################


def test_mc():

    value_mc = asian_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
        seed,
        accrued_avg,
    ).value

    assert round(value_mc, 3) == 7.121


SPOT_BUMP = 1.0e-4
VOL_BUMP = 0.01
RATE_BUMP = 1.0e-4

METHODS = [
    AsianOptionValuationTypes.KEMNA_VORST,
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    AsianOptionValuationTypes.CURRAN,
]


def make_market(value_dt, volatility=0.20):
    stock_price = 100.0

    discount_curve = FlatDiscountCurve(
        value_dt,
        0.05,
    )

    dividend_curve = FlatDiscountCurve(
        value_dt,
        0.02,
    )

    model = BlackScholes(volatility)

    return (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )


def make_option():
    return EquityAsianOption(
        Date(1, 1, 2015),
        Date(1, 1, 2016),
        100.0,
        OptionTypes.EUROPEAN_CALL,
        100,
    )


# ============================================================================
# DELTA
# ============================================================================


@pytest.mark.parametrize("method", METHODS)
def test_asian_delta_before_averaging(method):

    value_dt = Date(1, 7, 2014)

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    value_up = option.value(
        value_dt,
        stock_price + SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    expected = (value_up - value) / SPOT_BUMP

    actual = option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-10,
        abs=1.0e-10,
    )


@pytest.mark.parametrize("method", METHODS)
def test_asian_delta_during_averaging(method):

    value_dt = Date(1, 7, 2015)
    accrued_average = 103.0

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
        accrued_average=accrued_average,
    )

    value_up = option.value(
        value_dt,
        stock_price + SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        method=method,
        accrued_average=accrued_average,
    )

    expected = (value_up - value) / SPOT_BUMP

    actual = option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
        accrued_average=accrued_average,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-10,
        abs=1.0e-10,
    )


# ============================================================================
# GAMMA
# ============================================================================


@pytest.mark.parametrize("method", METHODS)
def test_asian_gamma_before_averaging(method):

    value_dt = Date(1, 7, 2014)

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    value_up = option.value(
        value_dt,
        stock_price + SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    value_down = option.value(
        value_dt,
        stock_price - SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    expected = (
        value_up
        - 2.0 * value
        + value_down
    ) / SPOT_BUMP**2

    actual = option.gamma(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        method=method,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-8,
        abs=1.0e-8,
    )


@pytest.mark.parametrize("method", METHODS)
def test_asian_gamma_during_averaging(method):

    value_dt = Date(1, 7, 2015)
    accrued_average = 103.0

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    kwargs = {
        "method": method,
        "accrued_average": accrued_average,
    }

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    value_up = option.value(
        value_dt,
        stock_price + SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    value_down = option.value(
        value_dt,
        stock_price - SPOT_BUMP,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    expected = (
        value_up
        - 2.0 * value
        + value_down
    ) / SPOT_BUMP**2

    actual = option.gamma(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-8,
        abs=1.0e-8,
    )


# ============================================================================
# VEGA
# ============================================================================


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "value_dt, accrued_average",
    [
        (Date(1, 7, 2014), None),
        (Date(1, 7, 2015), 103.0),
    ],
)
def test_asian_vega(method, value_dt, accrued_average):

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    kwargs = {
        "method": method,
        "accrued_average": accrued_average,
    }

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    bumped_model = BlackScholes(
        model.volatility + VOL_BUMP
    )

    value_up = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        bumped_model,
        **kwargs,
    )

    # EquityOption.vega() currently reports the value
    # change for a 1 vol-point bump, not dV/dsigma.
    expected = value_up - value

    actual = option.vega(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-10,
        abs=1.0e-10,
    )


# ============================================================================
# VANNA
# ============================================================================


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "value_dt, accrued_average",
    [
        (Date(1, 7, 2014), None),
        (Date(1, 7, 2015), 103.0),
    ],
)
def test_asian_vanna(method, value_dt, accrued_average):

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    kwargs = {
        "method": method,
        "accrued_average": accrued_average,
    }

    delta = option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    bumped_model = BlackScholes(
        model.volatility + SPOT_BUMP
    )

    delta_up = option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        bumped_model,
        **kwargs,
    )

    expected = (
        delta_up - delta
    ) / SPOT_BUMP

    actual = option.vanna(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-8,
        abs=1.0e-8,
    )


# ============================================================================
# RHO
# ============================================================================


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "value_dt, accrued_average",
    [
        (Date(1, 7, 2014), None),
        (Date(1, 7, 2015), 103.0),
    ],
)
def test_asian_rho(method, value_dt, accrued_average):

    option = make_option()

    (
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    ) = make_market(value_dt)

    kwargs = {
        "method": method,
        "accrued_average": accrued_average,
    }

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    value_up = option.value(
        value_dt,
        stock_price,
        discount_curve.bump_parallel(RATE_BUMP),
        dividend_curve,
        model,
        **kwargs,
    )

    expected = (
        value_up - value
    ) / RATE_BUMP

    actual = option.rho(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        **kwargs,
    )

    assert actual == pytest.approx(
        expected,
        rel=1.0e-8,
        abs=1.0e-8,
    )
