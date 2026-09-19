# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.fx.fx_barrier_option import FXBarrierOption
from financepy.models.black_scholes import BlackScholes
from financepy.utils.global_types import GBMNumericalSchemeTypes
from financepy.utils.global_types import ProcessTypes
from financepy.utils.global_types import FXBarrierTypes


value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
currency_pair = "USDJPY"
volatility = 0.20
dom_interest_rate = 0.05
for_interest_rate = 0.02

notional = 100.0
notional_currency = "USD"

drift = dom_interest_rate - for_interest_rate
scheme = GBMNumericalSchemeTypes.ANTITHETIC
process_type = ProcessTypes.GBM_PROCESS
domestic_curve = FlatDiscountCurve(value_dt, dom_interest_rate)
foreign_curve = FlatDiscountCurve(value_dt, for_interest_rate)
model = BlackScholes(volatility)
num_obs_per_year = 100

b = 105.0
k = 100.0

########################################################################################


def test_down_and_out_call():

    spot_fx_rate = 50
    opt_type = FXBarrierTypes.DOWN_AND_OUT_CALL

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 0.0000
    assert round(delta, 4) == 0.0000
    assert round(vega, 4) == 0.0000
    assert round(theta, 4) == 0.0000


########################################################################################


def test_down_and_in_call():

    spot_fx_rate = 100
    opt_type = FXBarrierTypes.DOWN_AND_IN_CALL

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 9.2270
    assert round(delta, 4) == 0.5869
    assert round(vega, 4) == 0.3791
    assert round(theta, 4) == -5.0924


########################################################################################


def test_up_and_out_call():

    spot_fx_rate = 50
    opt_type = FXBarrierTypes.UP_AND_OUT_CALL

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 0.0003
    assert round(delta, 4) == 0.0001
    assert round(vega, 4) == 0.0002
    assert round(theta, 4) == -0.0015


########################################################################################


def test_up_and_in_call():

    spot_fx_rate = 100
    opt_type = FXBarrierTypes.UP_AND_IN_CALL

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 9.2067
    assert round(delta, 4) == 0.5900
    assert round(vega, 4) == 0.3811
    assert round(theta, 4) == -5.1228


########################################################################################


def test_up_and_out_put():

    spot_fx_rate = 50
    opt_type = FXBarrierTypes.UP_AND_OUT_PUT

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 46.1145
    assert round(delta, 4) == -0.9796
    assert round(vega, 4) == 0.0013
    assert round(theta, 4) == 3.7654


########################################################################################


def test_up_and_in_put():

    spot_fx_rate = 100
    opt_type = FXBarrierTypes.UP_AND_IN_PUT

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 2.7444
    assert round(delta, 4) == 0.2102
    assert round(vega, 4) == 0.2680
    assert round(theta, 4) == -2.3502


########################################################################################


def test_down_and_out_put():

    spot_fx_rate = 50
    opt_type = FXBarrierTypes.DOWN_AND_OUT_PUT

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 0.0000
    assert round(delta, 4) == 0.0000
    assert round(vega, 4) == 0.0000
    assert round(theta, 4) == 0.0000


########################################################################################


def test_down_and_in_put():

    spot_fx_rate = 100
    opt_type = FXBarrierTypes.DOWN_AND_IN_PUT

    barrier_option = FXBarrierOption(
        expiry_dt,
        k,
        currency_pair,
        opt_type,
        b,
        num_obs_per_year,
        notional,
        notional_currency,
    )

    value = barrier_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    delta = barrier_option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    vega = barrier_option.vega(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    theta = barrier_option.theta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert round(value, 4) == 6.3301
    assert round(delta, 4) == -0.3933
    assert round(vega, 4) == 0.3791
    assert round(theta, 4) == -2.2964
