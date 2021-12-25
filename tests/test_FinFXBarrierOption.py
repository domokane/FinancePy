###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.fx.fx_barrier_option import FXBarrierOption
from financepy.products.fx.fx_barrier_option import FinFXBarrierTypes
from financepy.models.black_scholes import BlackScholes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes


valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
currency_pair = "USDJPY"
volatility = 0.20
dom_interest_rate = 0.05
forInterestRate = 0.02

notional = 100.0
notional_currency = "USD"

drift = dom_interest_rate - forInterestRate
scheme = FinGBMNumericalScheme.ANTITHETIC
process_type = ProcessTypes.GBM
dom_discount_curve = DiscountCurveFlat(valuation_date, dom_interest_rate)
for_discount_curve = DiscountCurveFlat(valuation_date, forInterestRate)
model = BlackScholes(volatility)
num_observations_per_year = 100

B = 105.0
K = 100.0


def test_DOWN_AND_OUT_CALL():
    spot_fx_rate = 50
    option_type = FinFXBarrierTypes.DOWN_AND_OUT_CALL

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 0.0000
    assert round(delta, 4) == 0.0000
    assert round(vega, 4) == 0.0000
    assert round(theta, 4) == 0.0000


def test_DOWN_AND_IN_CALL():
    spot_fx_rate = 100
    option_type = FinFXBarrierTypes.DOWN_AND_IN_CALL

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 9.2270
    assert round(delta, 4) == 0.5868
    assert round(vega, 4) == 0.3791
    assert round(theta, 4) == -5.0924


def test_UP_AND_OUT_CALL():
    spot_fx_rate = 50
    option_type = FinFXBarrierTypes.UP_AND_OUT_CALL

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 0.0003
    assert round(delta, 4) == 0.0001
    assert round(vega, 4) == 0.0002
    assert round(theta, 4) == -0.0015


def test_UP_AND_IN_CALL():
    spot_fx_rate = 100
    option_type = FinFXBarrierTypes.UP_AND_IN_CALL

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 9.2067
    assert round(delta, 4) == 0.5900
    assert round(vega, 4) == 0.3811
    assert round(theta, 4) == -5.1229


def test_UP_AND_OUT_PUT():
    spot_fx_rate = 50
    option_type = FinFXBarrierTypes.UP_AND_OUT_PUT

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 46.1145
    assert round(delta, 4) == -0.9796
    assert round(vega, 4) == 0.0013
    assert round(theta, 4) == 3.7654


def test_UP_AND_IN_PUT():
    spot_fx_rate = 100
    option_type = FinFXBarrierTypes.UP_AND_IN_PUT

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 2.7444
    assert round(delta, 4) == 0.2102
    assert round(vega, 4) == 0.2680
    assert round(theta, 4) == -2.3502


def test_DOWN_AND_OUT_PUT():
    spot_fx_rate = 50
    option_type = FinFXBarrierTypes.DOWN_AND_OUT_PUT

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 0.0000
    assert round(delta, 4) == 0.0000
    assert round(vega, 4) == 0.0000
    assert round(theta, 4) == 0.0000


def test_DOWN_AND_IN_PUT():
    spot_fx_rate = 100
    option_type = FinFXBarrierTypes.DOWN_AND_IN_PUT

    barrier_option = FXBarrierOption(expiry_date,
                                     K,
                                     currency_pair,
                                     option_type,
                                     B,
                                     num_observations_per_year,
                                     notional,
                                     notional_currency)

    value = barrier_option.value(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    delta = barrier_option.delta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    vega = barrier_option.vega(valuation_date,
                               spot_fx_rate,
                               dom_discount_curve,
                               for_discount_curve,
                               model)

    theta = barrier_option.theta(valuation_date,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)

    assert round(value, 4) == 6.3301
    assert round(delta, 4) == -0.3934
    assert round(vega, 4) == 0.3791
    assert round(theta, 4) == -2.2964
