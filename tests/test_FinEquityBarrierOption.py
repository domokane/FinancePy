###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.products.equity.equity_barrier_option import EquityBarrierTypes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes


valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
stock_price = 80
B = 110.0
K = 100.0
option_type = EquityBarrierTypes.DOWN_AND_OUT_CALL

drift = interest_rate - dividend_yield
scheme = FinGBMNumericalScheme.NORMAL
process_type = ProcessTypes.GBM

discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

model = BlackScholes(volatility)

num_observations_per_year = 100


def test_down_and_out_call():
    option_type = EquityBarrierTypes.DOWN_AND_OUT_CALL
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 0.0000

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 0.0000


def test_down_and_in_call():
    option_type = EquityBarrierTypes.DOWN_AND_IN_CALL
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 1.5307

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 1.5718


def test_up_and_out_call():
    option_type = EquityBarrierTypes.UP_AND_OUT_CALL
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 0.1789

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 0.1706


def test_up_and_in_call():
    option_type = EquityBarrierTypes.UP_AND_IN_CALL
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 1.3519

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 1.4426


def test_up_and_out_put():
    option_type = EquityBarrierTypes.UP_AND_OUT_PUT
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 18.1445

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 17.8602


def test_up_and_in_put():
    option_type = EquityBarrierTypes.UP_AND_IN_PUT
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 0.0933

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 0.0940


def test_down_and_out_put():
    option_type = EquityBarrierTypes.DOWN_AND_OUT_PUT
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 0.0000

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 0.0000


def test_down_and_in_put():
    option_type = EquityBarrierTypes.DOWN_AND_IN_PUT
    option = EquityBarrierOption(
        expiry_date, K, option_type, B, num_observations_per_year)

    value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

    assert round(value, 4) == 18.2378

    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

    assert round(value_mc, 4) == 17.9996
