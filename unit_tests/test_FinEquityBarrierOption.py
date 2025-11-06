# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.products.equity.equity_barrier_option import EquityBarrierTypes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 80.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
b = 110.0
k = 100.0
opt_type = EquityBarrierTypes.DOWN_AND_OUT_CALL
notional = 1.0

drift = interest_rate - dividend_yield
scheme = FinGBMNumericalScheme.NORMAL_SCHEME
process_type = ProcessTypes.GBM_PROCESS

discount_curve = DiscountCurveFlat(value_dt, interest_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

model = BlackScholes(volatility)

num_obs_per_year = 100

########################################################################################


def test_down_and_out_call():

    opt_type = EquityBarrierTypes.DOWN_AND_OUT_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 0.0000
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 0.0000


########################################################################################


def test_down_and_in_call():

    opt_type = EquityBarrierTypes.DOWN_AND_IN_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 1.5307
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 1.5483


########################################################################################


def test_up_and_out_call():

    opt_type = EquityBarrierTypes.UP_AND_OUT_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 0.1789
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 0.1572


########################################################################################


def test_up_and_in_call():

    opt_type = EquityBarrierTypes.UP_AND_IN_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 1.3519
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 1.3234


########################################################################################


def test_up_and_out_put():

    opt_type = EquityBarrierTypes.UP_AND_OUT_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 18.1445
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 18.2171


########################################################################################


def test_up_and_in_put():

    opt_type = EquityBarrierTypes.UP_AND_IN_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 0.0933
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 0.0937


########################################################################################


def test_down_and_out_put():

    opt_type = EquityBarrierTypes.DOWN_AND_OUT_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 0.0000
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 0.0000


########################################################################################


def test_down_and_in_put():

    opt_type = EquityBarrierTypes.DOWN_AND_IN_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 18.2378
    t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR
    model_params = (stock_price, drift, volatility, scheme)
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 18.2778


########################################################################################

########################################################################################

if __name__ == "__main__":
    test_down_and_in_put()
    test_down_and_in_call()
    test_up_and_in_call()
    test_up_and_in_put()
    test_down_and_out_put()
    test_down_and_out_call()
    test_up_and_out_call()
    test_up_and_out_put()
