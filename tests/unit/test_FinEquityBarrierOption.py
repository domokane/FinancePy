# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.products.equity.equity_barrier_option import BarrierTypes
from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import GBMNumericalSchemeTypes


value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 80.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
b = 110.0
k = 100.0
opt_type = BarrierTypes.DOWN_AND_OUT_CALL
notional = 1.0

drift = interest_rate - dividend_yield
scheme = GBMNumericalSchemeTypes.NORMAL
process_type = ProcessTypes.GBM_PROCESS

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

model = BlackScholes(volatility)

num_obs_per_year = 100

########################################################################################


def test_down_and_out_call():

    opt_type = BarrierTypes.DOWN_AND_OUT_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 0.000

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 0.000


########################################################################################


def test_down_and_in_call():

    opt_type = BarrierTypes.DOWN_AND_IN_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 1.531

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 1.551


########################################################################################


def test_up_and_out_call():

    opt_type = BarrierTypes.UP_AND_OUT_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 0.179

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 0.164


########################################################################################


def test_up_and_in_call():

    opt_type = BarrierTypes.UP_AND_IN_CALL
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 1.352

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 1.356


########################################################################################


def test_up_and_out_put():

    opt_type = BarrierTypes.UP_AND_OUT_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 18.145

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 18.133


########################################################################################


def test_up_and_in_put():

    opt_type = BarrierTypes.UP_AND_IN_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 3) == 0.093

    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 3) == 0.099


########################################################################################


def test_down_and_out_put():

    opt_type = BarrierTypes.DOWN_AND_OUT_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 0.0000
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 0.0000


########################################################################################


def test_down_and_in_put():

    opt_type = BarrierTypes.DOWN_AND_IN_PUT
    option = EquityBarrierOption(expiry_dt, k, opt_type, b, num_obs_per_year)

    value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    assert round(value, 4) == 18.2378
    value_mc = option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(value_mc, 4) == 18.2472


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
