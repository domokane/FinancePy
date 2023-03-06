from financepy.models.finite_difference import black_scholes_finite_difference
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.products.equity.equity_binomial_tree import EquityTreePayoffTypes
from financepy.products.equity.equity_binomial_tree import EquityTreeExerciseTypes
from financepy.products.equity.equity_binomial_tree import EquityBinomialTree

import numpy as np
from pytest import approx


def test_black_scholes_finite_difference():
    s0 = 1
    r = 0.04
    dividend_yield = 0.07
    sigma = 0.2

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(30, 12, 2020)
    strike = 1.025
    dig = 0
    smooth = 0

    theta = 0.5
    wind = 0
    num_std = 5
    num_t = 50
    num_s = 200
    update = 0
    num_pr = 1

    option_type = OptionTypes.EUROPEAN_CALL
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                             num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07939664662902503)
    
    smooth = True
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07945913698961202)
    smooth = 0

    dig = 1
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2153451094307548)

    #smooth dig
    smooth = 1
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.22078914857802928)
    smooth = 0
    dig = 0

    option_type = OptionTypes.EUROPEAN_PUT
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2139059947533305)

    option_type = OptionTypes.AMERICAN_PUT
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2165916613669189)

    option_type = OptionTypes.AMERICAN_CALL
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.10259475990431438)
    option_type = OptionTypes.EUROPEAN_CALL

    wind = 1
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07834108133101789)

    wind = 2
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.08042112779963827)

    wind = -1
    _, v = black_scholes_finite_difference(s0, r, dividend_yield, sigma, expiry_date, valuation_date, strike, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.08042112779963827)
    wind = 0


def test_european_call():
    stock_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.40

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
    num_steps = 100
    strike_price = 50.0
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    option_type = OptionTypes.EUROPEAN_CALL
    params = np.array([1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=risk_free_rate,
                                           dividend_yield=dividend_yield, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False, num_pr=1)
    tree = EquityBinomialTree()
    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)  # price, delta, gamma, theta
    assert v == approx(value[0], abs=1e-1)


def test_european_put():
    stock_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.40

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
    num_steps = 100
    strike_price = 50.0
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    option_type = OptionTypes.EUROPEAN_PUT
    params = np.array([-1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=risk_free_rate,
                                           dividend_yield=dividend_yield, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False, num_pr=1)
    tree = EquityBinomialTree()
    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)  # price, delta, gamma, theta
    assert v == approx(value[0], abs=1e-1)


def test_call_option():
    expiry_date = Date(1, 7, 2015)
    strike_price = 100.0
    option_type = OptionTypes.EUROPEAN_CALL
    call_option = EquityVanillaOption(
        expiry_date, strike_price, option_type)

    valuation_date = Date(1, 1, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    # Call option
    v0 = call_option.value(valuation_date, stock_price,
                           discount_curve, dividend_curve, model)

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=interest_rate,
                                           dividend_yield=dividend_yield, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=100.0, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False, num_pr=1)
    assert v == approx(v0, 1e-1)


def test_put_option():
    expiry_date = Date(1, 7, 2015)
    strike_price = 100.0
    option_type = OptionTypes.EUROPEAN_PUT
    put_option = EquityVanillaOption(
        expiry_date, strike_price, option_type)

    valuation_date = Date(1, 1, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.1
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    # Call option
    v0 = put_option.value(valuation_date, stock_price,
                          discount_curve, dividend_curve, model)

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=interest_rate,
                                           dividend_yield=dividend_yield, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=100.0, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False, num_pr=1)

    assert v == approx(v0, 1e-1)
