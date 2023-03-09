from financepy.models.finite_difference import (
    black_scholes_finite_difference, dx, dxx, solve_tridiagonal_matrix, band_matrix_multiplication)
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
    """
    Compare the output of black_schole_finite_difference to kBlack::fdRunner from
    https://github.com/domokane/CompFin/blob/main/Week%204/xladdin/Utility/kBlack.cpp

    Results are identical to 3dp. The residual is due to differences interest and
    dividend rates determined from the discount and dividend curves.
    """
    s0 = 1
    r = 0.04
    dividend_yield = 0.07
    sigma = 0.2

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(30, 12, 2020)
    discount_curve = DiscountCurveFlat(valuation_date, r)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
    strike = 1.025
    dig = 0
    smooth = 0

    theta = 0.5
    wind = 0
    num_std = 5
    num_t = 50
    num_s = 200
    update = 0

    option_type = OptionTypes.EUROPEAN_CALL
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.07939664662902503, abs=1e-3)

    smooth = True
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.07945913698961202, abs=1e-3)
    smooth = 0

    dig = 1
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.2153451094307548, abs=1e-3)

    #smooth dig
    smooth = 1
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.22078914857802928, abs=1e-3)
    smooth = 0
    dig = 0

    option_type = OptionTypes.EUROPEAN_PUT
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.2139059947533305, abs=1e-3)

    option_type = OptionTypes.AMERICAN_PUT
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.2165916613669189, abs=1e-3)

    option_type = OptionTypes.AMERICAN_CALL
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.10259475990431438, abs=1e-3)
    option_type = OptionTypes.EUROPEAN_CALL

    wind = 1
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.07834108133101789, abs=1e-3)

    wind = 2
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.08042112779963827, abs=1e-3)

    wind = -1
    _, v = black_scholes_finite_difference(s0, sigma, expiry_date, valuation_date, strike, discount_curve,
                                           dividend_curve, dig, option_type, smooth, theta, wind,
                                           num_std, num_t, num_s, update)
    assert v == approx(0.08042112779963827, abs=1e-3)
    wind = 0


def test_european_call():
    """
    Check finite difference method gives similar result to binomial tree
    """
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

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)
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
    """
    Check finite difference method gives similar result to binomial tree
    """
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

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)
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


def test_american_call():
    """
    Check finite difference method gives similar result to binomial tree
    """
    stock_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.05
    volatility = 0.40

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
    num_steps = 100
    strike_price = 50.0
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    option_type = OptionTypes.AMERICAN_CALL
    params = np.array([1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)
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


def test_american_put():
    """
    Check finite difference method gives similar result to binomial tree
    """
    stock_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.05
    volatility = 0.40

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
    num_steps = 100
    strike_price = 50.0
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    option_type = OptionTypes.AMERICAN_PUT
    params = np.array([-1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)
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
    """
    Check finite difference method gives similar result to BlackScholes model
    """
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

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=100.0, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)
    assert v == approx(v0, 1e-1)


def test_put_option():
    """
    Check finite difference method gives similar result to BlackScholes model
    """
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

    _, v = black_scholes_finite_difference(stock_price=stock_price, sigma=volatility,
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=100.0, discount_curve=discount_curve,
                                           dividend_curve=dividend_curve, digital=0,
                                           option_type=option_type, smooth=0, theta=0.5, wind=0,
                                           num_std=5, num_steps=50, num_samples=200, update=False)

    assert v == approx(v0, 1e-1)


def test_dx():
    np.testing.assert_array_equal(dx([0, 1, 2, 3, 4, 5], wind=0),
                                  np.array([[-1.,  0.,  1.],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-1.,  1.,  0.]]))
    np.testing.assert_array_almost_equal(dx([0, 1, 1.5, 3, 5, 10], wind=0),
                                         np.array([[-1., 0., 1.],
                                                   [-0.33333333, -1., 1.33333333],
                                                   [-1.5, 1.33333333, 0.16666667],
                                                   [-0.38095238, 0.16666667, 0.21428571],
                                                   [-0.35714286, 0.3, 0.05714286],
                                                   [-0.2, 0.2, 0.]]), decimal=3)


def test_dxx():
    np.testing.assert_array_equal(dxx([1, 1.5, 2, 2.5, 3]),
                                  np.array([[0.,  0.,  0.],
                                            [4., -8.,  4.],
                                            [4., -8.,  4.],
                                            [4., -8.,  4.],
                                            [0.,  0.,  0.]]))


def test_solve_tridiagonal_matrix():
    M = np.array([
        [0, 1, 1, 1],
        [-2, -2, -2, -2],
        [1, 1, 1, 0]]
    ).T

    r = np.array([1, 1, 1, 1]) * 0.04
    u = solve_tridiagonal_matrix(M, r)

    np.testing.assert_array_equal(u, np.array([-0.08, -0.12, -0.12, -0.08]))


def test_band_matrix_multiplication():
    M = np.array([
        [0, 1, 1, 1],
        [-2, -2, -2, -2],
        [1, 1, 1, 0]]
    ).T
    u = np.array([-0.08, -0.12, -0.12, -0.08])
    np.testing.assert_array_almost_equal(band_matrix_multiplication(M, 1, 1, u), np.array([0.04] * 4))
