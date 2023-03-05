from pytest import approx

from financepy.models.finite_difference import black_scholes_finite_difference, PUT_CALL, exercise_type
from financepy.utils.global_vars import gDaysInYear

def test_black_scholes_finite_difference():
    s0 = 1
    r = 0.04
    mu = -0.03
    sigma = 0.2

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(30, 12, 2020)
    strike = 1.025
    dig = 0
    pc = PUT_CALL.CALL.value
    ea = exercise_type.EUROPEAN.value
    smooth = 0

    theta = 0.5
    wind = 0
    num_std = 5
    num_t = 50
    num_s = 200
    update = 0
    num_pr = 1

    # European call
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                             num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07939664662902503)
    
    # smooth
    smooth = 1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07945913698961202)
    smooth = 0

    # dig
    dig = 1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2153451094307548)

    #smooth dig
    smooth = 1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.22078914857802928)
    smooth = 0
    dig = 0

    # European put
    pc = PUT_CALL.PUT.value
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2139059947533305)

    # American put
    ea = exercise_type.AMERICAN.value
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2165916613669189)

    # American call
    pc = PUT_CALL.CALL.value
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.10259475990431438)
    ea = exercise_type.EUROPEAN.value

    # wind=1
    wind = 1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07834108133101789)

    # wind=2
    wind = 2
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.08042112779963827)

    # wind=-1
    wind = -1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.08042112779963827)
    wind = 0


from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.products.equity.equity_binomial_tree import EquityTreePayoffTypes
from financepy.products.equity.equity_binomial_tree import EquityTreeExerciseTypes
from financepy.products.equity.equity_binomial_tree import EquityBinomialTree
import numpy as np

stock_price = 50.0
risk_free_rate = 0.06
dividend_yield = 0.00
volatility = 0.40

valuation_date = Date(1, 1, 2016)
expiry_date = Date(1, 1, 2021)
model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)
num_steps = 100
strike_price = 50.0



def test_european_call():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=risk_free_rate,
                                           mu=0, sigma=np.sqrt(volatility),
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, dig=0,
                                           pc=PUT_CALL.CALL.value, exercise=exercise, smooth=0, theta=0.5, wind=0,
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
    assert v == round(value[0], 4)

def test_european_put():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1.0, strike_price])

    _, v = black_scholes_finite_difference(stock_price=stock_price, risk_free_rate=risk_free_rate,
                                           mu=0, sigma=np.sqrt(volatility),
                                           expiry_date=expiry_date, valuation_date=valuation_date,
                                           strike_price=strike_price, dig=0,
                                           pc=PUT_CALL.PUT.value, exercise=exercise, smooth=0, theta=0.5, wind=0,
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
    assert v == round(value[0], 4)

