from pytest import approx

from financepy.models.finite_difference import black_scholes_finite_difference, PUT_CALL, AMER_EURO
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
    ea = AMER_EURO.EURO.value
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
    ea = AMER_EURO.AMER.value
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.2165916613669189)

    # American call
    pc = PUT_CALL.CALL.value
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.10259475990431438)
    ea = AMER_EURO.EURO.value

    # wind=1
    wind = 1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.07834108133101789)

    # wind=-1
    wind = -1
    _, v = black_scholes_finite_difference(s0, r, mu, sigma, expiry_date, valuation_date, strike, dig, pc, ea, smooth, theta, wind,
                                           num_std, num_t, num_s, update, num_pr)
    assert v == approx(0.08042112779963827)
    wind = 0
