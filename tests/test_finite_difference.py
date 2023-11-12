from financepy.models.finite_difference import (
    black_scholes_fd, dx, dxx, solve_tridiagonal_matrix, band_matrix_multiplication)
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.global_vars import gDaysInYear
from financepy.models.equity_crr_tree import crr_tree_val_avg

import numpy as np
from pytest import approx


def test_black_scholes_fd():
    """
    Compare the output of black_schole_finite_difference to kBlack::fdRunner from
    https://github.com/domokane/CompFin/blob/main/Week%204/xladdin/Utility/kBlack.cpp

    Results are identical to 3dp. The residual is due to differences interest and
    dividend rates determined from the discount and dividend curves.
    """
    s0 = 1
    r = 0.04
    dividend_yield = 0.07
    volatility = 0.2

    time_to_expiry = 5
    strike = 1.025
    dig = False
    smooth = False

    theta = 0.5
    wind = 0
    num_std = 5
    num_t = 50
    num_s = 200
    update = False

    option_type = OptionTypes.EUROPEAN_CALL
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.07939664662902503, abs=1e-3)

    smooth = True
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.07945913698961202, abs=1e-3)
    smooth = 0

    dig = 1
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.2153451094307548, abs=1e-3)

    #smooth dig
    smooth = 1
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.22078914857802928, abs=1e-3)
    smooth = 0
    dig = 0

    option_type = OptionTypes.EUROPEAN_PUT
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.2139059947533305, abs=1e-3)

    option_type = OptionTypes.AMERICAN_PUT
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.2165916613669189, abs=1e-3)

    option_type = OptionTypes.AMERICAN_CALL
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.10259475990431438, abs=1e-3)
    option_type = OptionTypes.EUROPEAN_CALL

    wind = 1
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.07834108133101789, abs=1e-3)

    wind = 2
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.08042112779963827, abs=1e-3)

    wind = -1
    v = black_scholes_fd(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                                        strike_price=strike, risk_free_rate=r,
                                        dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                                        smooth=smooth, theta=theta, wind=wind,
                                        num_std=num_std, num_time_steps=num_t, num_samples=num_s, update=update)
    assert v == approx(0.08042112779963827, abs=1e-3)
    wind = 0


def test_european_call():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.40

    value_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    num_steps_per_year = 20000
    strike_price = 50.0
    option_type = OptionTypes.EUROPEAN_CALL

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=strike_price, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=5, num_time_steps=5000, num_samples=10000, update=False)
    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v == approx(value['value'], abs=1e-3)


def test_european_put():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.40

    value_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    num_steps_per_year = 20000
    strike_price = 50.0
    option_type = OptionTypes.EUROPEAN_PUT

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=strike_price, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=5, num_time_steps=2500, num_samples=10000, update=False)
    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v == approx(value['value'], abs=1e-3)


def test_american_call():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.05
    volatility = 0.40

    value_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    num_steps_per_year = 20000
    strike_price = 50.0
    option_type = OptionTypes.AMERICAN_CALL

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=strike_price, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=6, num_time_steps=2500, num_samples=10000, update=False)

    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v == approx(value['value'], abs=1e-3)


def test_american_put():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.05
    volatility = 0.40

    value_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2021)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    num_steps_per_year = 20000
    strike_price = 50.0
    option_type = OptionTypes.AMERICAN_PUT

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=strike_price, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=5, num_time_steps=2500, num_samples=10000, update=False)
    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v == approx(value['value'], abs=1e-3)


def test_call_option():
    """
    Check finite difference method gives similar result to BlackScholes model
    """
    expiry_date = Date(1, 7, 2015)
    strike_price = 100.0
    option_type = OptionTypes.EUROPEAN_CALL
    call_option = EquityVanillaOption(
        expiry_date, strike_price, option_type)

    value_date = Date(1, 1, 2015)
    spot_price = 100
    volatility = 0.30
    risk_free_rate = 0.05
    dividend_yield = 0.01
    model = BlackScholes(volatility)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    discount_curve = DiscountCurveFlat(value_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(value_date, dividend_yield)

    # Call option
    v0 = call_option.value(value_date, spot_price,
                           discount_curve, dividend_curve, model)

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=100.0, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=5, num_time_steps=2500, num_samples=10000, update=False)
    assert v == approx(v0, 1e-5)


def test_put_option():
    """
    Check finite difference method gives similar result to BlackScholes model
    """
    expiry_date = Date(1, 7, 2015)
    strike_price = 100.0
    option_type = OptionTypes.EUROPEAN_PUT
    put_option = EquityVanillaOption(
        expiry_date, strike_price, option_type)

    value_date = Date(1, 1, 2015)
    spot_price = 100
    volatility = 0.30
    risk_free_rate = 0.05
    dividend_yield = 0.1
    model = BlackScholes(volatility)
    time_to_expiry = (expiry_date - value_date) / gDaysInYear
    discount_curve = DiscountCurveFlat(value_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(value_date, dividend_yield)

    # Call option
    v0 = put_option.value(value_date, spot_price,
                          discount_curve, dividend_curve, model)

    v = black_scholes_fd(spot_price=spot_price, volatility=volatility,
                                        time_to_expiry=time_to_expiry,
                                        strike_price=100.0, risk_free_rate=risk_free_rate,
                                        dividend_yield=dividend_yield, digital=0,
                                        option_type=option_type, smooth=0, theta=0.5, wind=0,
                                        num_std=5, num_time_steps=2500, num_samples=10000, update=False)

    assert v == approx(v0, 1e-5)


def test_dx():
    np.testing.assert_array_equal(dx([0, 1, 2, 3, 4, 5], wind=0),
                                  np.array([[0.,  -1.,  1.],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-0.5,  0.,  0.5],
                                            [-1.,  1.,  0.]]))
    np.testing.assert_array_almost_equal(dx([0, 1, 1.5, 3, 5, 10], wind=0),
                                         np.array([[0., -1, 1.],
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
