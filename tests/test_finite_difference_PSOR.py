from financepy.models.finite_difference_PSOR import black_scholes_fd_PSOR
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.global_vars import gDaysInYear
from financepy.models.equity_crr_tree import crr_tree_val_avg

from pytest import approx


def test_black_scholes_fd_PSOR():
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

    option_type = OptionTypes.EUROPEAN_CALL
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.07939664662902503, abs=1e-1)

    smooth = True
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.07945913698961202, abs=1e-1)
    smooth = 0

    dig = 1
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.2153451094307548, abs=1e-1)

    #smooth dig
    smooth = 1
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.22078914857802928, abs=1e-1)
    smooth = 0
    dig = 0

    option_type = OptionTypes.EUROPEAN_PUT
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.2139059947533305, abs=1e-1)

    option_type = OptionTypes.AMERICAN_PUT
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.2165916613669189, abs=1e-1)

    option_type = OptionTypes.AMERICAN_CALL
    v = black_scholes_fd_PSOR(spot_price=s0, volatility=volatility, time_to_expiry=time_to_expiry,
                              strike_price=strike, risk_free_rate=r,
                              dividend_yield=dividend_yield, digital=dig, option_type=option_type,
                              smooth=smooth, theta=theta,
                              )
    assert v == approx(0.10259475990431438, abs=1e-1)
    option_type = OptionTypes.EUROPEAN_CALL


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
    strike_price = 50.0
    option_type = OptionTypes.EUROPEAN_CALL
    num_steps_per_year = 2000

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=strike_price, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0, theta=0.5,
                              )
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
    num_steps_per_year = 2000
    strike_price = 50.0
    option_type = OptionTypes.EUROPEAN_PUT

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=strike_price, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0, theta=0.5,
                              )
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
    num_steps_per_year = 2000
    strike_price = 50.0
    option_type = OptionTypes.AMERICAN_CALL

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=strike_price, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0, theta=0.5,
                              num_samples=5000
                              )
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
    num_steps_per_year = 2000
    strike_price = 50.0
    option_type = OptionTypes.AMERICAN_PUT

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=strike_price, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0,
                              num_samples=5000
                              )
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

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=100.0, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0, theta=0.5,
                              )
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

    v = black_scholes_fd_PSOR(spot_price=spot_price, volatility=volatility,
                              time_to_expiry=time_to_expiry,
                              strike_price=100.0, risk_free_rate=risk_free_rate,
                              dividend_yield=dividend_yield, digital=0,
                              option_type=option_type, smooth=0, theta=0.5,
                              )

    assert v == approx(v0, 1e-5)
