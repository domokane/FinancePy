from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.global_vars import gDaysInYear
from financepy.models.equity_crr_tree import crr_tree_val_avg
from financepy.models.equity_lsmc import equity_lsmc, FIT_TYPES
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

from pytest import approx


def test_american_call():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    strike_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.40
    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    # TODO MAKE WORK WITH EUROPEAN OPTIONS
    option_type = OptionTypes.AMERICAN_CALL

    time_to_expiry = (expiry_date - valuation_date) / gDaysInYear
    num_steps_per_year = 500
    num_paths = 50_000
    poly_degree = 5

    v_ls = equity_lsmc(spot_price, risk_free_rate, dividend_yield, volatility, num_steps_per_year, num_paths,
                       time_to_expiry, option_type.value, strike_price, poly_degree, FIT_TYPES.LAGUERRE.value, 0, 0)

    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             50, #num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v_ls == approx(value['value'], abs=1e-1)


def test_american_put():
    """
    Check finite difference method gives similar result to binomial tree
    """
    spot_price = 50.0
    strike_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.00
    volatility = 0.20
    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    # TODO MAKE WORK WITH EUROPEAN OPTIONS
    option_type = OptionTypes.AMERICAN_PUT

    time_to_expiry = (expiry_date - valuation_date) / gDaysInYear
    num_steps_per_year = 500
    num_paths = 50_000
    poly_degree = 5

    v_ls = equity_lsmc(spot_price, risk_free_rate, dividend_yield, volatility, num_steps_per_year, num_paths,
                       time_to_expiry, option_type.value, strike_price, poly_degree, FIT_TYPES.LAGUERRE.value, 0, 0)

    value = crr_tree_val_avg(spot_price,
                             risk_free_rate,  # continuously compounded
                             dividend_yield,  # continuously compounded
                             volatility,  # Black scholes volatility
                             50,  # num_steps_per_year,
                             time_to_expiry,
                             option_type.value,
                             strike_price)
    assert v_ls == approx(value['value'], abs=1e-1)


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
    spot_price = 100
    volatility = 0.20
    risk_free_rate = 0.05
    dividend_yield = 0.01
    num_steps_per_year = 500
    num_paths = 50_000
    poly_degree = 5
    model = BlackScholes(volatility)
    time_to_expiry = (expiry_date - valuation_date) / gDaysInYear
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    # Call option
    v0 = call_option.value(valuation_date, spot_price,
                           discount_curve, dividend_curve, model)

    v_ls = equity_lsmc(spot_price, risk_free_rate, dividend_yield, volatility, num_steps_per_year, num_paths,
                       time_to_expiry, option_type.value, strike_price, poly_degree, FIT_TYPES.LAGUERRE.value, 0, 0)
    assert v_ls == approx(v0, 1e-1)


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
    spot_price = 100
    volatility = 0.20
    risk_free_rate = 0.05
    dividend_yield = 0.1
    num_steps_per_year = 500
    num_paths = 50_000
    poly_degree = 5
    model = BlackScholes(volatility)
    time_to_expiry = (expiry_date - valuation_date) / gDaysInYear
    discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    # Call option
    v0 = put_option.value(valuation_date, spot_price,
                          discount_curve, dividend_curve, model)
    v_ls = equity_lsmc(spot_price, risk_free_rate, dividend_yield, volatility, num_steps_per_year, num_paths,
                       time_to_expiry, option_type.value, strike_price, poly_degree, FIT_TYPES.LAGUERRE.value, 0, 0)

    assert v_ls == approx(v0, 1e-1)
