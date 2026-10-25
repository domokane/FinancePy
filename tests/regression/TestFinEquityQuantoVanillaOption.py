##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

import add_fp_to_path

import math

from financepy.products.equity.equity_quanto_vanilla_option import (
    EquityQuantoVanillaOption,
)
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_quanto_vanilla_option():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    strike_price = 125.0

    domestic_rate = 0.05
    foreign_rate = 0.03
    dividend_rate = 0.02

    equity_volatility = 0.25
    fx_volatility = 0.12
    correlation = 0.50

    num_options = 100.0
    quanto_fx_rate = 1.0

    domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)
    foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

    model = BlackScholes(equity_volatility)

    test_cases.header(
        "OPTION",
        "VALUE",
        "DELTA",
        "GAMMA",
        "EQ VEGA",
        "FX VEGA",
        "CORR DELTA",
        "THETA",
        "DOM RHO",
        "FOR RHO",
    )

    for opt_type in (
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ):

        option = EquityQuantoVanillaOption(
            expiry_dt,
            strike_price,
            opt_type,
            num_options,
            quanto_fx_rate,
        )

        args = (
            value_dt,
            stock_price,
            domestic_curve,
            foreign_curve,
            dividend_curve,
            model,
            fx_volatility,
            correlation,
        )

        test_cases.print(
            str(opt_type),
            option.value(*args),
            option.delta(*args),
            option.gamma(*args),
            option.vega(*args),
            option.fx_vega(*args),
            option.correlation_delta(*args),
            option.theta(*args),
            option.domestic_rho(*args),
            option.foreign_rho(*args),
        )


########################################################################################


def test_zero_correlation():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    strike_price = 125.0

    domestic_rate = 0.05
    foreign_rate = 0.03
    dividend_rate = 0.02

    model = BlackScholes(0.25)

    domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)
    foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

    option = EquityQuantoVanillaOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    quanto_value = option.value(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        dividend_curve,
        model,
        0.12,
        0.0,
    )

    # With zero correlation, effective dividend yield is:
    # q_quanto = q + rd - rf

    adjusted_dividend_rate = (
        dividend_rate + domestic_rate - foreign_rate
    )

    adjusted_dividend_curve = FlatDiscountCurve(
        value_dt, adjusted_dividend_rate
    )

    from financepy.products.equity.equity_vanilla_option import (
        EquityVanillaOption,
    )

    vanilla_option = EquityVanillaOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    expected_value = vanilla_option.value(
        value_dt,
        stock_price,
        domestic_curve,
        adjusted_dividend_curve,
        model,
    )

    assert math.isclose(
        quanto_value,
        expected_value,
        rel_tol=1.0e-10,
    )

    test_cases.header("ZERO CORR QUANTO", "EXPECTED BS")
    test_cases.print(quanto_value, expected_value)


########################################################################################


def test_put_call_parity():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    strike_price = 125.0

    domestic_rate = 0.05
    foreign_rate = 0.03
    dividend_rate = 0.02

    sigma = 0.25
    fx_volatility = 0.12
    correlation = 0.50

    domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)
    foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

    model = BlackScholes(sigma)

    call = EquityQuantoVanillaOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    put = EquityQuantoVanillaOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_PUT,
    )

    args = (
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        dividend_curve,
        model,
        fx_volatility,
        correlation,
    )

    call_value = call.value(*args)
    put_value = put.value(*args)

    t = (expiry_dt - value_dt) / 365.0

    # Use actual curve discount factors and the same maturity
    # convention as the option's pricing model.
    from financepy.utils.helpers import option_years

    t = option_years(value_dt, expiry_dt)

    rd = domestic_curve.zero_rate_cc(expiry_dt)
    rf = foreign_curve.zero_rate_cc(expiry_dt)
    q = dividend_curve.zero_rate_cc(expiry_dt)

    q_quanto = q + rd - rf + correlation * sigma * fx_volatility

    expected = (
        stock_price * math.exp(-q_quanto * t)
        - strike_price * math.exp(-rd * t)
    )

    assert math.isclose(
        call_value - put_value,
        expected,
        rel_tol=1.0e-8,
        abs_tol=1.0e-8,
    )

    test_cases.header("CALL-PUT", "PARITY RHS")
    test_cases.print(call_value - put_value, expected)


########################################################################################


def test_quanto_greeks():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    strike_price = 125.0

    domestic_curve = FlatDiscountCurve(value_dt, 0.05)
    foreign_curve = FlatDiscountCurve(value_dt, 0.03)
    dividend_curve = FlatDiscountCurve(value_dt, 0.02)

    equity_volatility = 0.25
    fx_volatility = 0.12
    correlation = 0.50

    option = EquityQuantoVanillaOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    def price(s, sigma, fx_vol, corr):
        return option.value(
            value_dt,
            s,
            domestic_curve,
            foreign_curve,
            dividend_curve,
            BlackScholes(sigma),
            fx_vol,
            corr,
        )

    args = (
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        dividend_curve,
        BlackScholes(equity_volatility),
        fx_volatility,
        correlation,
    )

    ds = 0.01
    dv = 0.0001
    dc = 0.0001

    # Delta
    fd_delta = (
        price(stock_price + ds, equity_volatility, fx_volatility, correlation)
        - price(stock_price - ds, equity_volatility, fx_volatility, correlation)
    ) / (2.0 * ds)

    # Gamma
    fd_gamma = (
        price(stock_price + ds, equity_volatility, fx_volatility, correlation)
        - 2.0 * price(stock_price, equity_volatility, fx_volatility, correlation)
        + price(stock_price - ds, equity_volatility, fx_volatility, correlation)
    ) / (ds * ds)

    # Total equity vega
    fd_vega = (
        price(stock_price, equity_volatility + dv, fx_volatility, correlation)
        - price(stock_price, equity_volatility - dv, fx_volatility, correlation)
    ) / (2.0 * dv)

    # FX vega
    fd_fx_vega = (
        price(stock_price, equity_volatility, fx_volatility + dv, correlation)
        - price(stock_price, equity_volatility, fx_volatility - dv, correlation)
    ) / (2.0 * dv)

    # Correlation sensitivity
    fd_corr_delta = (
        price(stock_price, equity_volatility, fx_volatility, correlation + dc)
        - price(stock_price, equity_volatility, fx_volatility, correlation - dc)
    ) / (2.0 * dc)

    analytical = (
        option.delta(*args),
        option.gamma(*args),
        option.vega(*args),
        option.fx_vega(*args),
        option.correlation_delta(*args),
    )

    numerical = (
        fd_delta,
        fd_gamma,
        fd_vega,
        fd_fx_vega,
        fd_corr_delta,
    )

    test_cases.header("GREEK", "ANALYTICAL", "FINITE DIFF", "ERROR")

    names = (
        "DELTA",
        "GAMMA",
        "EQUITY VEGA",
        "FX VEGA",
        "CORR DELTA",
    )

    for name, a, n in zip(names, analytical, numerical):

        error = a - n

        test_cases.print(name, a, n, error)

        assert math.isclose(
            a,
            n,
            rel_tol=1.0e-4,
            abs_tol=1.0e-5,
        )


########################################################################################


def test_quanto_fx_rate():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    domestic_curve = FlatDiscountCurve(value_dt, 0.05)
    foreign_curve = FlatDiscountCurve(value_dt, 0.03)
    dividend_curve = FlatDiscountCurve(value_dt, 0.02)

    option_1 = EquityQuantoVanillaOption(
        expiry_dt,
        125.0,
        OptionTypes.EUROPEAN_CALL,
        num_options=100.0,
        quanto_fx_rate=1.0,
    )

    option_2 = EquityQuantoVanillaOption(
        expiry_dt,
        125.0,
        OptionTypes.EUROPEAN_CALL,
        num_options=100.0,
        quanto_fx_rate=1.5,
    )

    args = (
        value_dt,
        130.0,
        domestic_curve,
        foreign_curve,
        dividend_curve,
        BlackScholes(0.25),
        0.12,
        0.50,
    )

    value_1 = option_1.value(*args)
    value_2 = option_2.value(*args)

    assert math.isclose(
        value_2,
        1.5 * value_1,
        rel_tol=1.0e-10,
    )

    test_cases.header("FX RATE 1", "FX RATE 1.5")
    test_cases.print(value_1, value_2)


########################################################################################


test_equity_quanto_vanilla_option()
test_zero_correlation()
test_put_call_parity()
test_quanto_greeks()
test_quanto_fx_rate()

test_cases.compare_test_cases()
