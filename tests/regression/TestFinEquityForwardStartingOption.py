
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

import numpy as np

from financepy.products.equity.equity_forward_start_option import (
    EquityForwardStartOption,
)
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_european_analytic import european_value
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date
from financepy.utils.helpers import option_years

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_forward_start_option():

    value_dt = Date(1, 1, 2015)
    start_dt = Date(1, 1, 2017)
    expiry_dt = Date(1, 1, 2018)

    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    test_cases.header(
        "TYPE",
        "START",
        "EXPIRY",
        "S",
        "VOL",
        "VALUE",
        "EXPECTED",
        "DIFF",
    )

    for opt_type in [
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ]:

        option = EquityForwardStartOption(
            start_dt,
            expiry_dt,
            opt_type,
        )

        value = option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        t = option_years(start_dt, expiry_dt)

        fwd_r = discount_curve.fwd_zero_rate_cc(
            start_dt, expiry_dt
        )

        fwd_q = dividend_curve.fwd_zero_rate_cc(
            start_dt, expiry_dt
        )

        dq_start = dividend_curve.df(start_dt)

        expected = stock_price * dq_start * european_value(
            1.0,
            t,
            1.0,
            fwd_r,
            fwd_q,
            volatility,
            opt_type.value,
        )

        diff = value - expected

        test_cases.print(
            opt_type,
            start_dt,
            expiry_dt,
            stock_price,
            volatility,
            value,
            expected,
            diff,
        )

        assert np.isclose(value, expected, rtol=1e-10, atol=1e-10)


########################################################################################


def test_equity_forward_start_stock_prices():

    value_dt = Date(1, 1, 2015)
    start_dt = Date(1, 1, 2017)
    expiry_dt = Date(1, 1, 2018)

    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    stock_prices = [50.0, 70.0, 85.0, 100.0, 120.0]

    test_cases.header(
        "TYPE",
        "S",
        "VALUE",
        "VALUE/S",
        "SCALE_DIFF",
    )

    for opt_type in [
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ]:

        option = EquityForwardStartOption(
            start_dt,
            expiry_dt,
            opt_type,
        )

        base_value = option.value(
            value_dt,
            100.0,
            discount_curve,
            dividend_curve,
            model,
        )

        for stock_price in stock_prices:

            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            expected = base_value * stock_price / 100.0
            diff = value - expected

            test_cases.print(
                opt_type,
                stock_price,
                value,
                value / stock_price,
                diff,
            )

            assert np.isclose(value, expected, rtol=1e-10, atol=1e-10)


########################################################################################


def test_equity_forward_start_volatilities():

    value_dt = Date(1, 1, 2015)
    start_dt = Date(1, 1, 2017)
    expiry_dt = Date(1, 1, 2018)

    stock_price = 85.0
    interest_rate = 0.035
    dividend_yield = 0.01

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    volatilities = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50]

    test_cases.header(
        "TYPE",
        "VOL",
        "VALUE",
    )

    for opt_type in [
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ]:

        option = EquityForwardStartOption(
            start_dt,
            expiry_dt,
            opt_type,
        )

        previous_value = None

        for volatility in volatilities:

            model = BlackScholes(volatility)

            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            test_cases.print(
                opt_type,
                volatility,
                value,
            )

            assert value >= 0.0

            if previous_value is not None:
                assert value >= previous_value - 1e-10

            previous_value = value


########################################################################################


def test_equity_forward_start_expiries():

    value_dt = Date(1, 1, 2015)

    start_dt = Date(1, 1, 2017)

    expiry_dts = [
        Date(1, 4, 2017),
        Date(1, 7, 2017),
        Date(1, 1, 2018),
        Date(1, 1, 2019),
    ]

    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    test_cases.header(
        "TYPE",
        "START",
        "EXPIRY",
        "VALUE",
    )

    for opt_type in [
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ]:

        for expiry_dt in expiry_dts:

            option = EquityForwardStartOption(
                start_dt,
                expiry_dt,
                opt_type,
            )

            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            test_cases.print(
                opt_type,
                start_dt,
                expiry_dt,
                value,
            )

            assert np.isfinite(value)
            assert value >= 0.0


########################################################################################


def test_equity_forward_start_put_call_parity():

    value_dt = Date(1, 1, 2015)
    start_dt = Date(1, 1, 2017)
    expiry_dt = Date(1, 1, 2018)

    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    call_option = EquityForwardStartOption(
        start_dt,
        expiry_dt,
        OptionTypes.EUROPEAN_CALL,
    )

    put_option = EquityForwardStartOption(
        start_dt,
        expiry_dt,
        OptionTypes.EUROPEAN_PUT,
    )

    call_value = call_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    put_value = put_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    dq_start = dividend_curve.df(start_dt)
    dq_end = dividend_curve.df(expiry_dt)

    df_start = discount_curve.df(start_dt)
    df_end = discount_curve.df(expiry_dt)

    # Forward-start put-call parity:
    #
    # C - P = S0 * [Q(0,T2) - Q(0,T1)*P(0,T2)/P(0,T1)]

    parity_value = stock_price * (
        dq_end - dq_start * df_end / df_start
    )

    diff = call_value - put_value - parity_value

    test_cases.header(
        "CALL",
        "PUT",
        "CALL-PUT",
        "PARITY",
        "DIFF",
    )

    test_cases.print(
        call_value,
        put_value,
        call_value - put_value,
        parity_value,
        diff,
    )

    assert np.isclose(diff, 0.0, atol=1e-10)


########################################################################################


test_equity_forward_start_option()
test_equity_forward_start_stock_prices()
test_equity_forward_start_volatilities()
test_equity_forward_start_expiries()
test_equity_forward_start_put_call_parity()

test_cases.compare_test_cases()
