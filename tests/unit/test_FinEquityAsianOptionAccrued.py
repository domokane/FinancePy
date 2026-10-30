# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
"""Tests for the Asian option during-averaging branch.

The pricing tests in test_FinEquityAsianOption.py all pass
accrued_avg = None, so only the pre-averaging path is covered there. This
module covers a value date that falls inside the averaging window, where the
accrued average has to be folded into the pricing.
"""

import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.utils.helpers import option_years
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.products.equity.equity_asian_option import AsianOptionValuationTypes
from financepy.utils.global_types import OptionTypes
from financepy.utils.error import FinError

STOCK_PRICE = 100.0
VOLATILITY = 0.20
INTEREST_RATE = 0.05
DIVIDEND_YIELD = 0.02
STRIKE = 100.0

START_AVERAGING_DT = Date(1, 1, 2015)
EXPIRY_DT = Date(1, 1, 2016)

METHODS = [
    AsianOptionValuationTypes.KEMNA_VORST,
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    AsianOptionValuationTypes.CURRAN,
]

# Accrued average observed part-way through the averaging period.
ACCRUED_AVERAGE = 103.0
VALUE_DT = Date(1, 7, 2015)


def make_option(opt_type=OptionTypes.EUROPEAN_CALL, strike=STRIKE, num_obs=100):
    return EquityAsianOption(
        START_AVERAGING_DT,
        EXPIRY_DT,
        strike,
        opt_type,
        num_obs,
    )


def make_market(value_dt, volatility=VOLATILITY):
    discount_curve = FlatDiscountCurve(value_dt, INTEREST_RATE)
    dividend_curve = FlatDiscountCurve(value_dt, DIVIDEND_YIELD)
    model = BlackScholes(volatility)
    return discount_curve, dividend_curve, model


########################################################################################
# Helpers reproducing the accrued-average convention independently of the library
########################################################################################


def averaging_times(value_dt):
    """Return (t_exp, t_avg, tau) as the library computes them."""
    t_exp = option_years(value_dt, EXPIRY_DT)
    t_avg = option_years(value_dt, START_AVERAGING_DT, fail=False)
    return t_exp, t_avg, t_exp - t_avg


def mc_asian_value_during_averaging(
    value_dt,
    accrued_average,
    opt_type,
    num_obs,
    num_paths,
    seed,
    volatility=VOLATILITY,
):
    """Monte Carlo value with the payoff defined here, not by the library.

    The average over the whole averaging window is

        A_final = (accrued_average * |t_avg| + A_remaining * t_exp) / tau

    where A_remaining is the mean of the observations still to come and
    tau = t_exp - t_avg is the length of the whole averaging window. This is
    the convention implied by the strike adjustment and multiplier the library
    applies, and it is checked against the analytic pricers below rather than
    assumed correct.
    """
    t_exp, t_avg, tau = averaging_times(value_dt)

    num_obs_remaining = max(1, int(t_exp * num_obs + 0.5))
    dt = t_exp / num_obs_remaining

    rng = np.random.default_rng(seed)
    z = rng.standard_normal((num_paths, num_obs_remaining))

    # Brownian motion accumulated over the observation grid: each observation
    # must be correlated with its neighbours, otherwise the average has far
    # too little variance and the option is badly underpriced.
    log_increments = (
        (INTEREST_RATE - DIVIDEND_YIELD - volatility**2 / 2.0) * dt
        + volatility * np.sqrt(dt) * z
    )
    log_s = np.log(STOCK_PRICE) + np.cumsum(log_increments, axis=1)
    avg_remaining = np.exp(log_s).mean(axis=1)

    avg_final = (
        accrued_average * abs(t_avg) + avg_remaining * t_exp
    ) / tau

    if opt_type == OptionTypes.EUROPEAN_CALL:
        payoff = np.maximum(avg_final - STRIKE, 0.0)
    else:
        payoff = np.maximum(STRIKE - avg_final, 0.0)

    discount_factor = np.exp(-INTEREST_RATE * t_exp)
    values = discount_factor * payoff
    return values.mean(), values.std() / np.sqrt(num_paths)


########################################################################################
# 1. Monte Carlo cross-check -- the strongest test, since the payoff is
#    defined here and shares no code with the analytic pricers.
########################################################################################


ARITHMETIC_METHODS = [
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    AsianOptionValuationTypes.CURRAN,
]


@pytest.mark.parametrize("method", ARITHMETIC_METHODS)
def test_asian_value_during_averaging_matches_monte_carlo(method):
    """TW and Curran approximate the arithmetic average, so both must sit on
    the Monte Carlo value of the arithmetic payoff."""
    discount_curve, dividend_curve, model = make_market(VALUE_DT)

    analytic = make_option().value(
        VALUE_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        method,
        ACCRUED_AVERAGE,
    )

    mc, stderr = mc_asian_value_during_averaging(
        VALUE_DT,
        ACCRUED_AVERAGE,
        OptionTypes.EUROPEAN_CALL,
        num_obs=100,
        num_paths=200_000,
        seed=1976,
    )

    # Three standard errors of a 200k-path run is ~1.5e-2 here, comfortably
    # tighter than the spread between the three analytic approximations.
    assert analytic == pytest.approx(mc, abs=3.0 * stderr + 5.0e-3)


def test_asian_geometric_is_lower_bound_during_averaging():
    """Kemna-Vorst prices the geometric average, which is a strict lower bound
    for the arithmetic average, so it must come in under the Monte Carlo value."""
    discount_curve, dividend_curve, model = make_market(VALUE_DT)

    geometric = make_option().value(
        VALUE_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.KEMNA_VORST,
        ACCRUED_AVERAGE,
    )

    mc, stderr = mc_asian_value_during_averaging(
        VALUE_DT,
        ACCRUED_AVERAGE,
        OptionTypes.EUROPEAN_CALL,
        num_obs=100,
        num_paths=200_000,
        seed=1976,
    )

    assert geometric < mc - 3.0 * stderr


########################################################################################
# 2. Put-call parity. The reference E[A_final] is a first-order expression, so
#    it carries a small truncation error; what matters is that the residual is
#    independent of the accrued average. If the accrued term were mishandled
#    the residual would move with it by an O(1) amount.
########################################################################################


def _parity_residual(value_dt, accrued_average, method, num_obs=100):
    discount_curve, dividend_curve, model = make_market(value_dt)
    t_exp, t_avg, tau = averaging_times(value_dt)

    call = make_option(OptionTypes.EUROPEAN_CALL, num_obs=num_obs).value(
        value_dt, STOCK_PRICE, discount_curve, dividend_curve, model,
        method, accrued_average,
    )
    put = make_option(OptionTypes.EUROPEAN_PUT, num_obs=num_obs).value(
        value_dt, STOCK_PRICE, discount_curve, dividend_curve, model,
        method, accrued_average,
    )

    # Mean observation time of the remaining fixings, as the library discretises it.
    num_obs_remaining = max(1, int(t_exp * num_obs + 0.5))
    dt = t_exp / num_obs_remaining
    mean_time = dt * (num_obs_remaining + 1) / 2.0

    expected_avg = (
        accrued_average * abs(t_avg)
        + STOCK_PRICE * np.exp((INTEREST_RATE - DIVIDEND_YIELD) * mean_time) * t_exp
    ) / tau

    parity = call - put
    forward = np.exp(-INTEREST_RATE * t_exp) * (expected_avg - STRIKE)
    return parity - forward


@pytest.mark.parametrize("method", ARITHMETIC_METHODS)
def test_asian_put_call_parity_during_averaging(method):
    for value_dt in [Date(1, 3, 2015), Date(1, 7, 2015), Date(1, 10, 2015)]:
        for accrued_average in [90.0, 103.0, 115.0]:
            residual = _parity_residual(value_dt, accrued_average, method)
            assert abs(residual) < 5.0e-3


def test_asian_parity_residual_is_independent_of_accrued_average():
    """The parity residual must not move when the accrued average changes.

    A mishandled accrued term shifts parity by an amount proportional to the
    accrued average, so a flat residual across a wide range of accrued averages
    is a direct check on that term.
    """
    residuals = [
        _parity_residual(VALUE_DT, accrued, AsianOptionValuationTypes.CURRAN)
        for accrued in [80.0, 90.0, 100.0, 105.0]
    ]

    spread = max(residuals) - min(residuals)
    assert spread < 1.0e-9


########################################################################################
# 3. Limit check: as the value date moves back to the start of the averaging
#    window the accrued average must stop having any effect.
########################################################################################


@pytest.mark.parametrize("method", METHODS)
def test_asian_accrued_average_has_no_effect_before_averaging(method):
    discount_curve, dividend_curve, model = make_market(START_AVERAGING_DT)

    with_accrued = make_option().value(
        START_AVERAGING_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        method,
        103.0,
    )

    without_accrued = make_option().value(
        START_AVERAGING_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        method,
        None,
    )

    assert with_accrued == pytest.approx(without_accrued, abs=1.0e-12)


########################################################################################
# 4. Monotonicity in the accrued average. A higher accrued average pulls the
#    final average up, so calls gain and puts lose.
########################################################################################


@pytest.mark.parametrize("method", METHODS)
def test_asian_value_is_monotone_in_accrued_average(method):
    discount_curve, dividend_curve, model = make_market(VALUE_DT)
    accrued_averages = [80.0, 90.0, 100.0, 110.0, 120.0]

    calls = [
        make_option(OptionTypes.EUROPEAN_CALL).value(
            VALUE_DT, STOCK_PRICE, discount_curve, dividend_curve, model,
            method, accrued,
        )
        for accrued in accrued_averages
    ]

    puts = [
        make_option(OptionTypes.EUROPEAN_PUT).value(
            VALUE_DT, STOCK_PRICE, discount_curve, dividend_curve, model,
            method, accrued,
        )
        for accrued in accrued_averages
    ]

    assert all(calls[i] < calls[i + 1] for i in range(len(calls) - 1))
    assert all(puts[i] > puts[i + 1] for i in range(len(puts) - 1))


########################################################################################
# 5. Error handling: inside the averaging window the accrued average is
#    required, and it has to be positive.
########################################################################################


@pytest.mark.parametrize("method", METHODS)
def test_asian_requires_accrued_average_inside_averaging(method):
    discount_curve, dividend_curve, model = make_market(VALUE_DT)

    with pytest.raises(Exception):
        make_option().value(
            VALUE_DT,
            STOCK_PRICE,
            discount_curve,
            dividend_curve,
            model,
            method,
            None,
        )


@pytest.mark.parametrize("method", METHODS)
def test_asian_rejects_non_positive_accrued_average(method):
    discount_curve, dividend_curve, model = make_market(VALUE_DT)

    with pytest.raises(FinError):
        make_option().value(
            VALUE_DT,
            STOCK_PRICE,
            discount_curve,
            dividend_curve,
            model,
            method,
            0.0,
        )


########################################################################################
# 6. The rescaling onto the remaining window must be algebraically exact.
#
#    This is the check that makes the negative-strike case a library bug rather
#    than a matter of convention: whatever accrued_average means, the two ways
#    of writing the payoff have to agree, and they do -- except where a negative
#    adjusted strike makes Kemna-Vorst and Turnbull-Wakeman fail outright.
########################################################################################


@pytest.mark.parametrize(
    "value_dt",
    [Date(1, 8, 2015), Date(1, 11, 2015), Date(1, 12, 2015), Date(28, 12, 2015)],
)
@pytest.mark.parametrize("avg_remaining", [90.0, 100.0, 110.0])
def test_asian_rescaling_is_algebraically_exact(value_dt, avg_remaining):
    t_exp, t_avg, tau = averaging_times(value_dt)

    k_adj = (STRIKE * tau + ACCRUED_AVERAGE * t_avg) / t_exp
    multiplier = t_exp / tau

    avg_final = (
        ACCRUED_AVERAGE * abs(t_avg) + avg_remaining * t_exp
    ) / tau

    payoff_direct = avg_final - STRIKE
    payoff_rescaled = multiplier * (avg_remaining - k_adj)

    assert payoff_direct == pytest.approx(payoff_rescaled, abs=1.0e-9)


########################################################################################
# 7. The reported bug: with an accrued average high enough and a value date late
#    enough in the averaging window, the adjusted strike turns negative and
#    Kemna-Vorst / Turnbull-Wakeman raise ZeroDivisionError while Curran, using
#    the same adjustment, returns the value the identity above predicts.
########################################################################################


NEGATIVE_STRIKE_VALUE_DT = Date(1, 12, 2015)
NEGATIVE_STRIKE_ACCRUED = 115.0


def test_asian_handles_negative_adjusted_strike():
    """An accrued average high enough late in the window drives the adjusted
    strike negative. Kemna-Vorst and Turnbull-Wakeman used to raise
    ZeroDivisionError here because log(s0 / k) went undefined.

    The option is guaranteed to finish in the money in this regime, so the
    value must come out finite and positive, and it must agree with Curran
    which already handled it.
    """
    discount_curve, dividend_curve, model = make_market(NEGATIVE_STRIKE_VALUE_DT)

    t_exp, t_avg, tau = averaging_times(NEGATIVE_STRIKE_VALUE_DT)
    k_adj = (STRIKE * tau + NEGATIVE_STRIKE_ACCRUED * t_avg) / t_exp
    assert k_adj < 0.0  # guards the premise of this test

    values = []
    for method in METHODS:
        value = make_option().value(
            NEGATIVE_STRIKE_VALUE_DT,
            STOCK_PRICE,
            discount_curve,
            dividend_curve,
            model,
            method,
            NEGATIVE_STRIKE_ACCRUED,
        )
        assert np.isfinite(value)
        assert value > 0.0
        values.append(value)

    # Kemna-Vorst prices the geometric average so it sits a little below the
    # two arithmetic pricers, but the spread has to stay small.
    assert max(values) - min(values) < 5.0e-3


def test_asian_put_is_worthless_on_negative_adjusted_strike():
    """Same regime, put side: the adjusted strike is negative and the average
    cannot fall below it, so the put has no value at all."""
    discount_curve, dividend_curve, model = make_market(NEGATIVE_STRIKE_VALUE_DT)

    for method in METHODS:
        value = make_option(OptionTypes.EUROPEAN_PUT).value(
            NEGATIVE_STRIKE_VALUE_DT,
            STOCK_PRICE,
            discount_curve,
            dividend_curve,
            model,
            method,
            NEGATIVE_STRIKE_ACCRUED,
        )
        assert value == pytest.approx(0.0, abs=1.0e-12)


def test_asian_arithmetic_pricers_agree_on_negative_adjusted_strike():
    """Turnbull-Wakeman and Curran must land on the same number.

    Both approximate the arithmetic average, and in this degenerate regime the
    payoff is a plain forward, so there is nothing left for them to disagree
    about. This is the tightest available pin on the fix.
    """
    discount_curve, dividend_curve, model = make_market(NEGATIVE_STRIKE_VALUE_DT)

    turnbull_wakeman = make_option().value(
        NEGATIVE_STRIKE_VALUE_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        NEGATIVE_STRIKE_ACCRUED,
    )

    curran = make_option().value(
        NEGATIVE_STRIKE_VALUE_DT,
        STOCK_PRICE,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.CURRAN,
        NEGATIVE_STRIKE_ACCRUED,
    )

    assert turnbull_wakeman == pytest.approx(curran, abs=1.0e-9)
