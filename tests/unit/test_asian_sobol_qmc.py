"""Scrambled Sobol quasi-Monte Carlo valuation of Asian options (issue #19)."""

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_asian_option_mc import asian_value_mc_sobol
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize(
    "option_type", (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT)
)
def test_sobol_matches_zero_volatility_discrete_average(option_type):
    t_avg, t_exp = 0.25, 1.0
    strike, spot, rate, yield_rate = 102.0, 100.0, 0.04, 0.01
    n = 12
    dt = (t_exp - t_avg) / n
    expected_average = np.mean([
        spot * np.exp((rate - yield_rate) * (t_avg + (j + 1) * dt))
        for j in range(n)
    ])
    expected_payoff = (
        max(expected_average - strike, 0.0)
        if option_type is OptionTypes.EUROPEAN_CALL
        else max(strike - expected_average, 0.0)
    )
    expected_price = np.exp(-rate * t_exp) * expected_payoff

    price, error = asian_value_mc_sobol(
        t_avg, t_exp, strike, 16, option_type.value,
        spot, rate, yield_rate, 0.0, 1024, 42,
    )
    assert price == pytest.approx(expected_price, abs=1e-11)
    assert error >= 0.0
    assert error < 1e-11


def test_sobol_replicated_seed_is_deterministic_and_has_finite_error():
    args = (0.25, 1.0, 100.0, 24, OptionTypes.EUROPEAN_CALL.value,
            100.0, 0.04, 0.01, 0.2, 4096, 123)
    value1, error1 = asian_value_mc_sobol(*args)
    value2, error2 = asian_value_mc_sobol(*args)
    assert (value1, error1) == (value2, error2)
    assert np.isfinite(value1) and value1 >= 0.0
    assert np.isfinite(error1) and error1 >= 0.0


def test_sobol_accrued_average_zero_volatility():
    t_avg, t_exp = -0.5, 0.5
    strike, accrued, spot = 100.0, 110.0, 95.0
    rate, yield_rate = 0.05, 0.02
    n = 6
    dt = t_exp / n
    remaining_avg = np.mean([
        spot * np.exp((rate - yield_rate) * (j + 1) * dt) for j in range(n)
    ])
    total_avg = 0.5 * accrued + 0.5 * remaining_avg
    expected = np.exp(-rate * t_exp) * max(total_avg - strike, 0.0)

    value, err = asian_value_mc_sobol(
        t_avg, t_exp, strike, 12, OptionTypes.EUROPEAN_CALL.value,
        spot, rate, yield_rate, 0.0, 1024, 42, accrued_average=accrued,
    )
    assert value == pytest.approx(expected, abs=1e-11)
    assert err < 1e-11


@pytest.mark.parametrize("invalid_paths", (31, 1000, 33, 0))
def test_sobol_requires_power_of_two_pairs_per_scramble(invalid_paths):
    with pytest.raises(FinError, match="paths|pairs"):
        asian_value_mc_sobol(
            0.0, 1.0, 100.0, 12, OptionTypes.EUROPEAN_CALL.value,
            100.0, 0.03, 0.01, 0.2, invalid_paths, 1,
        )


def test_sobol_requires_accrued_average_when_started():
    with pytest.raises(FinError, match="accrued|Accrued"):
        asian_value_mc_sobol(
            -0.2, 0.5, 100.0, 12, OptionTypes.EUROPEAN_CALL.value,
            100.0, 0.03, 0.01, 0.2, 1024, 1,
        )


def test_equity_asian_option_sobol_public_api():
    value_dt = Date(1, 1, 2025)
    option = EquityAsianOption(
        Date(1, 4, 2025), Date(1, 1, 2026), 100.0,
        OptionTypes.EUROPEAN_CALL, 24,
    )
    discount = FlatDiscountCurve(value_dt, 0.04)
    dividend = FlatDiscountCurve(value_dt, 0.01)
    model = BlackScholes(0.2)

    result = option.value_mc_sobol(
        value_dt, 100.0, discount, dividend, model,
        num_paths=4096, seed=123, num_scrambles=8,
    )
    assert np.isfinite(result.value)
    assert result.value > 0.0
    assert np.isfinite(result.std_err)
    assert result.std_err >= 0.0
    assert result == option.value_mc_sobol(
        value_dt, 100.0, discount, dividend, model, 4096, 123,
    )
