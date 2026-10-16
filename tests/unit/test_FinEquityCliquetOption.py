# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest
from financepy.utils.error import FinError
from financepy.models.bs_european_analytic import european_value

from financepy.utils.global_types import (
    OptionTypes,
    CliquetTypes,
    PaymentTimingTypes,
)
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_cliquet_option import EquityCliquetOption


@pytest.mark.parametrize("value_dt", [Date(4, 1, 2026), Date(1, 7, 2025)])
def test_cliquet_rejects_dates_before_first_reset(value_dt):
    """A pre-inception stock level must not be treated as an observed reset."""
    curve = FlatDiscountCurve(value_dt, .03)
    option = EquityCliquetOption(Date(5, 1, 2026), Date(4, 1, 2027), OptionTypes.EUROPEAN_CALL, FrequencyTypes.ANNUAL)
    with pytest.raises(FinError, match="before start date"):
        option.value(value_dt, 100.0, curve, curve, BlackScholes(.2))


def test_one_period_cliquet_at_inception_matches_vanilla():
    """At the first reset the observed stock sets the one-period strike."""
    start = Date(5, 1, 2026)
    expiry = Date(4, 1, 2027)
    option = EquityCliquetOption(start, expiry, OptionTypes.EUROPEAN_CALL, FrequencyTypes.ANNUAL)
    curve = FlatDiscountCurve(start, .03)
    dividend = FlatDiscountCurve(start, .01)
    expected = european_value(100.0, (expiry - start) / 365.0, 100.0, .03, .01, .2, OptionTypes.EUROPEAN_CALL.value)
    assert option.value(start, 100.0, curve, dividend, BlackScholes(.2)) == pytest.approx(expected)

########################################################################################


def assert_close(value, expected, tol=2.0e-3):
    assert np.isclose(value, expected, atol=tol), (
        f"value={value:.10f}, expected={expected:.10f}, " f"diff={value - expected:.10f}"
    )


def test_equity_cliquet_option():

    start_dt = Date(1, 1, 2014)
    final_expiry_dt = Date(1, 1, 2017)
    freq_type = FrequencyTypes.QUARTERLY
    opt_type = OptionTypes.EUROPEAN_CALL

    cliquet_option = EquityCliquetOption(start_dt, final_expiry_dt, opt_type, freq_type)

    value_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    v = cliquet_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    # 12 quarterly forward-start options from 1 Jan 2015. The previous expected
    # value of 34.531 included a spurious one-day option (1 Jan to 2 Jan 2014)
    # produced by regenerating the schedule from its adjusted termination date.
    assert_close(v, 34.136)


def test_cliquet_validation_caps_and_floors():
    """Verify input validation for local and global caps/floors."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)

    # local_floor > local_cap
    with pytest.raises(FinError, match="local_floor must be less than or equal to local_cap"):
        EquityCliquetOption(
            start,
            expiry,
            OptionTypes.EUROPEAN_CALL,
            FrequencyTypes.QUARTERLY,
            local_floor=0.05,
            local_cap=0.02,
        )

    # global_floor > global_cap
    with pytest.raises(FinError, match="global_floor must be less than or equal to global_cap"):
        EquityCliquetOption(
            start,
            expiry,
            OptionTypes.EUROPEAN_CALL,
            FrequencyTypes.QUARTERLY,
            payoff_timing_type=PaymentTimingTypes.MATURITY,
            global_floor=0.10,
            global_cap=0.05,
        )

    # global_floor with PERIODIC payment timing
    with pytest.raises(
        FinError,
        match="global_floor and global_cap are only supported when payoff_timing_type is PaymentTimingTypes.MATURITY",
    ):
        EquityCliquetOption(
            start,
            expiry,
            OptionTypes.EUROPEAN_CALL,
            FrequencyTypes.QUARTERLY,
            payoff_timing_type=PaymentTimingTypes.PERIODIC,
            global_floor=0.0,
        )

    # Analytical valuation does not support global floor
    cliquet_global = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        FrequencyTypes.QUARTERLY,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        global_floor=0.0,
    )
    curve = FlatDiscountCurve(start, 0.05)
    model = BlackScholes(0.20)
    with pytest.raises(
        FinError,
        match="Analytical valuation does not support global floor or global cap",
    ):
        cliquet_global.value(start, 100.0, curve, curve, model)


def test_cliquet_local_floor_and_cap_call_analytic_vs_mc():
    """Verify analytical valuation matches Monte Carlo for a Call Cliquet with local floor and cap."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)
    freq = FrequencyTypes.QUARTERLY

    curve = FlatDiscountCurve(start, 0.05)
    dividend = FlatDiscountCurve(start, 0.02)
    model = BlackScholes(0.20)
    stock_price = 100.0

    # Test both PERIODIC and MATURITY
    for timing in (PaymentTimingTypes.PERIODIC, PaymentTimingTypes.MATURITY):
        cliquet = EquityCliquetOption(
            start,
            expiry,
            OptionTypes.EUROPEAN_CALL,
            freq,
            payoff_type=CliquetTypes.RETURN,
            payoff_timing_type=timing,
            local_floor=-0.02,
            local_cap=0.05,
        )

        v_analytic = cliquet.value(start, stock_price, curve, dividend, model)
        mc_res = cliquet.value_mc(
            start,
            stock_price,
            curve,
            dividend,
            model,
            num_paths=200_000,
            seed=1234,
        )

        z_score = abs(v_analytic - mc_res.value) / mc_res.std_err
        assert z_score < 3.0, (
            f"Call timing {timing}: analytic={v_analytic:.6f}, "
            f"mc={mc_res.value:.6f}, z={z_score:.2f}"
        )


def test_cliquet_local_floor_and_cap_put_analytic_vs_mc():
    """Verify analytical valuation matches Monte Carlo for a Put Cliquet with local floor and cap."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)
    freq = FrequencyTypes.QUARTERLY

    curve = FlatDiscountCurve(start, 0.05)
    dividend = FlatDiscountCurve(start, 0.02)
    model = BlackScholes(0.20)
    stock_price = 100.0

    cliquet = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_PUT,
        freq,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        local_floor=-0.01,
        local_cap=0.06,
    )

    v_analytic = cliquet.value(start, stock_price, curve, dividend, model)
    mc_res = cliquet.value_mc(
        start,
        stock_price,
        curve,
        dividend,
        model,
        num_paths=200_000,
        seed=5678,
    )

    z_score = abs(v_analytic - mc_res.value) / mc_res.std_err
    assert z_score < 3.0, (
        f"Put: analytic={v_analytic:.6f}, "
        f"mc={mc_res.value:.6f}, z={z_score:.2f}"
    )


def test_cliquet_globally_floored_locally_capped_mc():
    """Verify Globally Floored Locally Capped (GFLC) Cliquet Monte Carlo pricing."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)
    freq = FrequencyTypes.MONTHLY

    curve = FlatDiscountCurve(start, 0.05)
    dividend = FlatDiscountCurve(start, 0.02)
    model = BlackScholes(0.20)
    stock_price = 100.0
    notional = 1_000_000.0

    # Contract specification from Issue #221:
    # 1Y maturity, monthly resets, local floor -2%, local cap +5%, global floor 0%
    cliquet_base = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        freq,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        notional=notional,
        local_floor=-0.02,
        local_cap=0.05,
        global_floor=0.0,
    )

    res_base = cliquet_base.value_mc(
        start, stock_price, curve, dividend, model, num_paths=100_000, seed=42
    )
    assert res_base.value > 0.0

    # With higher global floor (e.g. 2%), price must be strictly higher
    cliquet_higher_floor = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        freq,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        notional=notional,
        local_floor=-0.02,
        local_cap=0.05,
        global_floor=0.02,
    )
    res_higher_floor = cliquet_higher_floor.value_mc(
        start, stock_price, curve, dividend, model, num_paths=100_000, seed=42
    )
    assert res_higher_floor.value > res_base.value

    # With higher local cap (e.g. 8%), price must be strictly higher
    cliquet_higher_cap = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        freq,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        notional=notional,
        local_floor=-0.02,
        local_cap=0.08,
        global_floor=0.0,
    )
    res_higher_cap = cliquet_higher_cap.value_mc(
        start, stock_price, curve, dividend, model, num_paths=100_000, seed=42
    )
    assert res_higher_cap.value > res_base.value


def test_cliquet_repr_includes_caps_floors():
    """Verify that __repr__ includes local and global caps and floors."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)
    cliquet = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        FrequencyTypes.QUARTERLY,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        local_floor=-0.02,
        local_cap=0.05,
        global_floor=0.0,
        global_cap=0.20,
    )
    rep = repr(cliquet)
    assert "LOCAL_FLOOR" in rep
    assert "LOCAL_CAP" in rep
    assert "GLOBAL_FLOOR" in rep
    assert "GLOBAL_CAP" in rep


def test_one_period_cliquet_with_cap_matches_bull_spread():
    """Verify 1-period capped call cliquet analytically matches an exact textbook bull call spread."""
    start = Date(5, 1, 2026)
    expiry = Date(4, 1, 2027)
    s0 = 100.0
    r = 0.03
    q = 0.01
    vol = 0.20
    t = (expiry - start) / 365.0

    curve = FlatDiscountCurve(start, r)
    dividend = FlatDiscountCurve(start, q)
    model = BlackScholes(vol)

    cliquet = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        FrequencyTypes.ANNUAL,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        local_cap=0.10,
    )
    v_cliquet = cliquet.value(start, s0, curve, dividend, model)

    c_atm = european_value(s0, t, s0, r, q, vol, OptionTypes.EUROPEAN_CALL.value)
    c_otm = european_value(s0, t, 1.10 * s0, r, q, vol, OptionTypes.EUROPEAN_CALL.value)
    expected_spread = (c_atm - c_otm) / s0

    assert v_cliquet == pytest.approx(expected_spread, rel=1e-10)


def test_cliquet_high_cap_converges_to_uncapped():
    """Verify that as local_cap -> inf, capped cliquet value converges to uncapped cliquet value."""
    start = Date(1, 1, 2026)
    expiry = Date(1, 1, 2027)
    curve = FlatDiscountCurve(start, 0.05)
    dividend = FlatDiscountCurve(start, 0.02)
    model = BlackScholes(0.20)
    s0 = 100.0

    c_uncapped = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        FrequencyTypes.QUARTERLY,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
    )
    v_uncapped = c_uncapped.value(start, s0, curve, dividend, model)

    c_high_cap = EquityCliquetOption(
        start,
        expiry,
        OptionTypes.EUROPEAN_CALL,
        FrequencyTypes.QUARTERLY,
        payoff_type=CliquetTypes.RETURN,
        payoff_timing_type=PaymentTimingTypes.MATURITY,
        local_cap=10.0,
    )
    v_high_cap = c_high_cap.value(start, s0, curve, dividend, model)

    assert v_high_cap == pytest.approx(v_uncapped, rel=1e-8)

