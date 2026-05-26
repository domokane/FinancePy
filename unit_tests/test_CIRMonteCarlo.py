# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.models.cir_montecarlo import (
    CIRNumericalScheme,
    CIRMonteCarlo,
    satisfies_feller,
    meanr,
    variancer,
    zero_price,
    draw,
    rate_path_mc,
    zero_price_mc,
)

R0 = 0.03
A = 0.50
B = 0.04
SIGMA = 0.10
T = 2.0
DT = 1.0 / 252.0
SEED = 42


def test_cir_monte_carlo_stores_inputs():

    model = CIRMonteCarlo(A, B, SIGMA)

    assert model._a == A
    assert model._b == B
    assert model._sigma == SIGMA


def test_cir_repr():

    model = CIRMonteCarlo(A, B, SIGMA)

    s = repr(model)

    assert "CIRMonteCarlo" in s
    assert "Sigma" in s
    assert "a" in s
    assert "b" in s


def test_satisfies_feller():

    assert satisfies_feller(0.50, 0.04, 0.10)
    assert not satisfies_feller(0.10, 0.01, 0.20)


def test_meanr_at_zero_time():

    assert np.isclose(meanr(R0, A, B, 0.0), R0, atol=1e-14)


def test_variancer_at_zero_time():

    assert np.isclose(variancer(R0, A, B, SIGMA, 0.0), 0.0, atol=1e-14)


def test_meanr_long_time_reverts_to_b():

    value = meanr(R0, A, B, 100.0)

    assert np.isclose(value, B, atol=1e-10)


def test_variancer_is_non_negative():

    value = variancer(R0, A, B, SIGMA, T)

    assert value >= 0.0


def test_zero_price_at_zero_time():

    assert zero_price(R0, A, B, SIGMA, 0.0) == 1.0


def test_zero_price_is_between_zero_and_one():

    value = zero_price(R0, A, B, SIGMA, T)

    assert 0.0 < value <= 1.0


def test_zero_price_zero_volatility_matches_deterministic_integral():

    sigma = 0.0

    value = zero_price(R0, A, B, sigma, T)

    integral_mean = B * T + (R0 - B) * (1.0 - np.exp(-A * T)) / A
    expected = np.exp(-integral_mean)

    assert np.isclose(value, expected, atol=1e-14)


def test_draw_is_non_negative():

    for _ in range(100):
        r = draw(R0, A, B, SIGMA, DT)
        assert r >= 0.0


@pytest.mark.parametrize(
    "scheme",
    [
        CIRNumericalScheme.EULER,
        CIRNumericalScheme.LOGNORMAL,
        CIRNumericalScheme.MILSTEIN,
        CIRNumericalScheme.KAHLJACKEL,
        CIRNumericalScheme.EXACT,
    ],
)
def test_rate_path_has_expected_length(scheme):

    t = 1.0
    dt = 0.30

    path = rate_path_mc(
        R0,
        A,
        B,
        SIGMA,
        t,
        dt,
        SEED,
        scheme.value,
    )

    expected_steps = int(np.ceil(t / dt))

    assert len(path) == expected_steps + 1
    assert path[0] == R0


@pytest.mark.parametrize(
    "scheme",
    [
        CIRNumericalScheme.LOGNORMAL,
        CIRNumericalScheme.EXACT,
    ],
)
def test_positive_schemes_produce_non_negative_paths(scheme):

    path = rate_path_mc(
        R0,
        A,
        B,
        SIGMA,
        T,
        DT,
        SEED,
        scheme.value,
    )

    assert np.all(path >= 0.0)


@pytest.mark.parametrize(
    "scheme",
    [
        CIRNumericalScheme.EULER,
        CIRNumericalScheme.LOGNORMAL,
        CIRNumericalScheme.MILSTEIN,
        CIRNumericalScheme.KAHLJACKEL,
        CIRNumericalScheme.EXACT,
    ],
)
def test_rate_path_is_reproducible_with_same_seed(scheme):

    path1 = rate_path_mc(
        R0,
        A,
        B,
        SIGMA,
        T,
        DT,
        SEED,
        scheme.value,
    )

    path2 = rate_path_mc(
        R0,
        A,
        B,
        SIGMA,
        T,
        DT,
        SEED,
        scheme.value,
    )

    assert np.allclose(path1, path2)


@pytest.mark.parametrize(
    "scheme",
    [
        CIRNumericalScheme.EULER,
        CIRNumericalScheme.LOGNORMAL,
        CIRNumericalScheme.MILSTEIN,
        CIRNumericalScheme.KAHLJACKEL,
        CIRNumericalScheme.EXACT,
    ],
)
def test_zero_price_mc_at_zero_time(scheme):

    value = zero_price_mc(
        R0,
        A,
        B,
        SIGMA,
        0.0,
        DT,
        100,
        SEED,
        scheme.value,
    )

    assert value == 1.0


@pytest.mark.parametrize(
    "scheme",
    [
        CIRNumericalScheme.LOGNORMAL,
        CIRNumericalScheme.EXACT,
    ],
)
def test_zero_price_mc_reasonable_against_analytic(scheme):

    analytic = zero_price(R0, A, B, SIGMA, 1.0)

    mc = zero_price_mc(
        R0,
        A,
        B,
        SIGMA,
        1.0,
        1.0 / 52.0,
        20000,
        SEED,
        scheme.value,
    )

    assert np.isclose(mc, analytic, atol=5.0e-3)


def test_exact_scheme_terminal_mean_matches_analytic_mean():

    num_paths = 20000
    terminal_rates = np.empty(num_paths)

    for i in range(num_paths):
        path = rate_path_mc(
            R0,
            A,
            B,
            SIGMA,
            1.0,
            1.0 / 12.0,
            SEED + i,
            CIRNumericalScheme.EXACT.value,
        )
        terminal_rates[i] = path[-1]

    expected = meanr(R0, A, B, 1.0)

    assert np.isclose(np.mean(terminal_rates), expected, atol=2.0e-3)


def test_exact_scheme_terminal_variance_matches_analytic_variance():

    num_paths = 20000
    terminal_rates = np.empty(num_paths)

    for i in range(num_paths):
        path = rate_path_mc(
            R0,
            A,
            B,
            SIGMA,
            1.0,
            1.0 / 12.0,
            SEED + i,
            CIRNumericalScheme.EXACT.value,
        )
        terminal_rates[i] = path[-1]

    expected = variancer(R0, A, B, SIGMA, 1.0)

    assert np.isclose(np.var(terminal_rates), expected, atol=2.0e-4)


@pytest.mark.parametrize(
    "bad_args",
    [
        (-0.01, A, B, SIGMA, T, DT),
        (R0, 0.0, B, SIGMA, T, DT),
        (R0, A, -0.01, SIGMA, T, DT),
        (R0, A, B, -0.01, T, DT),
        (R0, A, B, SIGMA, -1.0, DT),
        (R0, A, B, SIGMA, T, 0.0),
    ],
)
def test_rate_path_rejects_bad_inputs(bad_args):

    r0, a, b, sigma, t, dt = bad_args

    with pytest.raises(ValueError):
        rate_path_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            SEED,
            CIRNumericalScheme.EXACT.value,
        )


def test_zero_price_mc_rejects_non_positive_num_paths():

    with pytest.raises(ValueError):
        zero_price_mc(
            R0,
            A,
            B,
            SIGMA,
            T,
            DT,
            0,
            SEED,
            CIRNumericalScheme.EXACT.value,
        )
