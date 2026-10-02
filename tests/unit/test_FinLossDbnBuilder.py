# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.models.gauss_copula_onefactor import loss_dbn_hetero_adj_binomial
from financepy.models.gauss_copula_onefactor import loss_dbn_recursion_gcd
from financepy.models.loss_dbn_builder import portfolio_gcd
from financepy.models.gauss_copula_onefactor import tranche_surv_prob_recursion
from financepy.utils.math import pair_gcd

########################################################################################


def test_fin_loss_dbn_builder():

    num_steps = 25

    num_credits = 125
    default_prob = 0.30
    loss_ratio = np.ones(num_credits)
    loss_units = np.ones(num_credits)

    beta_results = [
        (0.0, [0.0, 0.0, 0.0, 0.0]),
        (0.1, [0.0, 0.0, 0.0, 0.0]),
        (0.2, [0.0, 0.0, 0.0002, 0.0007]),
        (0.3, [0.0024, 0.0137, 0.0452, 0.112]),
        (0.4, [0.1385, 0.4565, 0.9554, 1.6197]),
        (0.5, [1.8407, 3.7058, 5.4378, 7.0117]),
        (0.6, [10.886, 13.7233, 15.1383, 15.91]),
        (0.7, [39.9844, 31.0789, 27.4473, 25.1626]),
        (0.8, [110.1321, 48.1134, 34.0898, 31.7425]),
    ]

    for beta, results in beta_results:

        default_probs = np.ones(num_credits) * default_prob
        beta_vector = np.ones(num_credits) * beta

        dbn1 = loss_dbn_recursion_gcd(
            num_credits, default_probs, loss_units, beta_vector, num_steps
        )
        assert [round(x * 1000, 4) for x in dbn1[:4]] == results

        dbn2 = loss_dbn_hetero_adj_binomial(
            num_credits, default_probs, loss_ratio, beta_vector, num_steps
        )
        assert [round(x * 1000, 4) for x in dbn2[:4]] == results


def test_pair_gcd_uses_integer_remainders_and_handles_zero():
    assert pair_gcd(12.0, 8.0) == 4.0
    assert pair_gcd(0.0, 8.0) == 8.0
    assert pair_gcd(8.0, 0.0) == 8.0
    assert pair_gcd(0.0, 0.0) == 0.0
    assert pair_gcd(-12.0, 8.0) == 4.0


@pytest.mark.parametrize(
    "losses, expected", [([0.12, 0.08], 0.04), ([0.0, 0.2, 0.3], 0.1)]
)
def test_portfolio_gcd_preserves_heterogeneous_loss_units(losses, expected):
    assert portfolio_gcd(np.array(losses)) == pytest.approx(expected)


def test_tranche_survival_handles_zero_loss_credits():
    result = tranche_surv_prob_recursion(
        0.0,
        0.25,
        2,
        np.array([0.9, 0.9]),
        np.array([1.0, 0.5]),
        np.array([0.0, 0.0]),
        100,
    )

    assert result == pytest.approx(0.9, abs=1e-3)


def test_tranche_survival_is_one_when_all_recoveries_are_full():
    result = tranche_surv_prob_recursion(
        0.0,
        0.25,
        2,
        np.array([0.9, 0.9]),
        np.array([1.0, 1.0]),
        np.array([0.0, 0.0]),
        100,
    )

    assert result == 1.0
