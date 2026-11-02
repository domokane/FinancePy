# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.models.gauss_copula_lhp import portfolio_cdf_lhp
from financepy.models.gauss_copula_lhp import prob_l_greater_than_k

num_credits = 100
default_prob = 0.03
recovery = 0.40
qvector = np.full(num_credits, 1.0 - default_prob)
recovery_rates = np.full(num_credits, recovery)


########################################################################################


def test_prob_l_greater_than_k_is_one_minus_the_loss_cdf():

    for beta in [0.1, 0.3, 0.5, 0.7, 0.9]:
        for k in np.linspace(0.001, 0.59, 25):
            cdf = portfolio_cdf_lhp(k, num_credits, qvector, recovery_rates, beta)
            prob = prob_l_greater_than_k(k, default_prob, recovery, beta)
            assert prob == pytest.approx(1.0 - cdf, abs=1e-12), (beta, k)


def test_prob_l_greater_than_k_limits():

    beta = 0.5

    # a loss is certain to exceed zero when defaults are possible
    assert prob_l_greater_than_k(0.0, default_prob, recovery, beta) == 1.0

    # the loss can never exceed the maximum loss 1 - r
    assert prob_l_greater_than_k(1.0 - recovery, default_prob, recovery, beta) == 0.0
    assert prob_l_greater_than_k(0.9, default_prob, recovery, beta) == 0.0

    # no defaults means no loss
    assert prob_l_greater_than_k(0.05, 0.0, recovery, beta) == 0.0

    # decreasing in k and close to its limits at the ends
    probs = [prob_l_greater_than_k(k, default_prob, recovery, beta) for k in np.linspace(1e-6, 0.6 - 1e-6, 50)]
    assert all(np.diff(probs) <= 1e-12)
    assert probs[0] == pytest.approx(1.0, abs=1e-3)
    assert probs[-1] == pytest.approx(0.0, abs=1e-3)


def test_prob_l_greater_than_k_known_value():

    # Reference case from the loss distribution: k = 5% of notional
    cdf = portfolio_cdf_lhp(0.05, num_credits, qvector, recovery_rates, 0.5)
    prob = prob_l_greater_than_k(0.05, default_prob, recovery, 0.5)
    assert prob == pytest.approx(1.0 - cdf, abs=1e-12)
    assert prob == pytest.approx(0.0859, abs=2e-4)
