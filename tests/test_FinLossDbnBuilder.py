###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.gauss_copula_onefactor import loss_dbn_hetero_adj_binomial
from financepy.models.gauss_copula_onefactor import loss_dbn_recursion_gcd
import numpy as np


def test_FinLossDbnBuilder():
    num_steps = 25

    num_credits = 125
    defaultProb = 0.30
    loss_ratio = np.ones(num_credits)
    lossUnits = np.ones(num_credits)

    beta_results = [
        (0.0, [0.0, 0.0, 0.0, 0.0]),
        (0.1, [0.0, 0.0, 0.0, 0.0]),
        (0.2, [0.0, 0.0, 0.0002, 0.0007]),
        (0.3, [0.0024, 0.0137, 0.0452, 0.112]),
        (0.4, [0.1385, 0.4565, 0.9554, 1.6197]),
        (0.5, [1.8407, 3.7058, 5.4378, 7.0117]),
        (0.6, [10.886, 13.7233, 15.1383, 15.91]),
        (0.7, [39.9844, 31.0789, 27.4473, 25.1626]),
        (0.8, [110.1321, 48.1134, 34.0898, 31.7425])
    ]

    for beta, results in beta_results:

        default_probs = np.ones(num_credits) * defaultProb
        beta_vector = np.ones(num_credits) * beta

        dbn1 = loss_dbn_recursion_gcd(num_credits,
                                      default_probs,
                                      lossUnits,
                                      beta_vector,
                                      num_steps)
        assert [round(x * 1000, 4) for x in dbn1[:4]] == results

        dbn2 = loss_dbn_hetero_adj_binomial(num_credits,
                                            default_probs,
                                            loss_ratio,
                                            beta_vector,
                                            num_steps)
        assert [round(x * 1000, 4) for x in dbn2[:4]] == results
