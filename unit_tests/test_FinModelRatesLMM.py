# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.models.lmm_mc import lmm_sticky_caplet_pricer
from financepy.models.lmm_mc import lmm_ratchet_caplet_pricer
from financepy.models.lmm_mc import lmm_simulate_fwds_mf
from financepy.models.lmm_mc import lmm_simulate_fwds_1f
from financepy.utils.helpers import check_vector_differences

########################################################################################


def test_hull_book_examples(capsys):
    """Examining examples on page 770 of Hull OFODS
    Last cap product has caplet starting in 10 years so we have to model
    the forward curve from time 0 out to 11 forwards, not 10 forwards.
    We have to model forward rates 0-1, 1-2, 2-3, ..., 10-11"""

    # We go out 11 periods because last caplet resets in 10 years
    num_fwds = 11
    dt = 1.00
    taus = np.array([dt] * num_fwds)
    seed = 438

    r = 0.05127
    fwd0 = np.zeros(num_fwds)
    for i in range(0, num_fwds):
        fwd0[i] = r

    num_paths = 500000
    spread = 0.0025  # basis points

    # HULL TABLE 32.1

    use_sobol = 1
    numeraire_index = 0

    # We need the volatility for the forward rates out to the one starting in
    # 10 years. So we have 11 elements. The one starting today has zero vol.
    num_factors = 1
    gammas1_f_list = [
        0.00,
        0.1550,
        0.2063674,
        0.1720986,
        0.1721993,
        0.1524579,
        0.1414779,
        0.1297711,
        0.1381053,
        0.135955,
        0.1339842,
    ]
    gammas1_f = np.array(gammas1_f_list)

    # One factor model
    fwds1_f = lmm_simulate_fwds_1f(
        num_fwds, num_paths, numeraire_index, fwd0, gammas1_f, taus, use_sobol, seed
    )

    ########################################################################################

    #    LMMPrintForwards(fwds1_f)

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1_f, taus)
        * 100.0
    )

    hull_ratchet_caplets1_f = [
        0.00,
        0.196,
        0.207,
        0.201,
        0.194,
        0.187,
        0.1890,
        0.172,
        0.167,
        0.160,
        0.153,
    ]

    hull_ratchet_caplets1_f = np.array(hull_ratchet_caplets1_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets1_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1_f, taus)
        * 100.0
    )

    hull_sticky_caplets1_f = [
        0.0,
        0.196,
        0.336,
        0.412,
        0.458,
        0.484,
        0.498,
        0.502,
        0.501,
        0.497,
        0.488,
    ]

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets1_f, 1e-2)

    num_factors = 1
    lambdas1_f_list = [
        [
            0.0,
            0.1550,
            0.2064,
            0.1721,
            0.1722,
            0.1525,
            0.1415,
            0.1298,
            0.1381,
            0.1360,
            0.1340,
        ]
    ]
    lambdas1_f = np.array(lambdas1_f_list)

    # One factor model
    fwds_mf = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas1_f,
        taus,
        use_sobol,
        seed,
    )

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds_mf, taus)
        * 100.0
    )

    hull_ratchet_caplets1_f = [
        0.0,
        0.196,
        0.207,
        0.201,
        0.194,
        0.187,
        0.1890,
        0.172,
        0.167,
        0.160,
        0.153,
    ]

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets1_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds_mf, taus)
        * 100.0
    )

    hull_sticky_caplets1_f = [
        0.00,
        0.196,
        0.336,
        0.412,
        0.458,
        0.484,
        0.498,
        0.502,
        0.501,
        0.497,
        0.488,
    ]

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets1_f, 1e-2)

    num_factors = 3
    lambdas3_f_list = [
        [
            0.00,
            0.1365,
            0.1928,
            0.1672,
            0.1698,
            0.1485,
            0.1395,
            0.1261,
            0.1290,
            0.1197,
            0.1097,
        ],
        [
            0.0,
            -0.0662,
            -0.0702,
            -0.0406,
            -0.0206,
            0.00,
            0.0169,
            0.0306,
            0.0470,
            0.0581,
            0.0666,
        ],
        [
            0.0,
            0.0319,
            0.0225,
            0.000,
            -0.0198,
            -0.0347,
            -0.0163,
            0.000,
            0.0151,
            0.0280,
            0.0384,
        ],
    ]
    lambdas3_f = np.array(lambdas3_f_list)

    # Three factor model
    fwds3_f = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas3_f,
        taus,
        use_sobol,
        seed,
    )

    hull_ratchet_caplets3_f = [
        0.00,
        0.194,
        0.207,
        0.205,
        0.198,
        0.193,
        0.189,
        0.180,
        0.174,
        0.168,
        0.162,
    ]

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3_f, taus)
        * 100.0
    )

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets3_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3_f, taus)
        * 100.0
    )

    hull_sticky_caplets3_f = [
        0.00,
        0.195,
        0.336,
        0.418,
        0.472,
        0.506,
        0.524,
        0.533,
        0.537,
        0.537,
        0.534,
    ]

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets3_f, 1e-2)

    captured = capsys.readouterr()

    assert captured.out == ""
    assert captured.err == ""
