# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import numpy as np

import add_fp_to_path

from financepy.models.cir_montecarlo import zero_price_mc, zero_price
from financepy.models.cir_montecarlo import CIRNumericalScheme
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_model_rates_cir():

    r0 = 0.05
    a = 0.20
    b = 0.05
    sigma = 0.20
    t = 5.0

    num_paths = 2000
    dt = 0.05
    seed = 1968

    test_cases.header(
        "MATURITY",
        "TIME",
        "FORMULA",
        "EULER",
        "LOGNORM",
        "MILSTEIN",
        "KJ",
        "EXACT",
    )

    for t in np.linspace(0, 10, 21):

        start = time.time()
        p = zero_price(r0, a, b, sigma, t)
        p_mc1 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.EULER.value,
        )
        p_mc2 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.LOGNORMAL.value,
        )
        p_mc3 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.MILSTEIN.value,
        )
        p_mc4 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.KAHLJACKEL.value,
        )
        p_mc5 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.EXACT.value,
        )
        end = time.time()
        elapsed = end - start
        test_cases.print(t, elapsed, p, p_mc1, p_mc2, p_mc3, p_mc4, p_mc5)


########################################################################################

test_fin_model_rates_cir()
test_cases.compare_test_cases()
