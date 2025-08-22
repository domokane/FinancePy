# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.models.cir_montecarlo import CIRNumericalScheme
from financepy.models.cir_montecarlo import zero_price_mc, zero_price
import numpy as np

r0 = 0.05
a = 0.20
b = 0.05
sigma = 0.20
t = 5.0

num_paths = 2000
dt = 0.05
seed = 1968

########################################################################################


def test_model_cir():

    p = zero_price(r0, a, b, sigma, t)
    p_mc1 = zero_price_mc(
        r0, a, b, sigma, t, dt, num_paths, seed, CIRNumericalScheme.EULER.value
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
        r0, a, b, sigma, t, dt, num_paths, seed, CIRNumericalScheme.EXACT.value
    )

    assert round(p, 4) == 0.7935
    assert round(p_mc1, 4) == 0.7913
    assert round(p_mc2, 4) == 0.7912
    assert round(p_mc3, 4) == 0.7914
    assert round(p_mc4, 4) == 0.7909
    assert round(p_mc5, 4) == 0.7922
