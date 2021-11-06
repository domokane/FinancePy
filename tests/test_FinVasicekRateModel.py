###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.vasicek_mc import zero_price, zero_price_mc


def test_FinModelRatesVasicek():
    r0 = 0.05
    a = 0.10
    b = 0.05
    sigma = 0.05

    num_paths = 1000
    dt = 0.02
    seed = 1968

    t = 0.0
    p_MC = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
    p_MC2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
    p = zero_price(r0, a, b, sigma, t)
    assert round(p, 4) == 1.0
    assert round(p_MC, 4) == 1.0
    assert round(p_MC2, 4) == 1.0

    t = 5.0
    p_MC = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
    p_MC2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
    p = zero_price(r0, a, b, sigma, t)
    assert round(p, 4) == 0.8077
    assert round(p_MC, 4) == 0.8042
    assert round(p_MC2, 4) == 0.8058

    t = 7.5
    p_MC = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
    p_MC2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
    p = zero_price(r0, a, b, sigma, t)
    assert round(p, 4) == 0.7626
    assert round(p_MC, 4) == 0.7567
    assert round(p_MC2, 4) == 0.7620

    t = 10.0
    p_MC = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
    p_MC2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
    p = zero_price(r0, a, b, sigma, t)
    assert round(p, 4) == 0.7483
    assert round(p_MC, 4) == 0.7392
    assert round(p_MC2, 4) == 0.7471
