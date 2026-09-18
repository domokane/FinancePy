########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
from numba import njit, prange


@njit(fastmath=True, parallel=True, cache=True)
def p_double_touch_bb_parallel(
    s0: float,
    L: float,
    U: float,
    mu: float,
    sigma: float,
    t_exp: float,
    steps: int,
    n_paths: int,
    seed: int,
) -> float:
    """
    Estimate P(touch before t_exp) for a double barrier using Brownian bridge.
    Parallel across paths; equality to a barrier counts as a hit.
    """
    if steps < 1:
        steps = 1

    logs0 = np.log(s0)
    lnL = np.log(L)
    lnU = np.log(U)

    dt = t_exp / steps
    nudt = (mu - 0.5 * sigma * sigma) * dt
    sigsdt = sigma * np.sqrt(dt)
    sig2_dt_inv = 1.0 / (sigma * sigma * dt)

    hits = 0

    # Each path uses its own seed offset for reproducibility
    for p in prange(n_paths):
        # Simple per-path reseeding; good enough in practice
        np.random.seed(seed + 1315423911 ^ (p * 2654435761))
        x = logs0
        hit = False

        # immediate hit at t=0?
        if x <= lnL or x >= lnU:
            hits += 1
            continue

        for _ in range(steps):
            z = np.random.randn()
            x_new = x + nudt + sigsdt * z

            # endpoint breaches (include equality as hit)
            if (x <= lnL) or (x >= lnU) or (x_new <= lnL) or (x_new >= lnU):
                hit = True
                break

            # Brownian-bridge hit probability between x and x_new
            # Lower barrier (a = lnL), valid when both endpoints strictly above lnL
            p_lower = 0.0
            if (x > lnL) and (x_new > lnL):
                t1 = (x - lnL) * (x_new - lnL)
                # guard tiny negatives due to rounding
                expo = -2.0 * t1 * sig2_dt_inv
                if expo < 0.0:
                    p_lower = np.exp(expo)
                else:
                    p_lower = 1.0  # if expo >= 0 due to numerical edge, force hit

            # Upper barrier (b = lnU), valid when both endpoints strictly below lnU
            p_upper = 0.0
            if (x < lnU) and (x_new < lnU):
                t2 = (lnU - x) * (lnU - x_new)
                expo = -2.0 * t2 * sig2_dt_inv
                if expo < 0.0:
                    p_upper = np.exp(expo)
                else:
                    p_upper = 1.0

            # Combined probability of hitting either barrier in (t, t+dt)
            # Assuming independence conditional on endpoints (standard approximation)
            p_any = 1.0 - (1.0 - p_lower) * (1.0 - p_upper)
            u = np.random.rand()
            if u < p_any:
                hit = True
                break

            x = x_new

        if hit:
            hits += 1

    return hits / n_paths


########################################################################################


@njit(fastmath=True, cache=True)
def barrier_pay_at_expiry_double_hit(s, k1, k2):
    """Pay $1 if the stock crosses the barrier H from above. PV payment."""
    num_paths, num_time_steps = s.shape
    hits = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            x = s[ip][it]
            if not (k1 < x < k2):
                hit_flag = 1
                break

        hits = hits + hit_flag

    p_hit = hits / num_paths
    return p_hit


########################################################################################


@njit(fastmath=True, cache=True)
def fast_double_no_touch_pricer(s0, L, U, K, t_exp, opt_type, r_d, r_f, sigma):

    df_d = np.exp(-r_d * t_exp)

    # Immediate touch
    if s0 <= L or s0 >= U:
        if opt_type == 2:  # DNT
            return 0.0
        else:  # DOT
            return K * df_d

    # Precompute constants
    Z = np.log(U / L)
    b = r_d - r_f
    sig2 = sigma * sigma
    term1 = 2.0 * b / sig2 - 1.0
    alpha = -0.5 * term1
    beta = -0.25 * term1 * term1 - 2.0 * r_d / sig2

    logSL = np.log(s0 / L)
    logSU = np.log(s0 / U)

    # Heuristic n_max from damping bound
    # exp(-0.5*sig2*(n*pi/Z)^2 * t) <= eps  => n >= ...
    eps = 1e-14
    base = (Z / (np.pi * sigma * np.sqrt(2.0 * t_exp))) * np.sqrt(np.log(1.0 / eps))
    n_max = int(base) + 5
    if n_max < 50:
        n_max = 50
    if n_max > 2000:
        n_max = 2000

    # (s0/L)^alpha and (s0/U)^alpha
    SL_alpha = np.exp(alpha * logSL)
    SU_alpha = np.exp(alpha * logSU)

    # Series sum
    c = 0.0

    for i in range(1, n_max + 1):
        # m = i * pi / Z
        m = (i * np.pi) / Z

        # (-1)^i without pow
        alt = -1.0 if (i & 1) else 1.0

        denom = (alpha * alpha) + (m * m)
        term_i_b = SL_alpha - alt * SU_alpha

        # sin(m * log(s0/L))
        s = np.sin(m * logSL)

        # exp damp: exp(-0.5*(m^2 - beta)*sig2*t)
        damp = np.exp(-0.5 * (m * m - beta) * sig2 * t_exp)

        # assemble term: 2*pi*i*K/Z^2 * (term_i_b/denom) * s * damp
        term = (2.0 * np.pi * i * K / (Z * Z)) * (term_i_b / denom) * s * damp

        c += term

    return c  # DNT price (PV)
