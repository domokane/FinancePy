"""
Exact Heston MC Greeks via AADC (Automatic Adjoint Differentiation).

This example records FinancePy's Heston EulerLog MC scheme on an AADC tape
and replays it for all Monte Carlo paths in a single batch call, producing
exact pathwise Greeks w.r.t. all 6 inputs (S0, v0, kappa, theta, sigma, rho)
at negligible extra cost compared to pricing alone.

Results (200K paths, 50 steps, Heston ATM call):

  FinancePy MC (1 price):               ~0.9s
  FinancePy MC + FD (13 runs, 6 Greeks): ~12s
  AADC batch (price + 6 Greeks):         ~0.35s    (34x faster)

All Greeks verified against finite differences on the same paths (ratio ~1.00).

Key AADC technique:
  - Random normals are marked with mark_as_input_no_diff() (no gradient needed)
  - An array of values (one per scenario) is passed for each random input
  - A single aadc.evaluate() call replays the tape for all scenarios in parallel
  - Model parameters are marked with mark_as_input() (gradient computed)

Requirements:
  pip install aadc    # free evaluation, matlogica.com/aadc
  pip install financepy
"""
import numpy as np
import time
import sys

sys.setrecursionlimit(10000)

from financepy.models.heston import Heston, HestonNumericalScheme
from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date

try:
    import aadc
    from aadc import idouble
except ImportError:
    raise ImportError(
        "This example requires the AADC package.\n"
        "Install: pip install aadc\n"
        "Free evaluation version, no licence key required.\n"
        "See: https://matlogica.com/aadc"
    )


# == Record Heston EulerLog on AADC tape =======================================

def record_heston_mc(S0, r, q, v0, kappa, theta, sigma, rho, T, K, n_steps):
    """Record one Heston EulerLog path + European call payoff on AADC tape.

    Returns (functions, param_args, z_args, payoff_result, n_z).
    """
    dt = T / n_steps
    n_z = 2 * n_steps  # z1 (vol) and z2 (stock) per step

    funcs = aadc.Functions()
    funcs.start_recording()

    # Model parameters — marked as differentiable inputs
    S_id = idouble(S0);    S_arg = S_id.mark_as_input()
    v0_id = idouble(v0);   v0_arg = v0_id.mark_as_input()
    kappa_id = idouble(kappa); kappa_arg = kappa_id.mark_as_input()
    theta_id = idouble(theta); theta_arg = theta_id.mark_as_input()
    sigma_id = idouble(sigma); sigma_arg = sigma_id.mark_as_input()
    rho_id = idouble(rho);     rho_arg = rho_id.mark_as_input()

    # Random normals — no gradient, but values change per scenario
    z_arr = aadc.array(np.random.randn(1, n_z))
    z_args = z_arr.mark_as_input_no_diff()

    # EulerLog scheme (matches FinancePy's get_paths with EULERLOG)
    rhohat = (idouble(1.0) - rho_id * rho_id) ** 0.5
    sdt = np.sqrt(dt)
    x = aadc.math.log(S_id)
    v = v0_id

    for t in range(n_steps):
        z_v = z_arr[0][t] * idouble(sdt)
        z_s = rho_id * z_arr[0][t] * idouble(sdt) + \
              rhohat * z_arr[0][n_steps + t] * idouble(sdt)

        vplus = aadc.iif(v > 0, v, idouble(0.0))
        rtvplus = vplus ** 0.5

        x = x + (idouble(r - q) - vplus * idouble(0.5)) * idouble(dt) + \
            rtvplus * z_s

        sigma2 = sigma_id * sigma_id
        v = v + kappa_id * (theta_id - vplus) * idouble(dt) + \
            sigma_id * rtvplus * z_v + \
            sigma2 * (z_arr[0][t] * z_arr[0][t] * idouble(dt) - idouble(dt)) * idouble(0.25)

    S_final = aadc.math.exp(x)
    payoff = aadc.iif(S_final > K, S_final - idouble(K), idouble(0.0)) * \
             idouble(np.exp(-r * T))
    payoff_res = payoff.mark_as_output()
    funcs.stop_recording()

    param_args = [S_arg, v0_arg, kappa_arg, theta_arg, sigma_arg, rho_arg]
    return funcs, param_args, z_args, payoff_res, n_z


# == Main ======================================================================

if __name__ == '__main__':
    # ── Setup ─────────────────────────────────────────────────────────────
    S0, r, q = 100.0, 0.05, 0.02
    v0, kappa, theta, sigma, rho = 0.04, 2.0, 0.04, 0.3, -0.7
    T, K = 1.0, 100.0
    n_steps = 50
    n_paths = 200000

    pnames = ['S0', 'v0', 'kappa', 'theta', 'sigma', 'rho']
    param_vals = [S0, v0, kappa, theta, sigma, rho]

    # ── FinancePy reference ───────────────────────────────────────────────
    print("=" * 60)
    print("Heston MC Greeks: AADC vs FinancePy")
    print("=" * 60)

    class SimpleOption:
        def __init__(self, expiry_dt, strike, opt_type):
            self.expiry_dt = expiry_dt
            self.strike_price = strike
            self.opt_type = opt_type

    value_dt = Date(1, 1, 2024)
    expiry_dt = Date(1, 1, 2025)
    option = SimpleOption(expiry_dt, K, OptionTypes.EUROPEAN_CALL)
    fp_model = Heston(v0=v0, kappa=kappa, theta=theta, sigma=sigma, rho=rho)

    t0 = time.time()
    fp_price = fp_model.value_mc(value_dt, option, S0, r, q,
                                  n_paths, n_steps, 42,
                                  HestonNumericalScheme.EULERLOG)
    t_fp = time.time() - t0
    print(f"FinancePy MC ({n_paths} paths): {fp_price:.4f}  ({t_fp:.2f}s)")

    fp_lewis = fp_model.value_lewis(value_dt, option, S0, r, q)
    print(f"FinancePy Lewis (exact):   {fp_lewis:.4f}")

    # ── AADC tape ─────────────────────────────────────────────────────────
    print(f"\nRecording Heston EulerLog on AADC tape ({n_steps} steps)...")
    np.random.seed(42)
    all_Z = np.random.randn(n_paths, 2 * n_steps)

    funcs, param_args, z_args, payoff_res, n_z = \
        record_heston_mc(S0, r, q, v0, kappa, theta, sigma, rho, T, K, n_steps)

    # Build inputs: scalar params + array randoms
    inputs = {param_args[j]: param_vals[j] for j in range(6)}
    for zi in range(n_z):
        inputs[z_args[0][zi]] = all_Z[:, zi].copy()

    workers = aadc.ThreadPool(4)

    # ── Batch evaluation ──────────────────────────────────────────────────
    t0 = time.time()
    result = aadc.evaluate(funcs, {payoff_res: param_args}, inputs, workers)
    t_aadc = time.time() - t0

    price = np.average(result[0][payoff_res])
    grads = {p: np.average(result[1][payoff_res][param_args[j]])
             for j, p in enumerate(pnames)}

    print(f"\nAADC batch ({n_paths} paths): {price:.4f}  ({t_aadc:.2f}s)")
    print(f"\nGreeks (exact pathwise, single reverse pass):")
    for p in pnames:
        print(f"  d/d({p:>6s}) = {grads[p]:+.6f}")

    # ── FD verification (same paths, same tape) ──────────────────────────
    print(f"\nFD verification (same paths):")
    h_vals = {'S0': 1.0, 'v0': 0.005, 'kappa': 0.1,
              'theta': 0.005, 'sigma': 0.03, 'rho': 0.01}
    for j, p in enumerate(pnames):
        h = h_vals[p]
        pv_up = list(param_vals); pv_up[j] += h
        pv_dn = list(param_vals); pv_dn[j] -= h
        inp_up = {param_args[i]: pv_up[i] for i in range(6)}
        inp_dn = {param_args[i]: pv_dn[i] for i in range(6)}
        for zi in range(n_z):
            inp_up[z_args[0][zi]] = all_Z[:, zi].copy()
            inp_dn[z_args[0][zi]] = all_Z[:, zi].copy()
        r_up = aadc.evaluate(funcs, {payoff_res: []}, inp_up, workers)
        r_dn = aadc.evaluate(funcs, {payoff_res: []}, inp_dn, workers)
        fd = (np.average(r_up[0][payoff_res]) - np.average(r_dn[0][payoff_res])) / (2*h)
        ratio = grads[p] / fd if abs(fd) > 1e-12 else float('nan')
        print(f"  d/d({p:>6s}): AD={grads[p]:+10.6f}  FD={fd:+10.6f}  ratio={ratio:.4f}")

    # ── Benchmark ─────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Benchmark")
    print("=" * 60)

    # AADC: price + 6 Greeks
    times_aadc = []
    for _ in range(3):
        t0 = time.time()
        aadc.evaluate(funcs, {payoff_res: param_args}, inputs, workers)
        times_aadc.append(time.time() - t0)
    t_aadc_med = sorted(times_aadc)[1]

    # FinancePy: 1 price
    times_fp = []
    for _ in range(3):
        t0 = time.time()
        fp_model.value_mc(value_dt, option, S0, r, q,
                          n_paths, n_steps, 42, HestonNumericalScheme.EULERLOG)
        times_fp.append(time.time() - t0)
    t_fp_med = sorted(times_fp)[1]

    # FinancePy FD: price + 6 Greeks (13 MC runs)
    t0 = time.time()
    fp_model.value_mc(value_dt, option, S0, r, q,
                      n_paths, n_steps, 42, HestonNumericalScheme.EULERLOG)
    for p in pnames:
        for sign in [+1, -1]:
            pv = dict(zip(pnames, param_vals))
            pv[p] += sign * h_vals[p]
            m = Heston(v0=pv['v0'], kappa=pv['kappa'], theta=pv['theta'],
                       sigma=pv['sigma'], rho=pv['rho'])
            m.value_mc(value_dt, option, pv['S0'], r, q,
                       n_paths, n_steps, 42, HestonNumericalScheme.EULERLOG)
    t_fp_fd = time.time() - t0

    print(f"  AADC (price + 6 Greeks): {t_aadc_med:.2f}s")
    print(f"  FinancePy (1 price):     {t_fp_med:.2f}s")
    print(f"  FinancePy + FD (13 runs):{t_fp_fd:.1f}s")
    print(f"  Speedup vs FD:           {t_fp_fd/t_aadc_med:.0f}x")
    print(f"  Speedup vs 1 price:      {t_fp_med/t_aadc_med:.1f}x")
