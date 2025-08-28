"""
Fast polynomial fit with optional Numba acceleration.
Set FINANCEPY_USE_NUMBA=1 to enable Numba (import & JIT deferred until first use).
"""

from __future__ import annotations
import os
from typing import Callable, Optional
import numpy as np

# --------------------------------------------------------------------
# Config: default OFF for fastest package import
# Turn on by: export FINANCEPY_USE_NUMBA=1
# --------------------------------------------------------------------
_USE_NUMBA_DEFAULT = os.getenv("FINANCEPY_USE_NUMBA", "0") == "1"

# Internal handle for (maybe) jitted callables
_fit_poly_impl: Optional[Callable[[np.ndarray, np.ndarray, int], np.ndarray]] = None
_eval_poly_impl: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None

# ----------------------- Pure NumPy implementations ------------------


def _coeff_mat_numpy(x: np.ndarray, deg: int) -> np.ndarray:
    # fast and stable: column order matches coefficients of increasing power
    # (we’ll reverse later to match your API)
    # vander with increasing=True gives [1, x, x^2, ...]
    return np.vander(x, N=deg + 1, increasing=True)


def _fit_poly_numpy(x: np.ndarray, y: np.ndarray, deg: int) -> np.ndarray:
    a = _coeff_mat_numpy(x, deg)
    # numpy lstsq: returns coeffs for increasing powers; reverse to highest-first
    p_inc, *_ = np.linalg.lstsq(a, y, rcond=None)
    return p_inc[::-1]


def _eval_polynomial_numpy(p: np.ndarray, x: np.ndarray) -> np.ndarray:
    # Horner with highest-order first (your API)
    result = np.zeros_like(x, dtype=np.result_type(p.dtype, x.dtype))
    for coeff in p:
        result = x * result + coeff
    return result


# ----------------------- Lazy Numba wiring ---------------------------


def _enable_numba() -> None:
    """Import numba and build jitted implementations on demand."""
    global _fit_poly_impl, _eval_poly_impl

    import numba as nb  # heavy import; done only once, on first enable

    # Numba-compatible pieces
    @nb.njit("f8[:,:](f8[:], i8)", cache=True)
    def _coeff_mat_numba(x, deg):
        mat = np.zeros((x.shape[0], deg + 1))
        c = np.ones_like(x)
        mat[:, 0] = c
        if deg >= 1:
            mat[:, 1] = x
        if deg > 1:
            for n in range(2, deg + 1):
                mat[:, n] = x**n
        return mat

    # NOTE: np.linalg.lstsq support in Numba varies by version.
    # If your Numba doesn’t support it natively, consider QR or normal equations.
    @nb.njit("f8[:](f8[:,:], f8[:])", cache=True)
    def _fit_x_numba(a, b):
        # Least squares: fall back to normal equations if needed:
        # x = solve(a.T @ a, a.T @ b)
        # Replace with lstsq if your Numba version supports it well.
        ata = a.T @ a
        atb = a.T @ b
        x = np.linalg.solve(ata, atb)
        return x

    @nb.njit("f8[:](f8[:], f8[:], i8)", cache=True)
    def _fit_poly_numba(x, y, deg):
        a = _coeff_mat_numba(x, deg)
        p_inc = _fit_x_numba(a, y)  # increasing order
        return p_inc[::-1]  # highest-first

    @nb.njit(fastmath=True, cache=True)
    def _eval_polynomial_numba(p, x):
        result = np.zeros_like(x)
        for coeff in p:
            result = x * result + coeff
        return result

    _fit_poly_impl = _fit_poly_numba
    _eval_poly_impl = _eval_polynomial_numba


def _ensure_impls():
    """Choose NumPy or Numba impls; import Numba lazily if requested."""
    global _fit_poly_impl, _eval_poly_impl
    if _fit_poly_impl is not None and _eval_poly_impl is not None:
        return
    if _USE_NUMBA_DEFAULT:
        try:
            _enable_numba()
            return
        except Exception:
            # Fall back silently to NumPy if Numba unavailable or unsupported op
            pass
    _fit_poly_impl = _fit_poly_numpy
    _eval_poly_impl = _eval_poly_numpy


# ----------------------- Public API ---------------------------------


def fit_poly(x: np.ndarray, y: np.ndarray, deg: int) -> np.ndarray:
    """
    Fit polynomial of degree `deg` to points (x, y).
    Returns coefficients with highest order first (same as original API).
    """
    _ensure_impls()
    return _fit_poly_impl(x.astype(np.float64), y.astype(np.float64), int(deg))


def eval_polynomial(p: np.ndarray, x: np.ndarray) -> np.ndarray:
    """
    Evaluate polynomial with coefficients `p` (highest-first) at points `x`.
    """
    _ensure_impls()
    p = np.asarray(p, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    return _eval_poly_impl(p, x)
