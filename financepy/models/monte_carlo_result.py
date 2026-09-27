# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 20:31:07 2026

@author: Dominic
"""

from dataclasses import dataclass


import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def standard_error(x):
    """Calculate the standard error of the sample mean."""

    n = x.size

    if n < 2:
        return 0.0

    mean = 0.0

    for i in range(n):
        mean += x[i]

    mean /= n

    variance = 0.0

    for i in range(n):
        diff = x[i] - mean
        variance += diff * diff

    variance /= n - 1

    return np.sqrt(variance / n)


@dataclass(frozen=True)
class MonteCarloResult:
    """Result from a Monte Carlo valuation."""

    value: float
    standard_error: float
