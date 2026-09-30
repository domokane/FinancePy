"""
Created on Sun Sep 27 20:31:07 2026

@author: Dominic
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class MCResult:
    """Result from a Monte Carlo valuation."""

    value: float
    std_err: float
