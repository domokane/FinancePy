# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 21:48:46 2026

@author: Dominic
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ...utils.date import Date


@dataclass(frozen=True)
class Cashflow:
    payment_dt: Date
    amount: Optional[float] = None  # None = not yet known (floating)

    accrual_start_dt: Optional[Date] = None
    accrual_end_dt: Optional[Date] = None
    year_frac: Optional[float] = None

    is_principal: bool = False
