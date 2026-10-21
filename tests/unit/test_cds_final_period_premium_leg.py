# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS

########################################################################################


def flat_hazard_curve(value_dt: Date, hazard_rate: float, recovery_rate: float):
    libor_curve = FlatDiscountCurve(
        value_dt, 0.03, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    curve = CDSCurve(value_dt, [], libor_curve, recovery_rate)
    curve.set_times(np.array([0.0, 50.0]))
    curve.set_qs(np.array([1.0, np.exp(-hazard_rate * 50.0)]))
    return curve


def test_cds_premium_leg_in_final_coupon_period():
    """The first premium uses the accrual factor of the current period. In the
    final coupon period there is no following period, so the premium leg used
    to read past the end of the accrual factors and return a negative clean
    RPV01 and par spread."""

    hazard_rate = 0.02
    recovery_rate = 0.4
    value_dt = Date(20, 10, 2026)
    curve = flat_hazard_curve(value_dt, hazard_rate, recovery_rate)

    cds = CDS(Date(20, 9, 2026), Date(20, 12, 2026), 0.01)
    dirty_rpv01, clean_rpv01 = cds.rpv01(value_dt, curve)

    # One premium period of 92 days (ACT/360) remains: the risky PV01 cannot
    # exceed its accrual factor and the clean value cannot be negative
    accrual_factor = cds.accrual_factors[0]
    assert 0.0 < clean_rpv01 < dirty_rpv01 < accrual_factor

    # On a flat hazard curve, the par spread of the stub is close to the par
    # spread of a 5y contract, (1 - R) * h up to accrual conventions
    par_spread = cds.par_spread(value_dt, curve, recovery_rate)
    long_cds = CDS(value_dt, Date(20, 12, 2031), 0.01)
    long_par_spread = long_cds.par_spread(value_dt, curve, recovery_rate)
    assert abs(par_spread / long_par_spread - 1.0) < 0.03
    assert abs(par_spread / ((1.0 - recovery_rate) * hazard_rate) - 1.0) < 0.03

