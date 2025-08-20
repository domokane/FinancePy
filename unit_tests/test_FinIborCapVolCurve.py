########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
import numpy as np


def test_FinCapVolCurve():
    # Reproduces example in Table 32.1 of Hull Book
    value_dt = Date(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    num_periods = 10
    capletDt = value_dt

    capVolDates.append(value_dt)

    for i in range(0, num_periods):
        capletDt = capletDt.add_tenor(capletVolTenor)
        capVolDates.append(capletDt)

    capVolatilities = [
        0.0,
        15.50,
        18.25,
        17.91,
        17.74,
        17.27,
        16.79,
        16.30,
        16.01,
        15.76,
        15.54,
    ]
    capVolatilities = np.array(capVolatilities) / 100.0

    day_count_type = DayCountTypes.ACT_ACT_ISDA
    vol_curve = IborCapVolCurve(
        value_dt, capVolDates, capVolatilities, day_count_type
    )

    dt = Date(1, 1, 2020)
    capFloorVol = vol_curve.cap_vol(dt)
    capFloorLetVol = vol_curve.caplet_vol(dt)
    assert round(capFloorVol * 100.0, 2) == 15.5
    assert round(capFloorLetVol * 100.0, 2) == 15.5

    dt = Date(1, 1, 2022)
    capFloorVol = vol_curve.cap_vol(dt)
    capFloorLetVol = vol_curve.caplet_vol(dt)
    assert round(capFloorVol * 100.0, 2) == 18.25
    assert round(capFloorLetVol * 100.0, 2) == 20.64

    dt = Date(1, 1, 2024)
    capFloorVol = vol_curve.cap_vol(dt)
    capFloorLetVol = vol_curve.caplet_vol(dt)
    assert round(capFloorVol * 100.0, 2) == 17.740
    assert round(capFloorLetVol * 100.0, 2) == 17.22

    dt = Date(1, 1, 2026)
    capFloorVol = vol_curve.cap_vol(dt)
    capFloorLetVol = vol_curve.caplet_vol(dt)
    assert round(capFloorVol * 100.0, 2) == 16.790
    assert round(capFloorLetVol * 100.0, 2) == 14.15

    dt = Date(1, 1, 2028)
    capFloorVol = vol_curve.cap_vol(dt)
    capFloorLetVol = vol_curve.caplet_vol(dt)
    assert round(capFloorVol * 100.0, 2) == 16.010
    assert round(capFloorLetVol * 100.0, 2) == 13.81
