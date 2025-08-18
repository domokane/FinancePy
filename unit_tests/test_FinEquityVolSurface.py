###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

from financepy.models.volatility_fns import VolFuncTypes
from financepy.utils.date import Date
from financepy.market.volatility.equity_vol_surface import EquityVolSurface
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.day_count import DayCountTypes


def test_equity_vol_surface():
    value_dt = Date(11, 1, 2021)

    stock_price = 3800.0  # Check

    expiry_dts = [
        Date(11, 2, 2021),
        Date(11, 3, 2021),
        Date(11, 4, 2021),
        Date(11, 7, 2021),
        Date(11, 10, 2021),
        Date(11, 1, 2022),
        Date(11, 1, 2023),
    ]

    strikes = np.array([3037, 3418, 3608, 3703, 3798, 3893, 3988, 4178, 4557])

    vol_surface = [
        [42.94, 31.30, 25.88, 22.94, 19.72, 16.90, 15.31, 17.54, 25.67],
        [37.01, 28.25, 24.19, 21.93, 19.57, 17.45, 15.89, 15.34, 21.15],
        [34.68, 27.38, 23.82, 21.85, 19.83, 17.98, 16.52, 15.31, 18.94],
        [31.41, 26.25, 23.51, 22.05, 20.61, 19.25, 18.03, 16.01, 15.90],
        [29.91, 25.58, 23.21, 22.01, 20.83, 19.70, 18.62, 16.63, 14.94],
        [29.26, 25.24, 23.03, 21.91, 20.81, 19.73, 18.69, 16.76, 14.63],
        [27.59, 24.33, 22.72, 21.93, 21.17, 20.43, 19.71, 18.36, 16.26],
    ]

    vol_surface = np.array(vol_surface)
    vol_surface = vol_surface / 100.0

    r = 0.020  # USD
    discount_curve = DiscountCurveFlat(value_dt, r)
    # discount_curve = DiscountCurveFlat(value_dt, r, dc_type=DayCountTypes.SIMPLE)

    q = 0.010  # USD
    dividend_curve = DiscountCurveFlat(value_dt, q)
    # dividend_curve = DiscountCurveFlat(value_dt, q, dc_type=DayCountTypes.SIMPLE)

    vol_function_type = VolFuncTypes.SVI

    equitySurface = EquityVolSurface(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        expiry_dts,
        strikes,
        vol_surface,
        vol_function_type,
    )

    expiry_dt = expiry_dts[0]
    delta = 0.10
    vol = equitySurface.vol_from_delta_date(delta, expiry_dt)
    assert round(vol[0], 4) == 0.1544
    assert round(vol[1], 4) == 4032.9156

    expiry_dt = expiry_dts[1]
    delta = 0.20
    vol = equitySurface.vol_from_delta_date(delta, expiry_dt)
    assert round(vol[0], 4) == 0.1555
    assert round(vol[1], 4) == 4019.3793

    expiry_dt = expiry_dts[6]
    delta = 0.90
    vol = equitySurface.vol_from_delta_date(delta, expiry_dt)
    assert (
        round(vol[0], 4) == 0.3498
    )  # 0.353 # 0.3498 VP TODO: had to rebase, not sure why. Investigate more. Interestingly the original numbers pass on github so restored them
    assert (
        round(vol[1], 4) == 2199.6665
    )  # 2190.7766 # 2199.6665 VP TODO: had to rebase, not sure why. Investigate more. Interestingly the original numbers pass on github so restored them
