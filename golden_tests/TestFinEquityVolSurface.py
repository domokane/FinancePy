# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
import matplotlib.pyplot as plt
from FinTestCases import FinTestCases, global_test_case_mode
from financepy.models.volatility_fns import VolFuncTypes
from financepy.utils.date import Date
from financepy.market.volatility.equity_vol_surface import EquityVolSurface
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
import numpy as np


test_cases = FinTestCases(__file__, global_test_case_mode)


PLOT_GRAPHS = False

# TODO: ADD LOGGING TO TEST CASES

########################################################################################


def test_equity_vol_surface(verbose_calibration):

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

    q = 0.010  # USD
    dividend_curve = DiscountCurveFlat(value_dt, q)

    vol_function_type = VolFuncTypes.SVI

    equity_surface = EquityVolSurface(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        expiry_dts,
        strikes,
        vol_surface,
        vol_function_type,
    )

    #    tol = 1e-4
    #    equity_surface.check_calibration(False, tol)

    if PLOT_GRAPHS:

        equity_surface.plot_vol_curves()

        plt.figure()

        mins = strikes[0] * 0.5
        maxs = strikes[-1] * 1.5

        dbns = equity_surface.implied_dbns(mins, maxs, 1000)

        for i in range(0, len(dbns)):
            expiry_dt_str = str(equity_surface._expiry_dts[i])
            plt.plot(dbns[i]._x, dbns[i]._densitydx, label=expiry_dt_str)
            plt.title(vol_function_type)
            plt.legend()
            print("SUM:", dbns[i].sum())

    deltas = np.linspace(0.10, 0.90, 9)

    test_cases.header("EXPIRY", "DELTA", "VOL", "STRIKE")
    for expiry_dt in expiry_dts:
        for delta in deltas:
            vol = equity_surface.vol_from_delta_date(delta, expiry_dt)
            test_cases.print(expiry_dt, delta, vol[0], vol[1])


########################################################################################

########################################################################################

if __name__ == "__main__":

    start = time.time()

    verbose_calibration = False

    test_equity_vol_surface(verbose_calibration)

    end = time.time()

    elapsed = end - start
    #    print("Elapsed Time:", elapsed)
    test_cases.compare_test_cases()
