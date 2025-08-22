# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
import numpy as np

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.volatility.swaption_vol_surface import SwaptionVolSurface
from financepy.utils.date import Date
from financepy.models.volatility_fns import VolFuncTypes
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)



plot_graphs = False

# TODO: ADD LOGGING TO TEST CASES


########################################################################################


def test_FinSwaptionVolSurface1():

    ###########################################################################

    # https://fr.mathworks.com/help/fininst/pricing-a-swaption-using-the-sabr-model.html

    value_dt = Date(12, 6, 2013)

    # These are 3M, 1Y, 2Y, 3Y, 4Y, 5Y, 7Y, 10Y
    exercise_dts = [
        Date(12, 9, 2013),
        Date(12, 6, 2014),
        Date(12, 6, 2015),
        Date(12, 6, 2016),
        Date(12, 6, 2017),
        Date(12, 6, 2018),
        Date(12, 6, 2020),
        Date(12, 6, 2023),
    ]

    # First dimension is the strike, then the expiry date
    market_volatilities = [
        [57.6, 53.7, 49.4, 45.6, 44.1, 41.1, 35.2, 32.0],
        [46.6, 46.9, 44.8, 41.6, 39.8, 37.4, 33.4, 31.0],
        [35.9, 39.3, 39.6, 37.9, 37.2, 34.7, 30.5, 28.9],
        [34.1, 36.5, 37.8, 36.6, 35.0, 31.9, 28.1, 26.6],
        [41.0, 41.3, 39.5, 37.8, 36.0, 32.6, 29.0, 26.0],
        [45.8, 43.4, 41.9, 39.2, 36.9, 33.2, 29.6, 26.3],
        [50.3, 46.9, 44.0, 40.0, 37.5, 33.8, 30.2, 27.3],
    ]

    market_volatilities = np.array(market_volatilities) / 100.0

    # First dimension is the strike, then the expiry date
    market_strikes = [
        [1.00, 1.25, 1.68, 2.00, 2.26, 2.41, 2.58, 2.62],
        [1.50, 1.75, 2.18, 2.50, 2.76, 2.91, 3.08, 3.12],
        [2.00, 2.25, 2.68, 3.00, 3.26, 3.41, 3.58, 3.62],
        [2.50, 2.75, 3.18, 3.50, 3.76, 3.91, 4.08, 4.12],
        [3.00, 3.25, 3.68, 4.00, 4.26, 4.41, 4.58, 4.62],
        [3.50, 3.75, 4.18, 4.50, 4.76, 4.91, 5.08, 5.12],
        [4.00, 4.25, 4.68, 5.00, 5.26, 5.41, 5.58, 5.62],
    ]

    market_strikes = np.array(market_strikes) / 100.0

    fwd_swap_rates = market_strikes[3]
    atm_vols = market_volatilities[3]

    rfr_rate = 0.020  # USD
    discount_curve = DiscountCurveFlat(value_dt, rfr_rate)

    div_rate = 0.010  # USD
    dividend_curve = DiscountCurveFlat(value_dt, div_rate)

    vol_function_type = VolFuncTypes.SABR_BETA_HALF

    swaption_surface = SwaptionVolSurface(
        value_dt,
        exercise_dts,
        fwd_swap_rates,
        market_strikes,
        market_volatilities,
        vol_function_type,
    )

    if 1 == 0:  # plot_graphs:

        tol = 1e-4
        swaption_surface.check_calibration(False, tol)

        swaption_surface.plot_vol_curves()

        # plt.figure()

        # mins = 0.5
        # maxs = 5.0

        # dbns = swaption_surface.implied_dbns(mins, maxs, 1000)

        # for i in range(0, len(dbns)):
        #     expiry_dt_str = str(equity_surface._expiry_dts[i])
        #     plt.plot(dbns[i]._x, dbns[i]._densitydx, label = expiry_dt_str)
        #     plt.title(vol_function_type)
        #     plt.legend()
        #     print("SUM:", dbns[i].sum())


########################################################################################

test_FinSwaptionVolSurface1()

########################################################################################

