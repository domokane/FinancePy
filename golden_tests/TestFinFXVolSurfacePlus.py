# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import sys

sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.models.volatility_fns import VolFuncTypes
from financepy.utils.date import Date
from financepy.market.volatility.fx_vol_surface_plus import FinFXDeltaMethod
from financepy.market.volatility.fx_vol_surface_plus import FinFXATMMethod
from financepy.market.volatility.fx_vol_surface_plus import FXVolSurfacePlus
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat


test_cases = FinTestCases(__file__, global_test_case_mode)


plot_graphs = False

# TODO: ADD LOGGING TO TEST CASES

########################################################################################


def test_fin_fx_mkt_vol_surface1(verbose_calibration):

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
        # print("EURUSD EXAMPLE CLARK")

        value_dt = Date(10, 4, 2020)

        for_name = "EUR"
        dom_name = "USD"
        for_cc_rate = 0.03460  # EUR
        dom_cc_rate = 0.02940  # USD

        domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
        foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

        currency_pair = for_name + dom_name
        spot_fx_rate = 1.3465

        tenors = ["1M", "2M", "3M", "6M", "1Y", "2Y"]
        atm_vols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        mkt_strangle_25d_vols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        rsk_reversal_25d_vols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]
        mkt_strangle_10d_vols = [2.433, 2.83, 3.228, 3.485, 3.806, 3.208]
        rsk_reversal_10d_vols = [
            -1.258,
            -1.297,
            -1.332,
            -1.408,
            -1.359,
            -1.208,
        ]

        notional_currency = for_name

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.SPOT_DELTA
        vol_function_type = VolFuncTypes.CLARK5
        alpha = 0.5  # FIT WINGS AT 10D if ALPHA = 1.0

        fx_market_plus = FXVolSurfacePlus(
            value_dt,
            spot_fx_rate,
            currency_pair,
            notional_currency,
            domestic_curve,
            foreign_curve,
            tenors,
            atm_vols,
            mkt_strangle_25d_vols,
            rsk_reversal_25d_vols,
            mkt_strangle_10d_vols,
            rsk_reversal_10d_vols,
            alpha,
            atm_method,
            delta_method,
            vol_function_type,
        )

        fx_market_plus.check_calibration(False)

        if 1 == 0:  # plot_graphs:

            fx_market_plus.plot_vol_curves()

            plt.figure()

            dbns = fx_market_plus.implied_dbns(0.5, 2.0, 1000)

            for i in range(0, len(dbns)):
                plt.plot(dbns[i]._x, dbns[i]._densitydx)
                plt.title(vol_function_type)
                print("SUM:", dbns[i].sum())


########################################################################################


def test_fin_fx_mkt_vol_surface2(verbose_calibration):

    # print("==============================================================")

    # Example from Book extract by Iain Clarke using Tables 3.3 and 3.4
    # print("EURJPY EXAMPLE CLARK")

    value_dt = Date(10, 4, 2020)

    for_name = "EUR"
    dom_name = "JPY"
    for_cc_rate = 0.0294  # EUR
    dom_cc_rate = 0.0171  # USD

    domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
    foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

    currency_pair = for_name + dom_name
    spot_fx_rate = 90.72

    tenors = ["1M", "2M", "3M", "6M", "1Y", "2Y"]
    atm_vols = [21.50, 20.50, 19.85, 18.00, 15.95, 14.009]
    mkt_strangle_25d_vols = [0.35, 0.325, 0.300, 0.225, 0.175, 0.100]
    rsk_reversal_25d_vols = [-8.350, -8.650, -8.950, -9.250, -9.550, -9.500]
    mkt_strangle_10d_vols = [3.704, 4.047, 4.396, 4.932, 5.726, 5.709]
    rsk_reversal_10d_vols = [
        -15.855,
        -16.467,
        -17.114,
        -17.882,
        -18.855,
        -18.217,
    ]
    alpha = 0.50  # Equally fit 10 and 25 Delta

    notional_currency = for_name

    atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ
    delta_method = FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ
    vol_function_type = VolFuncTypes.CLARK5

    fx_market_plus = FXVolSurfacePlus(
        value_dt,
        spot_fx_rate,
        currency_pair,
        notional_currency,
        domestic_curve,
        foreign_curve,
        tenors,
        atm_vols,
        mkt_strangle_25d_vols,
        rsk_reversal_25d_vols,
        mkt_strangle_10d_vols,
        rsk_reversal_10d_vols,
        alpha,
        atm_method,
        delta_method,
        vol_function_type,
    )

    #        fx_market_plus.check_calibration(True)

    if plot_graphs:
        fx_market_plus.plot_vol_curves()

        plt.figure()

        dbns = fx_market_plus.implied_dbns(30, 120, 1000)

        for i in range(0, len(dbns)):
            plt.plot(dbns[i]._x, dbns[i]._densitydx)
            plt.title(vol_function_type)
            print("SUM:", dbns[i].sum())


########################################################################################


def test_fin_fx_mkt_vol_surface3(verbose_calibration):

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 4.4 and 4.5
        # where we examine the calibration to a full surface in Chapter 4

        value_dt = Date(10, 4, 2020)

        for_name = "EUR"
        dom_name = "USD"
        for_cc_rate = 0.03460  # EUR
        dom_cc_rate = 0.02940  # USD

        domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
        foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

        currency_pair = for_name + dom_name
        spot_fx_rate = 1.3465

        tenors = ["1Y", "2Y"]
        atm_vols = [18.250, 17.677]
        mkt_strangle_25d_vols = [0.95, 0.85]
        rsk_reversal_25d_vols = [-0.60, -0.562]
        mkt_strangle_10d_vols = [3.806, 3.208]
        rsk_reversal_10d_vols = [-1.359, -1.208]

        notional_currency = for_name

        # I HAVE NOT YET MADE DELTA METHOD A VECTOR FOR EACH TERM AS I WOULD
        # NEED TO DO AS DESCRIBED IN CLARK PAGE 70

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.FORWARD_DELTA  # THIS IS DIFFERENT
        vol_function_type = VolFuncTypes.CLARK5
        alpha = 0.5  # FIT WINGS AT 10D if ALPHA = 1.0

        fx_market_plus = FXVolSurfacePlus(
            value_dt,
            spot_fx_rate,
            currency_pair,
            notional_currency,
            domestic_curve,
            foreign_curve,
            tenors,
            atm_vols,
            mkt_strangle_25d_vols,
            rsk_reversal_25d_vols,
            mkt_strangle_10d_vols,
            rsk_reversal_10d_vols,
            alpha,
            atm_method,
            delta_method,
            vol_function_type,
        )

        fx_market_plus.check_calibration(False)

        if 1 == 0:  # plot_graphs:

            fx_market_plus.plot_vol_curves()

            plt.figure()

            dbns = fx_market_plus.implied_dbns(0.5, 2.0, 1000)

            for i in range(0, len(dbns)):
                plt.plot(dbns[i]._x, dbns[i]._densitydx)
                plt.title(vol_function_type)
                print("SUM:", dbns[i].sum())

        # Test interpolation

        years = [1.0, 1.5, 2.0]
        dates = value_dt.add_years(years)

        strikes = np.linspace(1.0, 2.0, 20)

        if 1 == 0:
            vol_surface = []
            for k in strikes:
                vol_smile = []
                for dt in dates:
                    vol = fx_market_plus.volatility_from_strike_dt(k, dt)
                    vol_smile.append(vol * 100.0)

                    print(k, dt, vol * 100.0)
                vol_surface.append(vol_smile)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            x, y = np.meshgrid(years, strikes)
            zs = np.array(vol_surface)
            z = zs.reshape(x.shape)

            ax.plot_surface(x, y, z)

            ax.set_xlabel("Years")
            ax.set_ylabel("Strikes")
            ax.set_zlabel("Volatility")

            plt.show()

        deltas = np.linspace(0.10, 0.90, 17)

        if 1 == 0:
            vol_surface = []
            for delta in deltas:
                vol_smile = []
                for dt in dates:
                    (vol, k) = fx_market_plus.vol_from_delta_date(delta, dt)
                    vol_smile.append(vol * 100.0)
                    print(delta, k, dt, vol * 100.0)

                vol_surface.append(vol_smile)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            x, y = np.meshgrid(years, deltas)
            zs = np.array(vol_surface)
            z = zs.reshape(x.shape)

            ax.plot_surface(x, y, z)

            ax.set_xlabel("Years")
            ax.set_ylabel("Delta")
            ax.set_zlabel("Volatility")
            plt.title("EURUSD Volatility Surface")
            plt.show()


########################################################################################


def test_fin_fx_mkt_vol_surface4(verbose_calibration):

    # Here I remove the 25D Vols

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
        # print("EURUSD EXAMPLE CLARK")

        value_dt = Date(10, 4, 2020)

        for_name = "EUR"
        dom_name = "USD"
        for_cc_rate = 0.03460  # EUR
        dom_cc_rate = 0.02940  # USD

        domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
        foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

        currency_pair = for_name + dom_name
        spot_fx_rate = 1.3465

        tenors = ["1M", "2M", "3M", "6M", "1Y", "2Y"]
        atm_vols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        mkt_strangle_25d_vols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        rsk_reversal_25d_vols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]
        mkt_strangle_10d_vols = [2.433, 2.83, 3.228, 3.485, 3.806, 3.208]
        rsk_reversal_10d_vols = [
            -1.258,
            -1.297,
            -1.332,
            -1.408,
            -1.359,
            -1.208,
        ]

        mkt_strangle_25d_vols = None
        rsk_reversal_25d_vols = None

        notional_currency = for_name

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.SPOT_DELTA
        vol_function_type = VolFuncTypes.CLARK
        alpha = 0.50  # FIT WINGS AT 10D if ALPHA = 1.0

        fx_market_plus = FXVolSurfacePlus(
            value_dt,
            spot_fx_rate,
            currency_pair,
            notional_currency,
            domestic_curve,
            foreign_curve,
            tenors,
            atm_vols,
            mkt_strangle_25d_vols,
            rsk_reversal_25d_vols,
            mkt_strangle_10d_vols,
            rsk_reversal_10d_vols,
            alpha,
            atm_method,
            delta_method,
            vol_function_type,
        )

        fx_market_plus.check_calibration(False)

        years = [1.0 / 12.0, 2.0 / 12.0, 0.25, 0.5, 1.0, 2.0]

        dates = value_dt.add_years(years)

        deltas = np.linspace(0.10, 0.90, 17)

        if 1 == 0:
            vol_surface = []
            for delta in deltas:
                vol_smile = []
                for dt in dates:
                    (vol, k) = fx_market_plus.vol_from_delta_date(delta, dt)
                    vol_smile.append(vol * 100.0)
                    print(delta, k, dt, vol * 100.0)

                vol_surface.append(vol_smile)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
            x, y = np.meshgrid(years, deltas)
            zs = np.array(vol_surface)
            z = zs.reshape(x.shape)

            ax.plot_surface(x, y, z)

            ax.set_xlabel("Years")
            ax.set_ylabel("Delta")
            ax.set_zlabel("Volatility")
            plt.title("EURUSD Volatility Surface")
            plt.show()


########################################################################################


def test_fin_fx_mkt_vol_surface5(verbose_calibration):

    # Here I remove the 10D Vols

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
        # print("EURUSD EXAMPLE CLARK")

        value_dt = Date(10, 4, 2020)

        for_name = "EUR"
        dom_name = "USD"
        for_cc_rate = 0.03460  # EUR
        dom_cc_rate = 0.02940  # USD

        domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
        foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

        currency_pair = for_name + dom_name
        spot_fx_rate = 1.3465

        tenors = ["1M", "2M", "3M", "6M", "1Y", "2Y"]
        atm_vols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        mkt_strangle_25d_vols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        rsk_reversal_25d_vols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]
        mkt_strangle_10d_vols = [2.433, 2.83, 3.228, 3.485, 3.806, 3.208]
        rsk_reversal_10d_vols = [
            -1.258,
            -1.297,
            -1.332,
            -1.408,
            -1.359,
            -1.208,
        ]

        mkt_strangle_10d_vols = None
        rsk_reversal_10d_vols = None

        notional_currency = for_name

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.SPOT_DELTA
        vol_function_type = VolFuncTypes.CLARK
        alpha = 0.50  # FIT WINGS AT 10D if ALPHA = 1.0

        fx_market_plus = FXVolSurfacePlus(
            value_dt,
            spot_fx_rate,
            currency_pair,
            notional_currency,
            domestic_curve,
            foreign_curve,
            tenors,
            atm_vols,
            mkt_strangle_25d_vols,
            rsk_reversal_25d_vols,
            mkt_strangle_10d_vols,
            rsk_reversal_10d_vols,
            alpha,
            atm_method,
            delta_method,
            vol_function_type,
        )

        fx_market_plus.check_calibration(False)


########################################################################################

########################################################################################

if __name__ == "__main__":

    start = time.time()

    verbose_calibration = False

    test_fin_fx_mkt_vol_surface1(verbose_calibration)
    test_fin_fx_mkt_vol_surface2(verbose_calibration)
    test_fin_fx_mkt_vol_surface3(verbose_calibration)
    test_fin_fx_mkt_vol_surface4(verbose_calibration)
    test_fin_fx_mkt_vol_surface5(verbose_calibration)

    end = time.time()

    elapsed = end - start
    #    print("Elapsed Time:", elapsed)
    test_cases.compare_test_cases()
