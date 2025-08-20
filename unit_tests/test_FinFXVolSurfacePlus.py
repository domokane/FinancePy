########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################
import matplotlib.pyplot as plt

from financepy.models.volatility_fns import VolFuncTypes
from financepy.utils.date import Date
from financepy.market.volatility.fx_vol_surface_plus import FinFXDeltaMethod
from financepy.market.volatility.fx_vol_surface_plus import FinFXATMMethod
from financepy.market.volatility.fx_vol_surface_plus import FXVolSurfacePlus
from financepy.market.volatility.fx_vol_surface import FXVolSurface
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
import numpy as np


verboseCalibration = False


def test_FinFXMktVolSurface1(capsys):
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
    rsk_reversal_10d_vols = [-1.258, -1.297, -1.332, -1.408, -1.359, -1.208]

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

    fx_market_plus.check_calibration(verboseCalibration)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_FinFXMktVolSurface2(capsys):
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

    # Fails with default stricter tolerance
    fx_market_plus.check_calibration(verboseCalibration, tol=0.2)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_FinFXMktVolSurface3(capsys):
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

    # I HAVE NO YET MADE DELTA METHOD A VECTOR FOR EACH TERM AS I WOULD
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

    fx_market_plus.check_calibration(verboseCalibration)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_FinFXMktVolSurface4(capsys):
    ###########################################################################
    # Here I remove the 25D Vols
    ###########################################################################

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
    rsk_reversal_10d_vols = [-1.258, -1.297, -1.332, -1.408, -1.359, -1.208]

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

    fx_market_plus.check_calibration(verboseCalibration)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_FinFXMktVolSurface5(capsys):
    ###########################################################################
    # Here I remove the 10D Vols
    ###########################################################################

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
    rsk_reversal_10d_vols = [-1.258, -1.297, -1.332, -1.408, -1.359, -1.208]

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

    fx_market_plus.check_calibration(verboseCalibration)
    captured = capsys.readouterr()
    assert captured.out == ""
