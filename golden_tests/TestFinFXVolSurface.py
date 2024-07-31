###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import matplotlib.pyplot as plt
from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.volatility_fns import VolFuncTypes
from financepy.utils.date import Date
from financepy.market.volatility.fx_vol_surface import FinFXDeltaMethod
from financepy.market.volatility.fx_vol_surface import FinFXATMMethod
from financepy.market.volatility.fx_vol_surface import FXVolSurface
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)


###############################################################################

PLOT_GRAPHS = False

###############################################################################


def test_FinFXMktVolSurface1(verboseCalibration):

    ###########################################################################

    if 1 == 1:

        # Example from Book extract by Iain Clarke using Tables 3.3 and 3.4
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

        tenors = ['1M', '2M', '3M', '6M', '1Y', '2Y']
        atm_vols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        mkt_strangle_25d_vols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        rsk_reversal_25d_vols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]

        notional_currency = for_name

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.SPOT_DELTA
        vol_functionType = VolFuncTypes.CLARK

        fx_market = FXVolSurface(value_dt,
                                spot_fx_rate,
                                currency_pair,
                                notional_currency,
                                domestic_curve,
                                foreign_curve,
                                tenors,
                                atm_vols,
                                mkt_strangle_25d_vols,
                                rsk_reversal_25d_vols,
                                atm_method,
                                delta_method,
                                vol_functionType)

        fx_market.check_calibration(verboseCalibration)

        # EXPLORE AND TEST DIFFERENT CATEGORICAL PARAMETERS
        # for atm_method in FinFXATMMethod:
        #     for delta_method in FinFXDeltaMethod:
        #         for vol_functionType in VolFuncTypes:

        #             fx_market = FinFXVolSurface(value_dt,
        #                                        spot_fx_rate,
        #                                        currency_pair,
        #                                        notional_currency,
        #                                        domestic_curve,
        #                                        foreign_curve,
        #                                        tenors,
        #                                        atm_vols,
        #                                        mkt_strangle_25d_vols,
        #                                        rsk_reversal_25d_vols,
        #                                        atm_method,
        #                                        delta_method,
        #                                        vol_functionType)

        #             fx_market.check_calibration(verboseCalibration)

        if PLOT_GRAPHS:

            fx_market.plot_vol_curves()

            dbns = fx_market.implied_dbns(0.00001, 5.0, 10000)

            for i in range(0, len(dbns)):
                plt.plot(dbns[i]._x, dbns[i]._densitydx)
                plt.title(vol_functionType)
                print("SUM:", dbns[i].sum())

    ###########################################################################


def test_FinFXMktVolSurface2(verboseCalibration):

    # print("==============================================================")

    # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
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

    tenors = ['1M', '2M', '3M', '6M', '1Y', '2Y']
    atm_vols = [21.50, 20.50, 19.85, 18.00, 15.95, 14.009]
    mkt_strangle_25d_vols = [0.35, 0.325, 0.300, 0.225, 0.175, 0.100]
    rsk_reversal_25d_vols = [-8.350, -8.650, -8.950, -9.250, -9.550, -9.500]

    notional_currency = for_name

    atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ
    delta_method = FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ

    fx_market = FXVolSurface(value_dt,
                            spot_fx_rate,
                            currency_pair,
                            notional_currency,
                            domestic_curve,
                            foreign_curve,
                            tenors,
                            atm_vols,
                            mkt_strangle_25d_vols,
                            rsk_reversal_25d_vols,
                            atm_method,
                            delta_method)

    fx_market.check_calibration(verboseCalibration)

    if PLOT_GRAPHS:
        fx_market.plot_vol_curves()

#   print("==================================================================")

#   ###########################################################################


def test_FinFXMktVolSurface3(verboseCalibration):

    # EURUSD Example from Paper by Uwe Wystup using Tables 4
    #        print("EURUSD EXAMPLE WYSTUP")

    value_dt = Date(20, 1, 2009)

    for_name = "EUR"
    dom_name = "USD"
    for_cc_rate = 0.020113  # EUR
    dom_cc_rate = 0.003525  # USD

    domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
    foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

    currency_pair = for_name + dom_name
    spot_fx_rate = 1.3088

    tenors = ['1M']
    atm_vols = [21.6215]
    mkt_strangle_25d_vols = [0.7375]
    rsk_reversal_25d_vols = [-0.50]

    notional_currency = for_name

    atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
    delta_method = FinFXDeltaMethod.SPOT_DELTA

    fx_market = FXVolSurface(value_dt,
                            spot_fx_rate,
                            currency_pair,
                            notional_currency,
                            domestic_curve,
                            foreign_curve,
                            tenors,
                            atm_vols,
                            mkt_strangle_25d_vols,
                            rsk_reversal_25d_vols,
                            atm_method,
                            delta_method)

    fx_market.check_calibration(verboseCalibration)

    if PLOT_GRAPHS:
        fx_market.plot_vol_curves()

    ###########################################################################


def test_FinFXMktVolSurface4(verboseCalibration):

    # USDJPY Example from Paper by Uwe Wystup using Tables 4
    #        print("USDJPY EXAMPLE WYSTUP")

    value_dt = Date(20, 1, 2009)

    for_name = "USD"
    dom_name = "JPY"
    for_cc_rate = 0.003525  # USD
    dom_cc_rate = 0.0042875  # JPY

    domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
    foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

    currency_pair = for_name + dom_name
    spot_fx_rate = 90.68

    tenors = ['1M']
    atm_vols = [21.00]
    mkt_strangle_25d_vols = [0.184]
    rsk_reversal_25d_vols = [-5.30]

    notional_currency = for_name

    atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
    delta_method = FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ

    fx_market = FXVolSurface(value_dt,
                            spot_fx_rate,
                            currency_pair,
                            notional_currency,
                            domestic_curve,
                            foreign_curve,
                            tenors,
                            atm_vols,
                            mkt_strangle_25d_vols,
                            rsk_reversal_25d_vols,
                            atm_method,
                            delta_method)

    fx_market.check_calibration(verboseCalibration)

    if PLOT_GRAPHS:
        fx_market.plot_vol_curves()

    #    test_cases.header("value", "delta")
    #    test_cases.print(value, delta)

###############################################################################


if __name__ == '__main__':

    start = time.time()

    verboseCalibration = False

    test_FinFXMktVolSurface1(verboseCalibration)
    test_FinFXMktVolSurface2(verboseCalibration)
    test_FinFXMktVolSurface3(verboseCalibration)
    test_FinFXMktVolSurface4(verboseCalibration)

    end = time.time()

    elapsed = end - start
#    print("Elapsed Time:", elapsed)
    test_cases.compareTestCases()
