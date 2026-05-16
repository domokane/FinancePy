# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.market.volatility.fx_vol_surface import FinFXDeltaMethod
from financepy.market.volatility.fx_vol_surface import FinFXATMMethod
from financepy.market.volatility.fx_vol_surface import FXVolSurface
from financepy.models.volatility_fns import vol_function_clark
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################

PLOT_GRAPHS = False


def test_fin_option_implied_dbn():

    if True:

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
        atm_vols = np.array([21.00, 21.00, 20.750, 19.400, 18.250, 17.677])
        mkt_strangle_25d_vols = np.array([0.65, 0.75, 0.85, 0.90, 0.95, 0.85])
        rsk_reversal_25d_vols = np.array([-0.20, -0.25, -0.30, -0.50, -0.60, -0.562])

        notional_currency = for_name

        atm_method = FinFXATMMethod.FWD_DELTA_NEUTRAL
        delta_method = FinFXDeltaMethod.SPOT_DELTA

        fx_market = FXVolSurface(
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
            atm_method,
            delta_method,
        )

        #        fx_market.check_calibration(True)

        if PLOT_GRAPHS:
            fx_market.plot_vol_curves()

        for i_tenor in range(0, len(fx_market.tenors)):

            f = fx_market.fwd[i_tenor]
            t_exp = fx_market.t_exp[i_tenor]

            start_fx = f * 0.05
            end_fx = f * 5.0

            num_steps = 10000
            d_fx = (end_fx - start_fx) / num_steps

            #            dom_df = domestic_curve.df_t(t_exp)
            #            for_df = foreign_curve.df_t(t_exp)
            #            r_d = -np.log(dom_df) / t_exp
            #            r_f = -np.log(for_df) / t_exp

            params = fx_market.parameters[i_tenor]

            strikes = []
            vols = []

            for i_k in range(0, num_steps):
                strike = start_fx + i_k * d_fx
                vol = vol_function_clark(params, f, strike, t_exp)
                strikes.append(strike)
                vols.append(vol)

            strikes = np.array(strikes)
            vols = np.array(vols)


########################################################################################

#            dbn = optionImpliedDbn(spot_fx_rate, t_exp, rd, rf, strikes, vols)
#            print("SUM:", dbn.sum())
#            plt.figure()
#            plt.plot(dbn._x, dbn._densitydx)


test_fin_option_implied_dbn()
test_cases.compare_test_cases()
