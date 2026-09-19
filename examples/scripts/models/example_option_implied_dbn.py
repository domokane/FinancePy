# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.utils.global_types import FXDeltaMethodTypes
from financepy.utils.global_types import FXATMMethodTypes
from financepy.market.volatility.fx_vol_surface import FXVolSurface
from financepy.models.volatility_fns import vol_function_clark
from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

# ============================================================================
# FINANCEPY EXAMPLES - Option Implied Dbn
# ============================================================================



########################################################################################

PLOT_GRAPHS = False




########################################################################################

# ============================================================================
# 1. FIN OPTION IMPLIED DBN
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN OPTION IMPLIED DBN")
print("=" * 78)

if True:

    # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
    # print("EURUSD EXAMPLE CLARK")

    value_dt = Date(10, 4, 2020)
    for_name = "EUR"
    dom_name = "USD"
    for_cc_rate = 0.03460  # EUR
    dom_cc_rate = 0.02940  # USD

    domestic_curve = FlatDiscountCurve(value_dt, dom_cc_rate)
    foreign_curve = FlatDiscountCurve(value_dt, for_cc_rate)

    currency_pair = for_name + dom_name
    spot_fx_rate = 1.3465

    tenors = ["1M", "2M", "3M", "6M", "1Y", "2Y"]
    atm_vols = np.array([21.00, 21.00, 20.750, 19.400, 18.250, 17.677])
    mkt_strangle_25d_vols = np.array([0.65, 0.75, 0.85, 0.90, 0.95, 0.85])
    rsk_reversal_25d_vols = np.array([-0.20, -0.25, -0.30, -0.50, -0.60, -0.562])

    notional_currency = for_name

    atm_method = FXATMMethodTypes.FWD_DELTA_NEUTRAL
    delta_method = FXDeltaMethodTypes.SPOT_DELTA

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

    if PLOT_GRAPHS:
        fx_market.plot_vol_curves()

    for i_tenor in range(0, len(fx_market.tenors)):

        f = fx_market.fwd[i_tenor]
        t_exp = fx_market.t_exp[i_tenor]

        start_fx = f * 0.05
        end_fx = f * 5.0

        num_steps = 10000
        d_fx = (end_fx - start_fx) / num_steps

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

