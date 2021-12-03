###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from financepy.models.volatility_fns import vol_function_clark
from financepy.market.volatility.fx_vol_surface import FXVolSurface
from financepy.market.volatility.fx_vol_surface import FinFXATMMethod
from financepy.market.volatility.fx_vol_surface import FinFXDeltaMethod
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinOptionImpliedDbn():

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
        # print("EURUSD EXAMPLE CLARK")

        valuation_date = Date(10, 4, 2020)

        forName = "EUR"
        domName = "USD"
        forCCRate = 0.03460  # EUR
        domCCRate = 0.02940  # USD

        dom_discount_curve = DiscountCurveFlat(valuation_date, domCCRate)
        for_discount_curve = DiscountCurveFlat(valuation_date, forCCRate)

        currency_pair = forName + domName
        spot_fx_rate = 1.3465

        tenors = ['1M', '2M', '3M', '6M', '1Y', '2Y']
        atm_vols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        marketStrangle25DeltaVols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        riskReversal25DeltaVols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]

        notional_currency = forName

        atmMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL
        deltaMethod = FinFXDeltaMethod.SPOT_DELTA

        fxMarket = FXVolSurface(valuation_date,
                                spot_fx_rate,
                                currency_pair,
                                notional_currency,
                                dom_discount_curve,
                                for_discount_curve,
                                tenors,
                                atm_vols,
                                marketStrangle25DeltaVols,
                                riskReversal25DeltaVols,
                                atmMethod,
                                deltaMethod)

#        fxMarket.check_calibration(True)

        PLOT_GRAPHS = False
        if PLOT_GRAPHS:
            fxMarket.plot_vol_curves()

        for iTenor in range(0, len(fxMarket._tenors)):

            F = fxMarket._F0T[iTenor]
            texp = fxMarket._texp[iTenor]

            startFX = F * 0.05
            endFX = F * 5.0

            num_steps = 10000
            dFX = (endFX - startFX) / num_steps

            domDF = dom_discount_curve._df(texp)
            forDF = for_discount_curve._df(texp)

            rd = -np.log(domDF) / texp
            rf = -np.log(forDF) / texp

            params = fxMarket._parameters[iTenor]

            strikes = []
            vols = []

            for iK in range(0, num_steps):
                strike = startFX + iK*dFX
                vol = vol_function_clark(params, F, strike, texp)
                strikes.append(strike)
                vols.append(vol)

            strikes = np.array(strikes)
            vols = np.array(vols)

#            dbn = optionImpliedDbn(spot_fx_rate, texp, rd, rf, strikes, vols)
#            print("SUM:", dbn.sum())
#            plt.figure()
#            plt.plot(dbn._x, dbn._densitydx)

###############################################################################


test_FinOptionImpliedDbn()
testCases.compareTestCases()
