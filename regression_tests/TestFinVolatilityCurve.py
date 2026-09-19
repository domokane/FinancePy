# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.market.volatility.equity_vol_curve import (
    EquityVolCurveExtrapTypes,
)
from financepy.market.volatility.equity_vol_curve import (
    EquityVolCurveInterpTypes,
)

from financepy.market.volatility.equity_vol_curve import (
    EquityVolCurveDeltaTypes,
)
from financepy.utils.global_types import OptionTypes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

########################################################################################

s_min = 50.0
s_max = 200
n_intervals = 100


def test_fin_volatility_curve():
    strikes = np.linspace(70, 130, 7)
    mkt_vols = np.array([0.23, 0.24, 0.267, 0.29, 0.31, 0.33, 0.35])

    interp_strikes = np.linspace(50, 150, 100)

    for interp_type in EquityVolCurveInterpTypes:

        for extrap_type in EquityVolCurveExtrapTypes:

            vol_curve = EquityVolCurve(
                strikes,
                mkt_vols,
                s=100.0,
                t_exp=0.5,
                r=0.04,
                q=0.02,
                interp=interp_type,
                extrap=extrap_type,
                smoothing=1.0e-4,
                derivative_bump=1.0e-4,
            )

            interp_vols = vol_curve.volatility(interp_strikes)

            delta1 = vol_curve.delta(strikes, OptionTypes.EUROPEAN_CALL, EquityVolCurveDeltaTypes.STICKY_STRIKE)
            delta2 = vol_curve.delta(strikes, OptionTypes.EUROPEAN_CALL, EquityVolCurveDeltaTypes.STICKY_MONEYNESS)
            title = str(interp_type) + " - " + str(extrap_type)
            pdf_strikes, pdf_values = vol_curve.pdf(s_min, s_max, n_intervals)

            if PLOT_GRAPHS:
                import matplotlib.pyplot as plt

                plt.figure()
                plt.plot(strikes, mkt_vols * 100, "o", label="Market")
                plt.plot(interp_strikes, interp_vols * 100, label="Interpolated")
                plt.xlabel("Strike")
                plt.ylabel("Volatility")
                plt.title(title)
                plt.legend()
                plt.grid()
                plt.show()

                plt.figure()
                plt.plot(strikes, delta1, "o-", label="Sticky Strike")
                plt.plot(strikes, delta2, "o-", label="Sticky Moneyness")
                plt.xlabel("Strike")
                plt.ylabel("Delta")
                plt.title(title)
                plt.legend()
                plt.grid()
                plt.show()

                plt.figure()
                plt.plot(pdf_strikes, pdf_values, label="Interpolated")
                plt.xlabel("Strike")
                plt.ylabel("Density g(K) dK")
                plt.title(title)
                plt.legend()
                plt.grid()
                plt.show()


########################################################################################

test_fin_volatility_curve()
test_cases.compare_test_cases()
