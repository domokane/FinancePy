# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import os
import datetime as dt

import numpy as np
import pandas as pd

import add_fp_to_path

from financepy.utils.date import Date, from_datetime
from financepy.products.bonds.bond import Bond
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.market.curves.curve_fits import CurveFitTypes
from financepy.market.curves import BondParametricDiscountCurve



SHOW_PLOTS = False

########################################################################################


def test_bond_parametric_discount_curve():

    path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    settle_dt = Date(19, 9, 2012)

    bonds = []
    clean_prices = []

    for _, bond in bond_dataframe.iterrows():

        date_string = bond["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bond["coupon"] / 100.0
        clean_price = bond["mid"]
        clean_prices.append(clean_price)
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
        bonds.append(bond)

    curve_fitters = [
        CurveFitTypes.CUBIC_POLYNOMIAL,
        CurveFitTypes.QUINTIC_POLYNOMIAL,
        CurveFitTypes.NELSON_SIEGEL,
        CurveFitTypes.NELSON_SIEGEL_SVENSSON,
        CurveFitTypes.BSPLINE,
    ]

    print("CURVE_FITTER", "Time", "DISCOUNT FACTOR", "ZERO_RATE", "FORWARD RATE")
    times = np.linspace(1e-6, 50, 100)

    for curve_fitter in curve_fitters:

        name = curve_fitter.name

        fitted_curve = BondParametricDiscountCurve(
            settle_dt,
            bonds,
            clean_prices,
            curve_fitter,
        )

        if SHOW_PLOTS:
            fitted_curve.plot_bond_yield_fit(name + " Bond Yield Fit")

        # PRINT TO TEST CASES FILE
        dfs = fitted_curve.df_t(times)
        fwds = fitted_curve.fwd_rate_inst_t(times)
        zeros = fitted_curve.zero_rate_cc_t(times)

        print("CURVE_FITTER", "TIME", "DISCOUNT FACTOR", "ZERO_RATE", "FORWARD RATE")
        for i, t in enumerate(times):
            print(name, t, dfs[i], zeros[i], fwds[i])

        rms_yield_err = fitted_curve.rms_yield_error()
        rms_price_err = fitted_curve.rms_price_error()

        print("CURVE_FITTER", "ERROR TYPE", "VALUE")
        print(name, "RMS_YIELD_ERROR", rms_yield_err)
        print(name, "RMS_PRICE_ERROR", rms_price_err)


########################################################################################

test_bond_parametric_discount_curve()
