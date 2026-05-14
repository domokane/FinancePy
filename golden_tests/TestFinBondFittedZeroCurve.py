# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os

import datetime as dt
import pandas as pd
import numpy as np

import add_fp_to_path

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date, from_datetime

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_fitted_zero_curve import BondFittedZeroCurve

from financepy.products.bonds.curve_fits import CurveFitPolynomial
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegel
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegelSvensson
from financepy.products.bonds.curve_fits import CurveFitBSpline

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

SHOW_PLOTS = False

########################################################################################


def test_bond_fitted_zero_curve():

    path = os.path.join(
        os.path.dirname(__file__), ".//data//gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

    freq_type = FrequencyTypes.SEMI_ANNUAL
    settle_dt = Date(19, 9, 2012)
    ex_div_days = 0

    dc_type = DayCountTypes.ACT_ACT_ICMA

    bonds = []
    dirty_prices = []

    for _, bond_row in bond_dataframe.iterrows():

        date_string = bond_row["maturity"]
        mat_dt_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_dt_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)

        coupon = bond_row["coupon"] / 100.0
        clean_price = bond_row["mid"]

        bond = Bond(
            issue_dt, maturity_dt, coupon, freq_type, dc_type, ex_div_days
        )

        accrued = bond.accrued_interest(settle_dt)

        dirty_price = clean_price + accrued

        bonds.append(bond)

        dirty_prices.append(dirty_price)

    curve_fitter = CurveFitPolynomial(4)

    exact_bond_curve = BondFittedZeroCurve(
        settle_dt, bonds, dirty_prices, curve_fitter
    )

    if SHOW_PLOTS == True:
        exact_bond_curve.plot_zero_rate("GBP Zero Rate Curve")
        exact_bond_curve.plot_fwd_rate("GBP Fwd Rate Curve")
        exact_bond_curve.plot_par_rate("GBP Par Rate Curve")

    total_err = 0.0

    for i_bond in range(0, len(bonds)):

        dirty_price_1 = bonds[i_bond].dirty_price_from_discount_curve(
            settle_dt, exact_bond_curve
        )

        dirty_price_2 = dirty_prices[i_bond]

        err = dirty_price_1 - dirty_price_2

        # print(i_bond, dirty_price_1, dirty_price_2, err)

        accrued = bond.accrued_interest(settle_dt)

        dirty_price = clean_price + accrued

        total_err += np.abs(err)

    # print("Total Error:", total_err)

    ###############################################################################


test_bond_fitted_zero_curve()

test_cases.compare_test_cases()
