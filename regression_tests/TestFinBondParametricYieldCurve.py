# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import datetime as dt

import pandas as pd

import add_fp_to_path

from financepy.utils.date import Date, from_datetime
from financepy.products.bonds.bond import Bond
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.market.curves import CurveFitSvensson
from financepy.market.curves import CurveFitNelsonSiegel
from financepy.market.curves import CurveFitBSpline
from financepy.market.curves import CurveFitPolynomial
from financepy.market.curves import BondParametricYieldCurve

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

SHOW_PLOTS = False

########################################################################################


def test_bond_parametric_yield_curve():

    path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    settle_dt = Date(19, 9, 2012)

    bonds = []
    ylds = []

    for _, bond in bond_dataframe.iterrows():

        date_string = bond["maturity"]
        mat_date_time = dt.datetime.strptime(date_string, "%d-%b-%y")
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt.d, maturity_dt.m, 2000)
        coupon = bond["coupon"] / 100.0
        clean_price = bond["mid"]
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
        yld = bond.yield_to_maturity(settle_dt, clean_price)
        bonds.append(bond)
        ylds.append(yld)

    curve_fitter = CurveFitPolynomial()
    fitted_curve1 = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fitter)

    if SHOW_PLOTS:
        fitted_curve1.plot("GBP Yield Curve")

    curve_fitter = CurveFitPolynomial(5)
    fitted_curve2 = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fitter)
    if SHOW_PLOTS:
        fitted_curve2.plot("GBP Yield Curve")

    curve_fitter = CurveFitNelsonSiegel()
    fitted_curve3 = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fitter)
    if SHOW_PLOTS:
        fitted_curve3.plot("GBP Yield Curve")

    curve_fitter = CurveFitSvensson()
    fitted_curve4 = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fitter)
    if SHOW_PLOTS:
        fitted_curve4.plot("GBP Yield Curve")

    curve_fitter = CurveFitBSpline()
    fitted_curve5 = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fitter)
    if SHOW_PLOTS:
        fitted_curve5.plot("GBP Yield Curve")

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("values", fitted_curve1.curve_fit.coeffs)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("values", fitted_curve2.curve_fit.coeffs)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("beta_1", fitted_curve3.curve_fit.beta_1)
    test_cases.print("beta_2", fitted_curve3.curve_fit.beta_2)
    test_cases.print("beta_3", fitted_curve3.curve_fit.beta_3)
    test_cases.print("tau", fitted_curve3.curve_fit.tau)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("beta_1", fitted_curve4.curve_fit.beta_1)
    test_cases.print("beta_2", fitted_curve4.curve_fit.beta_2)
    test_cases.print("beta_3", fitted_curve4.curve_fit.beta_3)
    test_cases.print("beta_4", fitted_curve4.curve_fit.beta_4)
    test_cases.print("tau_1", fitted_curve4.curve_fit.tau_1)
    test_cases.print("tau_2", fitted_curve4.curve_fit.tau_2)

    maturity_dt = Date(19, 9, 2030)
    interp_yield = fitted_curve5.interp_yield(maturity_dt)
    test_cases.print(maturity_dt, interp_yield)


########################################################################################

test_bond_parametric_yield_curve()
test_cases.compare_test_cases()
