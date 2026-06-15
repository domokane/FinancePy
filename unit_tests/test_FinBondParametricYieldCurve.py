# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import datetime as dt

import pandas as pd

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date, from_datetime
from financepy.products.bonds.bond import Bond
from financepy.market.curves.bond_parametric_yield_curve import BondParametricYieldCurve
from financepy.market.curves.curve_fits import CurveFitPolynomial
from financepy.market.curves.curve_fits import CurveFitBSpline
from financepy.market.curves.curve_fits import CurveFitNelsonSiegel
from financepy.market.curves.curve_fits import CurveFitSvensson

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

########################################################################################


def test_poly():

    curve_fit_method = CurveFitPolynomial(5)
    fitted_curve = BondParametricYieldCurve(settle_dt, bonds, ylds, curve_fit_method)

    coeffs = fitted_curve.curve_fit.coeffs

    assert round(coeffs[0]*1e7, 4) == 4.1581
    assert round(coeffs[1], 4) == 0.0641
    assert round(coeffs[2], 4) == 0.2034
    assert round(coeffs[3], 4) == -0.7883
    assert round(coeffs[4], 4) == 0.8984
    assert round(coeffs[5], 4) == -0.3454


########################################################################################


def test_nelson_siegel():

    curve_fit_method = CurveFitNelsonSiegel()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    assert round(fitted_curve.curve_fit.beta_1, 3) == -0.094
    assert round(fitted_curve.curve_fit.beta_2, 3) == 0.092
    assert round(fitted_curve.curve_fit.beta_3, 3) == 0.259
    assert round(fitted_curve.curve_fit.tau, 3) == 0.755


########################################################################################


def test_svensson():

    curve_fit_method = CurveFitSvensson()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    assert round(fitted_curve.curve_fit.beta_1, 4) == -0.6137
    assert round(fitted_curve.curve_fit.beta_2, 4) == +0.6121
    assert round(fitted_curve.curve_fit.beta_3, 4) == 0.6353
    assert round(fitted_curve.curve_fit.beta_4, 4) == 5.0000
    assert round(fitted_curve.curve_fit.tau_1, 4) == 0.9529
    assert round(fitted_curve.curve_fit.tau_2, 4) == 25.4461


########################################################################################


def test_interp_yield():

    curve_fit_method = CurveFitBSpline()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    mat_dt = Date(19, 9, 2030)
    interp_yield = fitted_curve.interp_yield(mat_dt)

    assert round(float(interp_yield), 8) == 0.02594837

# test_poly()
