###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################
import os
import pandas as pd
import datetime as dt

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date, from_datetime
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_yield_curve import BondYieldCurve
from financepy.products.bonds.curve_fits import CurveFitPolynomial
from financepy.products.bonds.curve_fits import CurveFitBSpline
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegel
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegelSvensson

path = os.path.join(os.path.dirname(__file__), "./data/gilt_bond_prices.txt")
bond_dataframe = pd.read_csv(path, sep="\t")
bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
settlement = Date(19, 9, 2012)

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
    yld = bond.yield_to_maturity(settlement, clean_price)
    bonds.append(bond)
    ylds.append(yld)


def test_poly():
    curveFitMethod = CurveFitPolynomial(5)
    fitted_curve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    coeffs = fitted_curve.curve_fit.coeffs
    assert round(coeffs[0] * 1e9, 4) == -1.4477
    assert round(coeffs[1] * 1e7, 4) == 1.7840
    assert round(coeffs[2] * 1e6, 4) == -7.4147
    assert round(coeffs[3] * 1e5, 4) == 9.0622
    assert round(coeffs[4] * 1e3, 4) == 1.3536
    assert round(coeffs[5] * 1e7, 4) == 4.1514


def test_nelson_siegel():
    curveFitMethod = CurveFitNelsonSiegel()
    fitted_curve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    assert round(fitted_curve.curve_fit.beta_1, 3) == -0.094
    assert round(fitted_curve.curve_fit.beta_2, 3) == 0.092
    assert round(fitted_curve.curve_fit.beta_3, 3) == 0.259
    assert round(fitted_curve.curve_fit.tau, 1) == 35.8


def test_nelson_siegel_svensson():
    curveFitMethod = CurveFitNelsonSiegelSvensson()
    fitted_curve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    assert round(fitted_curve.curve_fit.beta_1, 4) == 0.0460
    assert round(fitted_curve.curve_fit.beta_2, 4) == -0.0433
    assert round(fitted_curve.curve_fit.beta_3, 4) == -0.0523
    assert round(fitted_curve.curve_fit.beta_4, 4) == -0.0376
    assert round(fitted_curve.curve_fit.tau_1, 3) == 3.177
    assert round(fitted_curve.curve_fit.tau_2, 4) == 100.0000


def test_interp_yield():
    curveFitMethod = CurveFitBSpline()
    fitted_curve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    maturity_dt = Date(19, 9, 2030)
    interp_yield = fitted_curve.interp_yield(maturity_dt)
    assert round(float(interp_yield), 8) == 0.02601858
