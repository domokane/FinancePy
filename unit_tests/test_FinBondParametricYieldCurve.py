# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import os
import datetime as dt
import numpy as np
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
    mean_err, max_err = fitted_curve.errors()
    print("Poly 5", mean_err, max_err)

    assert mean_err < 7.0
    assert max_err < 17.0


########################################################################################


def test_nelson_siegel():

    curve_fit_method = CurveFitNelsonSiegel()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    mean_err, max_err = fitted_curve.errors()
    print("NS", mean_err, max_err)

    assert mean_err < 9.0
    assert max_err < 30.0


########################################################################################


def test_svensson():

    curve_fit_method = CurveFitSvensson()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    mean_err, max_err = fitted_curve.errors()
    print("Svenson", mean_err, max_err)

    assert mean_err < 3.26
    assert max_err < 9.090

########################################################################################


def test_interp_yield():

    curve_fit_method = CurveFitBSpline()
    fitted_curve = BondParametricYieldCurve(
        settle_dt, bonds, ylds, curve_fit_method)

    mat_dt = Date(19, 9, 2030)
    interp_yield = fitted_curve.interp_yield(mat_dt)

    assert round(float(interp_yield), 4) == 0.0260

