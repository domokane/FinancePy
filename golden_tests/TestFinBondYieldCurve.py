###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import datetime as dt
import os

import sys

sys.path.append("..")

import pandas as pd

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegelSvensson
from financepy.products.bonds.curve_fits import CurveFitNelsonSiegel
from financepy.products.bonds.curve_fits import CurveFitBSpline
from financepy.products.bonds.curve_fits import CurveFitPolynomial
from financepy.products.bonds.bond_yield_curve import BondYieldCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
###############################################################################


def test_BondYieldCurve():

    ###########################################################################

    path = os.path.join(
        os.path.dirname(__file__), "./data/gilt_bond_prices.txt"
    )
    bond_dataframe = pd.read_csv(path, sep="\t")
    bond_dataframe["mid"] = 0.5 * (
        bond_dataframe["bid"] + bond_dataframe["ask"]
    )

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

    ###############################################################################

    curve_fitter = CurveFitPolynomial()
    fitted_curve1 = BondYieldCurve(settlement, bonds, ylds, curve_fitter)
    #    fitted_curve1.display("GBP Yield Curve")

    curve_fitter = CurveFitPolynomial(5)
    fitted_curve2 = BondYieldCurve(settlement, bonds, ylds, curve_fitter)
    #    fitted_curve2.display("GBP Yield Curve")

    curve_fitter = CurveFitNelsonSiegel()
    fitted_curve3 = BondYieldCurve(settlement, bonds, ylds, curve_fitter)
    #    fitted_curve3.display("GBP Yield Curve")

    curve_fitter = CurveFitNelsonSiegelSvensson()
    fitted_curve4 = BondYieldCurve(settlement, bonds, ylds, curve_fitter)
    #    fitted_curve4.display("GBP Yield Curve")

    curve_fitter = CurveFitBSpline()
    fitted_curve5 = BondYieldCurve(settlement, bonds, ylds, curve_fitter)
    #    fitted_curve5.display("GBP Yield Curve")

    ###############################################################################

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

    ###############################################################################

    maturity_dt = Date(19, 9, 2030)
    interp_yield = fitted_curve5.interp_yield(maturity_dt)
    test_cases.print(maturity_dt, interp_yield)


###############################################################################


test_BondYieldCurve()
test_cases.compareTestCases()
