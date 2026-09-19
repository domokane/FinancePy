# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import os
import datetime as dt
import pandas as pd


from financepy.utils.date import Date, from_datetime
from financepy.products.bonds.bond import Bond
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.market.curves import CurveFitTypes
from financepy.market.curves import BondParametricYieldCurve

# ============================================================================
# FINANCEPY EXAMPLES - Bond
# ============================================================================


SHOW_PLOTS = False

########################################################################################




########################################################################################

# ============================================================================
# 1. BOND PARAMETRIC YIELD CURVE
# ============================================================================
# What this section demonstrates:
# Solves for the yield that reproduces the observed bond price. This checks the inverse relationship between price and yield.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND PARAMETRIC YIELD CURVE")
print("=" * 78)

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

fit_type = CurveFitTypes.CUBIC_POLYNOMIAL
fitted_curve1 = BondParametricYieldCurve(settle_dt, bonds, ylds, fit_type)

# print(fitted_curve1.errors())

if SHOW_PLOTS:
    fitted_curve1.plot("GBP Yield Curve")

fit_type = CurveFitTypes.QUINTIC_POLYNOMIAL
fitted_curve2 = BondParametricYieldCurve(settle_dt, bonds, ylds, fit_type)
if SHOW_PLOTS:
    fitted_curve2.plot("GBP Yield Curve")

# print(fitted_curve2.errors())

fit_type = CurveFitTypes.NELSON_SIEGEL
fitted_curve3 = BondParametricYieldCurve(settle_dt, bonds, ylds, fit_type)
if SHOW_PLOTS:
    fitted_curve3.plot("GBP Yield Curve")

# print(fitted_curve3.errors())

fit_type = CurveFitTypes.NELSON_SIEGEL_SVENSSON
fitted_curve4 = BondParametricYieldCurve(settle_dt, bonds, ylds, fit_type)
if SHOW_PLOTS:
    fitted_curve4.plot("GBP Yield Curve")

# print(fitted_curve4.errors())

fit_type = CurveFitTypes.BSPLINE
fitted_curve5 = BondParametricYieldCurve(settle_dt, bonds, ylds, fit_type)

# print(fitted_curve5.errors())

if SHOW_PLOTS:
    fitted_curve5.plot("GBP Yield Curve")

print("PARAMETER", "VALUE")
print("values", fitted_curve1.curve_fit.coeffs)

print("PARAMETER", "VALUE")
print("values", fitted_curve2.curve_fit.coeffs)

print("PARAMETER", "VALUE")
print("beta_1", fitted_curve3.curve_fit.beta_1)
print("beta_2", fitted_curve3.curve_fit.beta_2)
print("beta_3", fitted_curve3.curve_fit.beta_3)
print("tau", fitted_curve3.curve_fit.tau)

print("PARAMETER", "VALUE")
print("beta_1", fitted_curve4.curve_fit.beta_1)
print("beta_2", fitted_curve4.curve_fit.beta_2)
print("beta_3", fitted_curve4.curve_fit.beta_3)
print("beta_4", fitted_curve4.curve_fit.beta_4)
print("tau_1", fitted_curve4.curve_fit.tau_1)
print("tau_2", fitted_curve4.curve_fit.tau_2)

maturity_dt = Date(19, 9, 2030)
interp_yield = fitted_curve5.interp_yield(maturity_dt)
print(maturity_dt, interp_yield)

