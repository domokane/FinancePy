###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import datetime as dt
import os

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.yield_curve_model import CurveFitNelsonSiegelSvensson
from financepy.products.bonds.yield_curve_model import CurveFitNelsonSiegel
from financepy.products.bonds.yield_curve_model import CurveFitBSpline
from financepy.products.bonds.yield_curve_model import CurveFitPolynomial
from financepy.products.bonds.yield_curve import BondYieldCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes



test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
###############################################################################


def test_BondYieldCurve():

    ###########################################################################

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    settlement = Date(19, 9, 2012)

    bonds = []
    ylds = []

    for _, bond in bondDataFrame.iterrows():

        date_string = bond['maturity']
        mat_date_time = dt.datetime.strptime(date_string, '%d-%b-%y')
        maturity_dt = from_datetime(mat_date_time)
        issue_dt = Date(maturity_dt._d, maturity_dt._m, 2000)
        coupon = bond['coupon']/100.0
        clean_price = bond['mid']
        bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)
        yld = bond.yield_to_maturity(settlement, clean_price)
        bonds.append(bond)
        ylds.append(yld)

###############################################################################

    curveFitMethod = CurveFitPolynomial()
    fitted_curve1 = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)
#    fitted_curve1.display("GBP Yield Curve")

    curveFitMethod = CurveFitPolynomial(5)
    fitted_curve2 = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)
#    fitted_curve2.display("GBP Yield Curve")

    curveFitMethod = CurveFitNelsonSiegel()
    fitted_curve3 = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)
#    fitted_curve3.display("GBP Yield Curve")

    curveFitMethod = CurveFitNelsonSiegelSvensson()
    fitted_curve4 = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)
#    fitted_curve4.display("GBP Yield Curve")

    curveFitMethod = CurveFitBSpline()
    fitted_curve5 = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)
#    fitted_curve5.display("GBP Yield Curve")

###############################################################################

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("values", fitted_curve1._curveFit._coeffs)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("values", fitted_curve2._curveFit._coeffs)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("beta1", fitted_curve3._curveFit._beta1)
    test_cases.print("beta2", fitted_curve3._curveFit._beta2)
    test_cases.print("beta3", fitted_curve3._curveFit._beta3)
    test_cases.print("tau", fitted_curve3._curveFit._tau)

    test_cases.header("PARAMETER", "VALUE")
    test_cases.print("beta1", fitted_curve4._curveFit._beta1)
    test_cases.print("beta2", fitted_curve4._curveFit._beta2)
    test_cases.print("beta3", fitted_curve4._curveFit._beta3)
    test_cases.print("beta4", fitted_curve4._curveFit._beta4)
    test_cases.print("tau1", fitted_curve4._curveFit._tau1)
    test_cases.print("tau2", fitted_curve4._curveFit._tau2)

###############################################################################

    maturity_dt = Date(19, 9, 2030)
    interpolated_yield = fitted_curve5.interpolated_yield(maturity_dt)
    test_cases.print(maturity_dt, interpolated_yield)

###############################################################################


test_BondYieldCurve()
test_cases.compareTestCases()
