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
from financepy.products.bonds.yield_curve import BondYieldCurve
from financepy.products.bonds.yield_curve_model import CurveFitPolynomial
from financepy.products.bonds.yield_curve_model import CurveFitBSpline
from financepy.products.bonds.yield_curve_model import CurveFitNelsonSiegel
from financepy.products.bonds.yield_curve_model import CurveFitNelsonSiegelSvensson

path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
bondDataFrame = pd.read_csv(path, sep='\t')
bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

freq_type = FrequencyTypes.SEMI_ANNUAL
accrual_type = DayCountTypes.ACT_ACT_ICMA
settlement = Date(19, 9, 2012)

bonds = []
ylds = []

for _, bond in bondDataFrame.iterrows():

    date_string = bond['maturity']
    matDatetime = dt.datetime.strptime(date_string, '%d-%b-%y')
    maturityDt = from_datetime(matDatetime)
    issueDt = Date(maturityDt._d, maturityDt._m, 2000)
    coupon = bond['coupon']/100.0
    clean_price = bond['mid']
    bond = Bond(issueDt, maturityDt, coupon, freq_type, accrual_type)
    yld = bond.yield_to_maturity(settlement, clean_price)
    bonds.append(bond)
    ylds.append(yld)


def test_poly():
    curveFitMethod = CurveFitPolynomial(5)
    fittedCurve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    coeffs = fittedCurve._curveFit._coeffs
    assert round(coeffs[0] * 1e9, 4) == -1.4477
    assert round(coeffs[1] * 1e7, 4) == 1.7840
    assert round(coeffs[2] * 1e6, 4) == -7.4147
    assert round(coeffs[3] * 1e5, 4) == 9.0622
    assert round(coeffs[4] * 1e3, 4) == 1.3536
    assert round(coeffs[5] * 1e7, 4) == 4.1514


def test_nelson_siegel():
    curveFitMethod = CurveFitNelsonSiegel()
    fittedCurve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    assert round(fittedCurve._curveFit._beta1, 3) == -0.094
    assert round(fittedCurve._curveFit._beta2, 3) == 0.092
    assert round(fittedCurve._curveFit._beta3, 3) == 0.259
    assert round(fittedCurve._curveFit._tau, 1) == 35.8


def test_nelson_siegel_svensson():
    curveFitMethod = CurveFitNelsonSiegelSvensson()
    fittedCurve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    assert round(fittedCurve._curveFit._beta1, 4) == 0.0460
    assert round(fittedCurve._curveFit._beta2, 4) == -0.0433
    assert round(fittedCurve._curveFit._beta3, 4) == -0.0523
    assert round(fittedCurve._curveFit._beta4, 4) == -0.0376
    assert round(fittedCurve._curveFit._tau1, 3) == 3.177
    assert round(fittedCurve._curveFit._tau2, 4) == 100.0000


def test_interpolated_yield():
    curveFitMethod = CurveFitBSpline()
    fittedCurve = BondYieldCurve(settlement, bonds, ylds, curveFitMethod)

    maturity_date = Date(19, 9, 2030)
    interpolated_yield = fittedCurve.interpolated_yield(maturity_date)
    assert round(float(interpolated_yield), 8) == 0.02601858
