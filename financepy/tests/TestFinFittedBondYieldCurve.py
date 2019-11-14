# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import datetime as dt

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond
from financepy.market.curves.FinFitBondYieldCurve import FinFitBondYieldCurve

from financepy.market.curves.FinCurveFitMethod import FinCurveFitMethodPolynomial
from financepy.market.curves.FinCurveFitMethod import FinCurveFitMethodNelsonSiegel
from financepy.market.curves.FinCurveFitMethod import FinCurveFitMethodNelsonSiegelSvensson
from financepy.market.curves.FinCurveFitMethod import FinCurveFitMethodBSpline

##########################################################################
# TODO: Inherit from FinCurve and add df method
# TODO: Add fitting optimiser to take in bonds and do a best fit
##########################################################################

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
##########################################################################


def test_FinFitBondYieldCurve():

    ###########################################################################

    import pandas as pd
    bondDataFrame = pd.read_csv('./data/giltbondprices.txt', sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    settlement = FinDate(2012, 9, 19)

    bonds = []
    ylds = []

    for index, bond in bondDataFrame.iterrows():

        dateString = bond['maturity']
        matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
        maturityDt = FinDate.fromDatetime(matDatetime)
        coupon = bond['coupon']/100.0
        cleanPrice = bond['mid']
        bond = FinBond(maturityDt, coupon, frequencyType, accrualType)
        yld = bond.yieldToMaturity(settlement, cleanPrice)
        bonds.append(bond)
        ylds.append(yld)

    ###########################################################################

    curveFitMethod = FinCurveFitMethodPolynomial()
    fittedCurve1 = FinFitBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
    fittedCurve1.display("GBP Yield Curve")

    curveFitMethod = FinCurveFitMethodPolynomial(5)
    fittedCurve2 = FinFitBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
    fittedCurve2.display("GBP Yield Curve")

    curveFitMethod = FinCurveFitMethodNelsonSiegel()
    fittedCurve3 = FinFitBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
    fittedCurve3.display("GBP Yield Curve")

    curveFitMethod = FinCurveFitMethodNelsonSiegelSvensson()
    fittedCurve4 = FinFitBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
    fittedCurve4.display("GBP Yield Curve")

    curveFitMethod = FinCurveFitMethodBSpline()
    fittedCurve5 = FinFitBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
    fittedCurve5.display("GBP Yield Curve")

    ###########################################################################

    print("beta1", fittedCurve3._curveFitMethod._beta1)
    print("beta2", fittedCurve3._curveFitMethod._beta2)
    print("beta3", fittedCurve3._curveFitMethod._beta3)
    print("tau", fittedCurve3._curveFitMethod._tau)

    print("beta1", fittedCurve4._curveFitMethod._beta1)
    print("beta2", fittedCurve4._curveFitMethod._beta2)
    print("beta3", fittedCurve4._curveFitMethod._beta3)
    print("beta4", fittedCurve4._curveFitMethod._beta4)
    print("tau1", fittedCurve4._curveFitMethod._tau1)
    print("tau2", fittedCurve4._curveFitMethod._tau2)

    ###########################################################################


    maturityDate = FinDate(2030, 9, 19)
    interpolatedYield = fittedCurve5.interpolatedYield(maturityDate)
    print(maturityDate, interpolatedYield)

    ###########################################################################

test_FinFitBondYieldCurve()
testCases.compareTestCases()
