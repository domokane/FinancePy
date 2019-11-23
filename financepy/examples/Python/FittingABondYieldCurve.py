###########################################################################

import datetime as dt
import pandas as pd

from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond
from financepy.market.curves.FinBondYieldCurve import FinBondYieldCurve
from financepy.market.curves.FinBondYieldCurveModel import *

###############################################################################
# FIRST LOAD UP SOME UK BOND DATA USING PANDAS

bondDataFrame = pd.read_csv('../tests/data/giltbondprices.txt', sep='\t')

# CALCULATE MID MARKET PRICES
bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

# SPECIFY UK BOND CONVENTIONS
frequencyType = FinFrequencyTypes.SEMI_ANNUAL
accrualType = FinDayCountTypes.ACT_ACT_ICMA
settlement = FinDate(2012, 9, 19)

bonds = []
ylds = []

# LOAD BONDS AND CREATE A VECTOR OF FINBOND AND THEIR CORRESPONDING YIELDS

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

# FIT THE BOND YIELDS TO A CUBIC POLYNOMIAL
curveFitMethod = FinCurveFitMethodPolynomial()
fittedCurve1 = FinBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
fittedCurve1.display("GBP Yield Curve")

# FIT THE BOND YIELDS TO A QUINTIC POLYNOMIAL
curveFitMethod = FinCurveFitMethodPolynomial(5)
fittedCurve2 = FinBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
fittedCurve2.display("GBP Yield Curve")

# FIT THE BONDS TO A NELSON-SIEGEL CURVE
curveFitMethod = FinCurveFitMethodNelsonSiegel()
fittedCurve3 = FinBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
fittedCurve3.display("GBP Yield Curve")

print("NELSON-SIEGEL PARAMETERS")
print("beta1", fittedCurve3._curveFitMethod._beta1)
print("beta2", fittedCurve3._curveFitMethod._beta2)
print("beta3", fittedCurve3._curveFitMethod._beta3)
print("tau", fittedCurve3._curveFitMethod._tau)

# FIT THE BONDS TO A NELSON-SIEGEL-SVENSSON CURVE
curveFitMethod = FinCurveFitMethodNelsonSiegelSvensson()
fittedCurve4 = FinBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
fittedCurve4.display("GBP Yield Curve")

print("NELSON-SIEGEL-SVENSSON PARAMETERS")
print("beta1", fittedCurve4._curveFitMethod._beta1)
print("beta2", fittedCurve4._curveFitMethod._beta2)
print("beta3", fittedCurve4._curveFitMethod._beta3)
print("beta4", fittedCurve4._curveFitMethod._beta4)
print("tau1", fittedCurve4._curveFitMethod._tau1)
print("tau2", fittedCurve4._curveFitMethod._tau2)

# FIT THE BONDS TO A B-SPLINE CURVE
curveFitMethod = FinCurveFitMethodBSpline()
fittedCurve5 = FinBondYieldCurve(settlement, bonds, ylds, curveFitMethod)
fittedCurve5.display("GBP Yield Curve")

# EXTRACT A YIELD FROM A FITTED YIELD CURVE
maturityDate = FinDate(2030, 9, 19)
interpolatedYield = fittedCurve5.interpolatedYield(maturityDate)
print(maturityDate, interpolatedYield)
