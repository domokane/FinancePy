###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import os
import datetime as dt

from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate, fromDatetime
from financepy.products.bonds.FinBond import FinBond
from financepy.products.bonds.FinBondZeroCurve import FinBondZeroCurve

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_FinBondZeroCurve():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    settlement = FinDate(19, 9, 2012)

    bonds = []
    cleanPrices = []

    for _, bondRow in bondDataFrame.iterrows():
        dateString = bondRow['maturity']
        matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
        maturityDt = fromDatetime(matDatetime)
        issueDt = FinDate(maturityDt._d, maturityDt._m, 2000)
        coupon = bondRow['coupon']/100.0
        cleanPrice = bondRow['mid']
        bond = FinBond(issueDt, maturityDt, coupon, freqType, accrualType)
        bonds.append(bond)
        cleanPrices.append(cleanPrice)

###############################################################################

    bondCurve = FinBondZeroCurve(settlement, bonds, cleanPrices)

    testCases.header("DATE", "ZERO RATE")

    for _, bond in bondDataFrame.iterrows():

        dateString = bond['maturity']
        matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
        maturityDt = fromDatetime(matDatetime)
        zeroRate = bondCurve.zeroRate(maturityDt)
        testCases.print(maturityDt, zeroRate)

    if plotGraphs:
        bondCurve.plot("BOND CURVE")

###############################################################################


test_FinBondZeroCurve()
testCases.compareTestCases()
