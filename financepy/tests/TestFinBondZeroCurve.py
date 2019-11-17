# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""


import sys
import datetime as dt

from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond
from financepy.market.curves.FinBondZeroCurve import FinBondZeroCurve

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

sys.path.append("..//..")

##########################################################################
##########################################################################


def test_FinBondZeroCurve():

    import pandas as pd
    bondDataFrame = pd.read_csv('./data/giltbondprices.txt', sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    settlement = FinDate(2012, 9, 19)

    bonds = []
    cleanPrices = []

    for index, bond in bondDataFrame.iterrows():

        dateString = bond['maturity']
        matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
        maturityDt = FinDate.fromDatetime(matDatetime)
        coupon = bond['coupon']/100.0
        cleanPrice = bond['mid']
        bond = FinBond(maturityDt, coupon, frequencyType, accrualType)
        bonds.append(bond)
        cleanPrices.append(cleanPrice)

###############################################################################

    bondCurve = FinBondZeroCurve(settlement, bonds, cleanPrices)

    testCases.header("DATE", "ZERO RATE")

    for index, bond in bondDataFrame.iterrows():

        dateString = bond['maturity']
        matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
        maturityDt = FinDate.fromDatetime(matDatetime)
        zeroRate = bondCurve.zeroRate(maturityDt)
        testCases.print(maturityDt, zeroRate)

#    bondCurve.display("BOND CURVE")

###############################################################################


test_FinBondZeroCurve()
testCases.compareTestCases()
