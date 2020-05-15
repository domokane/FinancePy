# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinDate import FinDate, fromDatetime
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes

#  from scipy import optimize
import datetime as dt

###############################################################################

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBondPortfolio():

    import pandas as pd
    bondDataFrame = pd.read_csv('./data/giltbondprices.txt', sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA

    settlement = FinDate(2012, 9, 19)

    testCases.header("DCTYPE", "MATDATE", "CPN", "PRICE", "ACCD", "YTM")

    for accrualType in FinDayCountTypes:

        for index, bond in bondDataFrame.iterrows():

            dateString = bond['maturity']
            matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
            maturityDt = fromDatetime(matDatetime)
            coupon = bond['coupon']/100.0
            cleanPrice = bond['mid']
            bond = FinBond(maturityDt, coupon, frequencyType, accrualType)

            ytm = bond.yieldToMaturity(settlement, cleanPrice)
            accd = bond._accruedInterest

            testCases.print(accrualType, maturityDt, coupon*100.0,
                            cleanPrice, accd, ytm*100.0)

##########################################################################


test_FinBondPortfolio()
testCases.compareTestCases()
