###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import datetime as dt

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate, fromDatetime
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinBondPortfolio():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freqType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA

    settlement = FinDate(19, 9, 2012)

    testCases.header("DCTYPE", "MATDATE", "CPN", "PRICE", "ACCD", "YTM")

    for accrualType in FinDayCountTypes:

        for _, bond in bondDataFrame.iterrows():

            dateString = bond['maturity']
            matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
            maturityDt = fromDatetime(matDatetime)
            issueDt = FinDate(maturityDt._d, maturityDt._m, 2000)
            coupon = bond['coupon']/100.0
            cleanPrice = bond['mid']
            bond = FinBond(issueDt, maturityDt, 
                           coupon, freqType, accrualType)

            ytm = bond.yieldToMaturity(settlement, cleanPrice)
            accd = bond._accruedInterest

            testCases.print(accrualType, maturityDt, coupon*100.0,
                            cleanPrice, accd, ytm*100.0)

##########################################################################


test_FinBondPortfolio()
testCases.compareTestCases()
