###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import datetime as dt

import sys
sys.path.append("..")

from financepy.utils.date import Date, fromDatetime
from financepy.products.bonds.bond import Bond
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_BondPortfolio():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA

    settlement = Date(19, 9, 2012)

    testCases.header("DCTYPE", "MATDATE", "CPN", "PRICE", "ACCD", "YTM")

    for accrual_type in DayCountTypes:

        for _, bond in bondDataFrame.iterrows():

            dateString = bond['maturity']
            matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
            maturityDt = fromDatetime(matDatetime)
            issueDt = Date(maturityDt._d, maturityDt._m, 2000)
            coupon = bond['coupon']/100.0
            clean_price = bond['mid']
            bond = Bond(issueDt, maturityDt,
                           coupon, freq_type, accrual_type)

            ytm = bond.yield_to_maturity(settlement, clean_price)
            accrued_interest= bond._accruedInterest

            testCases.print(accrual_type, maturityDt, coupon*100.0,
                            clean_price, accrued_interest, ytm*100.0)

##########################################################################


test_BondPortfolio()
testCases.compareTestCases()
