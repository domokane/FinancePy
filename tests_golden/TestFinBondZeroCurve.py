###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.zero_curve import BondZeroCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
import datetime as dt
import os
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

###############################################################################


def test_BondZeroCurve():

    import pandas as pd
    path = os.path.join(os.path.dirname(__file__), './data/giltBondPrices.txt')
    bondDataFrame = pd.read_csv(path, sep='\t')
    bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    settlement = Date(19, 9, 2012)

    bonds = []
    clean_prices = []

    for _, bondRow in bondDataFrame.iterrows():
        date_string = bondRow['maturity']
        matDatetime = dt.datetime.strptime(date_string, '%d-%b-%y')
        maturityDt = from_datetime(matDatetime)
        issueDt = Date(maturityDt._d, maturityDt._m, 2000)
        coupon = bondRow['coupon']/100.0
        clean_price = bondRow['mid']
        bond = Bond(issueDt, maturityDt, coupon, freq_type, accrual_type)
        bonds.append(bond)
        clean_prices.append(clean_price)

###############################################################################

    bondCurve = BondZeroCurve(settlement, bonds, clean_prices)

    testCases.header("DATE", "ZERO RATE")

    for _, bond in bondDataFrame.iterrows():

        date_string = bond['maturity']
        matDatetime = dt.datetime.strptime(date_string, '%d-%b-%y')
        maturityDt = from_datetime(matDatetime)
        zero_rate = bondCurve.zero_rate(maturityDt)
        testCases.print(maturityDt, zero_rate)

    if plotGraphs:
        bondCurve.plot("BOND CURVE")

###############################################################################


test_BondZeroCurve()
testCases.compareTestCases()
