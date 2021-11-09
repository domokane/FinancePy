###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.bonds.zero_curve import BondZeroCurve
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
import datetime as dt
import os


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

bondCurve = BondZeroCurve(settlement, bonds, clean_prices)


def test_zero_curve():

    maturityDt = Date(7, 3, 2013)
    zero_rate = bondCurve.zero_rate(maturityDt)
    assert round(zero_rate, 4) == 0.0022

    maturityDt = Date(7, 9, 2019)
    zero_rate = bondCurve.zero_rate(maturityDt)
    assert round(zero_rate, 4) == 0.0128

    maturityDt = Date(22, 1, 2060)
    zero_rate = bondCurve.zero_rate(maturityDt)
    assert round(zero_rate, 4) == 0.0351
