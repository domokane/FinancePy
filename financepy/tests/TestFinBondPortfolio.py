# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond, FinBondAccruedTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes

#  from scipy import optimize
import datetime as dt

###############################################################################

import pandas as pd
bondDataFrame = pd.read_csv('./data/giltbondprices.txt', sep='\t')
bondDataFrame['mid'] = 0.5*(bondDataFrame['bid'] + bondDataFrame['ask'])

frequencyType = FinFrequencyTypes.SEMI_ANNUAL
accrualType = FinBondAccruedTypes.ACT_365

ytms = []
yearsToMaturities = []
fullPrices = []

settlement = FinDate(2012, 9, 19)

for index, bond in bondDataFrame.iterrows():

    dateString = bond['maturity']
    matDatetime = dt.datetime.strptime(dateString, '%d-%b-%y')
    maturityDt = FinDate.fromDatetime(matDatetime)
    coupon = bond['coupon']/100.0
    cleanPrice = bond['mid']
    bond = FinBond(maturityDt, coupon, frequencyType, accrualType)

    ytm = bond.yieldToMaturity(settlement, cleanPrice)
    accd = bond.accruedInterest(settlement)
    accd_days = bond.accruedDays(settlement)

    print("%5d %18s %8.4f %10.4f %6.0f %10.4f %8.4f" %
          (index, maturityDt, coupon*100.0,
           cleanPrice, accd_days, accd, ytm*100.0))
