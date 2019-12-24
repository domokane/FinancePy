# -*- coding: utf-8 -*-

import sys
sys.path.append("..//..")

import numpy as np
from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.finutils.FinHelperFunctions import printTree
from financepy.products.bonds.FinBondOption import FinBondOption, FinBondOptionTypes
from financepy.models.FinHullWhiteRateModel import FinHullWhiteRateModel
from financepy.models.FinBlackKarasinskiRateModel import FinBlackKarasinskiRateModel

import matplotlib.pyplot as plt
import time


###############################################################################P

def test_FinBondOption():

    settlementDate = FinDate(1, 12, 2019)

    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    tmat = (maturityDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dfs = np.exp(-0.05*times)
    discountCurve = FinDiscountCurve(settlementDate, times, dfs)

    expiryDate = settlementDate.addTenor("18m")
    strikePrice = 105.0
    face = 100.0
    optionType = FinBondOptionTypes.AMERICAN_CALL

    price = bond.fullPriceFromDiscountCurve(settlementDate, discountCurve)
    print("Fixed Income Price:", price)

    for strikePrice in [90,95,100,105,110]:

        sigma = 0.01
        a = 0.1

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinHullWhiteRateModel(a, sigma)
        v = bondOption.value(settlementDate, discountCurve, model)
        print("HW", strikePrice, v)

    for strikePrice in [90,95,100,105,110]:

        sigma = 0.20
        a = 0.05

        bondOption = FinBondOption(bond, expiryDate, strikePrice, face, optionType)
        model = FinBlackKarasinskiRateModel(a, sigma)
        v = bondOption.value(settlementDate, discountCurve, model)
        print("BK", strikePrice, v)

###############################################################################


test_FinBondOption()
