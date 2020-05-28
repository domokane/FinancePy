# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.market.volatility.FinVolatilityCurve import FinVolatilityCurve
from financepy.products.equity.FinEquityVarianceSwap import FinEquityVarianceSwap
from financepy.market.curves.FinFlatCurve import FinFlatCurve

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def volSkew(K, atmVol, atmK, skew):
    v = atmVol + skew * (K-atmK)
    return v

###############################################################################


def test_FinEquityVarianceSwap():

    startDate = FinDate(2018, 3, 20)
    tenor = "3M"
    strike = 0.3*0.3

    volSwap = FinEquityVarianceSwap(startDate, tenor, strike)

    valuationDate = FinDate(2018, 3, 20)
    stockPrice = 100.0
    dividendYield = 0.0
    maturityDate = startDate.addMonths(3)

    atmVol = 0.20
    atmK = 100.0
    skew = -0.02/5.0  # defined as dsigma/dK
    strikes = np.linspace(50.0, 135.0, 18)
    vols = volSkew(strikes, atmVol, atmK, skew)
    volCurve = FinVolatilityCurve(valuationDate, maturityDate, strikes, vols)

    strikeSpacing = 5.0
    numCallOptions = 10
    numPutOptions = 10
    r = 0.05
    discountCurve = FinFlatCurve(valuationDate, r)

    useForward = False

    testCases.header("LABEL", "VALUE")

    k1 = volSwap.fairStrike(valuationDate, stockPrice, dividendYield,
                            volCurve, numCallOptions, numPutOptions,
                            strikeSpacing, discountCurve, useForward)

    testCases.print("REPLICATION VARIANCE:", k1)

    # volSwap.print()

    k2 = volSwap.fairStrikeApprox(valuationDate, stockPrice, strikes, vols)
    testCases.print("DERMAN SKEW APPROX for K:", k2)

##########################################################################


test_FinEquityVarianceSwap()
testCases.compareTestCases()
