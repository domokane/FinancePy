# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from math import sqrt
import numpy as np

from financepy.finutils.FinDate import FinDate
from financepy.market.volatility.FinVolatilityCurve import FinVolatilityCurve
from financepy.products.equities.FinVarianceSwap import FinVarianceSwap
import sys
sys.path.append("..//..")


def volSkew(K, atmVol, atmK, skew):
    v = atmVol + skew * (K-atmK)
    return v


startDate = FinDate(2018, 3, 20)
tenor = "3M"
strike = 0.3*0.3

volSwap = FinVarianceSwap(startDate, tenor, strike)

valuationDate = FinDate(2018, 3, 20)
stockPrice = 100.0
dividendYield = 0.0

atmVol = 0.20
atmK = 100.0
skew = -0.02/5.0  # defined as dsigma/dK
strikes = np.linspace(50.0, 135.0, 18)
volatilities = volSkew(strikes, atmVol, atmK, skew)
volCurve = FinVolatilityCurve(valuationDate, strikes, volatilities)

print(strikes)
print(volatilities)

strikeSpacing = 5.0
numCallOptions = 10
numPutOptions = 10
r = 0.05
useForward = False

k1 = volSwap.fairStrike(valuationDate, stockPrice, dividendYield,
                        volCurve, numCallOptions, numPutOptions,
                        strikeSpacing, r, useForward)
print("REPLICATION VARIANCE:", k1)

volSwap.print()

k2 = volSwap.fairStrikeApprox(valuationDate, stockPrice, strikes, volatilities)
print("DERMAN SKEW APPROX for K:", k2)
