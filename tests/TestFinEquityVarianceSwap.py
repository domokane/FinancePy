###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.market.volatility.FinEquityVolCurve import FinEquityVolCurve
from financepy.products.equity.FinEquityVarianceSwap import FinEquityVarianceSwap
from financepy.market.discount.curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def volSkew(K, atmVol, atmK, skew):
    v = atmVol + skew * (K-atmK)
    return v

###############################################################################


def test_FinEquityVarianceSwap():

    start_date = Date(20, 3, 2018)
    tenor = "3M"
    strike = 0.3*0.3

    volSwap = FinEquityVarianceSwap(start_date, tenor, strike)

    valuation_date = Date(20, 3, 2018)
    stock_price = 100.0
    dividendYield = 0.0
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    maturity_date = start_date.addMonths(3)

    atmVol = 0.20
    atmK = 100.0
    skew = -0.02/5.0  # defined as dsigma/dK
    strikes = np.linspace(50.0, 135.0, 18)
    vols = volSkew(strikes, atmVol, atmK, skew)
    volCurve = FinEquityVolCurve(valuation_date, maturity_date, strikes, vols)

    strikeSpacing = 5.0
    numCallOptions = 10
    numPutOptions = 10
    r = 0.05

    discount_curve = DiscountCurveFlat(valuation_date, r)

    useForward = False

    testCases.header("LABEL", "VALUE")

    k1 = volSwap.fairStrike(valuation_date, stock_price, dividendCurve,
                            volCurve, numCallOptions, numPutOptions,
                            strikeSpacing, discount_curve, useForward)

    testCases.print("REPLICATION VARIANCE:", k1)

    k2 = volSwap.fairStrikeApprox(valuation_date, stock_price, strikes, vols)
    testCases.print("DERMAN SKEW APPROX for K:", k2)

##########################################################################


test_FinEquityVarianceSwap()
testCases.compareTestCases()
