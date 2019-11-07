# -*- coding: utf-8 -*-

# TODO
import time

from financepy.products.bonds.FinConvertibleBond import FinConvertibleBond
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.market.curves.FinFlatCurve import FinCompoundingMethods
from financepy.products.bonds.FinBond import FinBondAccruedTypes


def test_FinConvertibleBond():

    settlementDate = FinDate(2003, 12, 31)
    startConvertDate = FinDate(2003, 12, 31)
    maturityDate = FinDate(2022, 3, 15)
    conversionRatio = 38.4615  # adjust for face
    coupon = 0.0575
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualBasis = FinBondAccruedTypes.ACT_365
    face = 1000.0

    callPrice = 1100
    callDates = [FinDate(2007, 3, 20),
                 FinDate(2012, 3, 15),
                 FinDate(2017, 3, 15)]
    callPrices = [callPrice, callPrice, callPrice]

    putPrice = 90
    putDates = [FinDate(2007, 3, 20),
                FinDate(2012, 3, 15),
                FinDate(2017, 3, 15)]
    putPrices = [putPrice, putPrice, putPrice]

    bond = FinConvertibleBond(maturityDate,
                              coupon,
                              frequencyType,
                              startConvertDate,
                              conversionRatio,
                              callDates,
                              callPrices,
                              putDates,
                              putPrices,
                              accrualBasis,
                              face)

    dividendDates = [FinDate(2007, 3, 20),
                     FinDate(2008, 3, 15),
                     FinDate(2009, 3, 15),
                     FinDate(2010, 3, 15),
                     FinDate(2011, 3, 15),
                     FinDate(2012, 3, 15),
                     FinDate(2013, 3, 15),
                     FinDate(2014, 3, 15),
                     FinDate(2015, 3, 15),
                     FinDate(2016, 3, 15),
                     FinDate(2017, 3, 15),
                     FinDate(2018, 3, 15),
                     FinDate(2019, 3, 15),
                     FinDate(2020, 3, 15),
                     FinDate(2021, 3, 15),
                     FinDate(2022, 3, 15)]

    dividendYields = [0.00] * 16
    stockPrice = 28.5
    stockVolatility = 0.370
    rate = 0.04
    discountCurve = FinFlatCurve(settlementDate,
                                 rate,
                                 FinCompoundingMethods.CONTINUOUS)
    creditSpread = 0.00
    recoveryRate = 0.40
    numStepsPerYear = 20

    print("NO CALLS OR PUTS")
    for numStepsPerYear in [1, 2, 5, 10, 20, 40, 80]:
        start = time.time()
        res = bond.value(settlementDate,
                         stockPrice,
                         stockVolatility,
                         dividendDates,
                         dividendYields,
                         discountCurve,
                         creditSpread,
                         recoveryRate,
                         numStepsPerYear)

        end = time.time()
        period = end - start
        print(period, numStepsPerYear, "PRICE:", res[0], "BULLET", res[1], "DTA", res[2])

    dividendYields = [0.02] * 16
    print("DIVIDENDS")
    for numStepsPerYear in [1, 2, 5, 10, 20, 40, 80]:
        start = time.time()
        res = bond.value(settlementDate,
                         stockPrice,
                         stockVolatility,
                         dividendDates,
                         dividendYields,
                         discountCurve,
                         creditSpread,
                         recoveryRate,
                         numStepsPerYear)
        end = time.time()
        period = end - start
        print(period, numStepsPerYear, "PRICE:", res[0], "BULLET", res[1], "DTA", res[2])


test_FinConvertibleBond()
