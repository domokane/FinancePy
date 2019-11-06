# -*- coding: utf-8 -*-

# TODO
from math import log

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

    callDates = [FinDate(2007, 3, 20)]
    callPrices = [11000.0]

    putDates = [FinDate(2007, 3, 20),
                FinDate(2012, 3, 15),
                FinDate(2017, 3, 15)]
    putPrices = [95.0, 95.0, 95.0]

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

    dividendYields = [0.01] * 16
    stockPrice = 28.5
    stockVolatility = 0.370
    rate = 0.04
    discountCurve = FinFlatCurve(settlementDate, rate,
                                 FinCompoundingMethods.CONTINUOUS)
    creditSpread = 0.01
    numStepsPerYear = 200

    res = bond.value(settlementDate,
                     stockPrice,
                     stockVolatility,
                     dividendDates,
                     dividendYields,
                     discountCurve,
                     creditSpread,
                     numStepsPerYear)

    print("PRICE:", res[0])
    print("BULLET:", res[1])
    print("DELTA:", res[2])


test_FinConvertibleBond()
