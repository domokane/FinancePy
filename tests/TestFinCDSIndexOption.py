# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

import time
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.credit.FinCDSIndexPortfolio import FinCDSIndexPortfolio
from financepy.products.credit.FinCDSIndexOption import FinCDSIndexOption
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.market.curves.FinCDSCurve import FinCDSCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
import numpy as np
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################


def buildLiborCurve(tradeDate):

    valuationDate = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360
    depos = []

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlementDate = valuationDate

    maturityDate = settlementDate.addMonths(12)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(24)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(36)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(48)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(60)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    liborCurve = FinLiborCurve(
        "USD_LIBOR", settlementDate, depos, fras, swaps)

    return liborCurve

##########################################################################


def buildFlatIssuerCurve(tradeDate, liborCurve, spread, recoveryRate):

    valuationDate = tradeDate.addDays(1)

    cdsMarketContracts = []

    maturityDate = FinDate(2010, 6, 29)
    cds = FinCDS(valuationDate, maturityDate, spread)
    cdsMarketContracts.append(cds)

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsMarketContracts,
                              liborCurve,
                              recoveryRate)

    return issuerCurve

##########################################################################


def test_fullPriceCDSIndexOption():

    tradeDate = FinDate(2007, 8, 1)
    stepInDate = tradeDate.addDays(1)
    valuationDate = stepInDate

    liborCurve = buildLiborCurve(tradeDate)

    maturity3Y = tradeDate.nextCDSDate(36)
    maturity5Y = tradeDate.nextCDSDate(60)
    maturity7Y = tradeDate.nextCDSDate(84)
    maturity10Y = tradeDate.nextCDSDate(120)

    f = open('.//data//CDX_NA_IG_S7_SPREADS.csv', 'r')
    data = f.readlines()
    issuerCurves = []

    for row in data[1:]:

        splitRow = row.split(",")
        creditName = splitRow[0]
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recoveryRate = float(splitRow[5])

        cds3Y = FinCDS(stepInDate, maturity3Y, spd3Y)
        cds5Y = FinCDS(stepInDate, maturity5Y, spd5Y)
        cds7Y = FinCDS(stepInDate, maturity7Y, spd7Y)
        cds10Y = FinCDS(stepInDate, maturity10Y, spd10Y)
        cdsContracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuerCurve = FinCDSCurve(valuationDate,
                                  cdsContracts,
                                  liborCurve,
                                  recoveryRate)

        issuerCurves.append(issuerCurve)

    ##########################################################################
    ##########################################################################

    indexUpfronts = [0.0, 0.0, 0.0, 0.0]
    indexMaturityDates = [FinDate(2009, 12, 20),
                          FinDate(2011, 12, 20),
                          FinDate(2013, 12, 20),
                          FinDate(2016, 12, 20)]
    indexRecovery = 0.40

    testCases.banner(
        "======================= CDS INDEX OPTION ==========================")

    indexCoupon = 0.004
    volatility = 0.50
    numCredits = 125
    expiryDate = FinDate(2008, 2, 1)
    maturityDate = FinDate(2011, 12, 20)
    notional = 10000.0
    tolerance = 1e-3

    testCases.header(
        "TIME",
        "STRIKE",
        "INDEX",
        "PAYER",
        "RECEIVER",
        "G(K)",
        "X",
        "EXPH",
        "ABPAY",
        "ABREC")

    for index in np.linspace(20, 60, 5):

        #######################################################################

        cdsContracts = []
        for dt in indexMaturityDates:
            cds = FinCDS(valuationDate, dt, index / 10000.0)
            cdsContracts.append(cds)

        indexCurve = FinCDSCurve(
            valuationDate,
            cdsContracts,
            liborCurve,
            indexRecovery)

        if 1 == 0:

            indexSpreads = [index / 10000.0] * 4

            adjustedIssuerCurves = FinCDSIndexPortfolio.hazardRateAdjustIntrinsic(
                valuationDate,
                issuerCurves,
                indexSpreads,
                indexUpfronts,
                indexMaturityDates,
                indexRecovery,
                tolerance)
        else:

            indexSpread = index / 10000.0
            issuerCurve = buildFlatIssuerCurve(tradeDate,
                                               liborCurve,
                                               indexSpread,
                                               indexRecovery)

            adjustedIssuerCurves = []
            for iCredit in range(0, 125):
                adjustedIssuerCurves.append(issuerCurve)

        #######################################################################

        for strike in np.linspace(20, 60, 5):

            start = time.time()

            option = FinCDSIndexOption(expiryDate,
                                       maturityDate,
                                       indexCoupon,
                                       strike / 10000.0,
                                       notional)

            v_pay_1, v_rec_1, strikeValue, mu, expH = option.valueAnderson(
                valuationDate, adjustedIssuerCurves, indexRecovery, volatility)
            end = time.time()
            elapsed = end - start

            end = time.time()

            v_pay_2, v_rec_2 = option.valueAdjustedBlack(valuationDate,
                                                         indexCurve,
                                                         indexRecovery,
                                                         liborCurve,
                                                         volatility)

            elapsed = end - start

            testCases.print(
                elapsed,
                strike,
                index,
                v_pay_1,
                v_rec_1,
                strikeValue,
                mu,
                expH,
                v_pay_2,
                v_rec_2)

##########################################################################


test_fullPriceCDSIndexOption()
testCases.compareTestCases()
