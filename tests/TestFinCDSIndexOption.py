###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import os
import time
import numpy as np

import sys
sys.path.append("..")

from financepy.products.credit.FinCDSIndexPortfolio import FinCDSIndexPortfolio
from financepy.products.credit.FinCDSIndexOption import FinCDSIndexOption
from financepy.products.credit.cds import FinCDS
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.products.credit.cds_curve import FinCDSCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


##########################################################################


def buildIborCurve(tradeDate):

    valuation_date = tradeDate.addDays(1)
    dcType = DayCountTypes.ACT_360
    depos = []

    depos = []
    fras = []
    swaps = []

    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    settlement_date = valuation_date

    maturity_date = settlement_date.addMonths(12)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(24)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(36)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(48)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(60)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve

##########################################################################


def buildFlatIssuerCurve(tradeDate, libor_curve, spread, recovery_rate):

    valuation_date = tradeDate.addDays(1)

    cdsMarketContracts = []

    maturity_date = Date(29, 6, 2010)
    cds = FinCDS(valuation_date, maturity_date, spread)
    cdsMarketContracts.append(cds)

    issuer_curve = FinCDSCurve(valuation_date,
                              cdsMarketContracts,
                              libor_curve,
                              recovery_rate)

    return issuer_curve

##########################################################################


def test_full_priceCDSIndexOption():

    tradeDate = Date(1, 8, 2007)
    step_in_date = tradeDate.addDays(1)
    valuation_date = step_in_date

    libor_curve = buildIborCurve(tradeDate)

    maturity3Y = tradeDate.nextCDSDate(36)
    maturity5Y = tradeDate.nextCDSDate(60)
    maturity7Y = tradeDate.nextCDSDate(84)
    maturity10Y = tradeDate.nextCDSDate(120)

    path = os.path.join(os.path.dirname(__file__), './/data//CDX_NA_IG_S7_SPREADS.csv')
    f = open(path, 'r')
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        creditName = splitRow[0]
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = FinCDS(step_in_date, maturity3Y, spd3Y)
        cds5Y = FinCDS(step_in_date, maturity5Y, spd5Y)
        cds7Y = FinCDS(step_in_date, maturity7Y, spd7Y)
        cds10Y = FinCDS(step_in_date, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = FinCDSCurve(valuation_date,
                                  cds_contracts,
                                  libor_curve,
                                  recovery_rate)

        issuer_curves.append(issuer_curve)

    ##########################################################################
    ##########################################################################

    indexUpfronts = [0.0, 0.0, 0.0, 0.0]
    indexMaturityDates = [Date(20, 12, 2009),
                          Date(20, 12, 2011),
                          Date(20, 12, 2013),
                          Date(20, 12, 2016)]
    indexRecovery = 0.40

    testCases.banner(
        "======================= CDS INDEX OPTION ==========================")

    indexCoupon = 0.004
    volatility = 0.50
    expiry_date = Date(1, 2, 2008)
    maturity_date = Date(20, 12, 2011)
    notional = 10000.0
    tolerance = 1e-6

    testCases.header(
        "TIME",
        "STRIKE",
        "INDEX",
        "PAY",
        "RECEIVER",
        "G(K)",
        "X",
        "EXPH",
        "ABPAY",
        "ABREC")

    for index in np.linspace(20, 60, 10):

        #######################################################################

        cds_contracts = []
        for dt in indexMaturityDates:
            cds = FinCDS(valuation_date, dt, index / 10000.0)
            cds_contracts.append(cds)

        index_curve = FinCDSCurve(valuation_date, cds_contracts,
                                 libor_curve, indexRecovery)

        if 1 == 1:

            indexSpreads = [index / 10000.0] * 4

            indexPortfolio = FinCDSIndexPortfolio()
            adjustedIssuerCurves = indexPortfolio.hazardRateAdjustIntrinsic(
                valuation_date,
                issuer_curves,
                indexSpreads,
                indexUpfronts,
                indexMaturityDates,
                indexRecovery,
                tolerance)
        else:

            indexSpread = index / 10000.0
            issuer_curve = buildFlatIssuerCurve(tradeDate,
                                               libor_curve,
                                               indexSpread,
                                               indexRecovery)

            adjustedIssuerCurves = []
            for iCredit in range(0, 125):
                adjustedIssuerCurves.append(issuer_curve)

        #######################################################################

        for strike in np.linspace(20, 60, 20):

            start = time.time()

            option = FinCDSIndexOption(expiry_date,
                                       maturity_date,
                                       indexCoupon,
                                       strike / 10000.0,
                                       notional)

            v_pay_1, v_rec_1, strikeValue, mu, expH = option.valueAnderson(
                valuation_date, adjustedIssuerCurves, indexRecovery, volatility)
            end = time.time()
            elapsed = end - start

            end = time.time()

            v_pay_2, v_rec_2 = option.valueAdjustedBlack(valuation_date,
                                                         index_curve,
                                                         indexRecovery,
                                                         libor_curve,
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


test_full_priceCDSIndexOption()
testCases.compareTestCases()
