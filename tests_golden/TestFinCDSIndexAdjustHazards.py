###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from os.path import dirname, join

import sys
sys.path.append("..")

from financepy.products.credit.FinCDSIndexPortfolio import FinCDSIndexPortfolio
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################


def buildIborCurve(tradeDate):

    valuationDate = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlementDate = valuationDate

    maturityDate = settlementDate.addMonths(12)
    swap1 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(24)
    swap2 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(36)
    swap3 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(48)
    swap4 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(60)
    swap5 = FinIborSwap(
        settlementDate,
        maturityDate,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    liborCurve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    return liborCurve

##########################################################################


def buildIssuerCurve(tradeDate, liborCurve):

    valuationDate = tradeDate.addDays(1)

    cdsMarketContracts = []

    cdsCoupon = 0.0048375
    maturityDate = FinDate(29, 6, 2010)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    recoveryRate = 0.40

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsMarketContracts,
                              liborCurve,
                              recoveryRate)

    return issuerCurve

##########################################################################


def test_performCDSIndexHazardRateAdjustment():

    tradeDate = FinDate(1, 8, 2007)
    stepInDate = tradeDate.addDays(1)
    valuationDate = stepInDate

    liborCurve = buildIborCurve(tradeDate)

    maturity3Y = tradeDate.nextCDSDate(36)
    maturity5Y = tradeDate.nextCDSDate(60)
    maturity7Y = tradeDate.nextCDSDate(84)
    maturity10Y = tradeDate.nextCDSDate(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, 'r')

    data = f.readlines()
    issuerCurves = []

    for row in data[1:]:

        splitRow = row.split(",")
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
    # Now determine the average spread of the index
    ##########################################################################

    cdsIndex = FinCDSIndexPortfolio()

    averageSpd3Y = cdsIndex.averageSpread(valuationDate,
                                          stepInDate,
                                          maturity3Y,
                                          issuerCurves) * 10000.0

    averageSpd5Y = cdsIndex.averageSpread(valuationDate,
                                          stepInDate,
                                          maturity5Y,
                                          issuerCurves) * 10000.0

    averageSpd7Y = cdsIndex.averageSpread(valuationDate,
                                          stepInDate,
                                          maturity7Y,
                                          issuerCurves) * 10000.0

    averageSpd10Y = cdsIndex.averageSpread(valuationDate,
                                           stepInDate,
                                           maturity10Y,
                                           issuerCurves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("AVERAGE SPD 3Y", averageSpd3Y)
    testCases.print("AVERAGE SPD 5Y", averageSpd5Y)
    testCases.print("AVERAGE SPD 7Y", averageSpd7Y)
    testCases.print("AVERAGE SPD 10Y", averageSpd10Y)
    testCases.banner(
        "===================================================================")

    ##########################################################################
    # Now determine the intrinsic spread of the index to the same maturity dates
    # As the single name CDS contracts
    ##########################################################################

    cdsIndex = FinCDSIndexPortfolio()

    intrinsicSpd3Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              maturity3Y,
                                              issuerCurves) * 10000.0

    intrinsicSpd5Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              maturity5Y,
                                              issuerCurves) * 10000.0

    intrinsicSpd7Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              maturity7Y,
                                              issuerCurves) * 10000.0

    intrinsicSpd10Y = cdsIndex.intrinsicSpread(valuationDate,
                                               stepInDate,
                                               maturity10Y,
                                               issuerCurves) * 10000.0

    ##########################################################################
    ##########################################################################

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD 3Y", intrinsicSpd3Y)
    testCases.print("INTRINSIC SPD 5Y", intrinsicSpd5Y)
    testCases.print("INTRINSIC SPD 7Y", intrinsicSpd7Y)
    testCases.print("INTRINSIC SPD 10Y", intrinsicSpd10Y)
    testCases.banner(
        "===================================================================")

    ##########################################################################
    ##########################################################################

    indexCoupons = [0.002, 0.0037, 0.0050, 0.0063]
    indexUpfronts = [0.0, 0.0, 0.0, 0.0]

    indexMaturityDates = [FinDate(20, 12, 2009),
                          FinDate(20, 12, 2011),
                          FinDate(20, 12, 2013),
                          FinDate(20, 12, 2016)]

    indexRecoveryRate = 0.40

    tolerance = 1e-6

    import time
    start = time.time()

    indexPortfolio = FinCDSIndexPortfolio()
    adjustedIssuerCurves = indexPortfolio.hazardRateAdjustIntrinsic(
        valuationDate,
        issuerCurves,
        indexCoupons,
        indexUpfronts,
        indexMaturityDates,
        indexRecoveryRate,
        tolerance)

    end = time.time()
    testCases.header("TIME")
    testCases.print(end - start)

#    numCredits = len(issuerCurves)
#    testCases.print("#","MATURITY","CDS_UNADJ","CDS_ADJ")
#    for m in range(0,numCredits):
#        for cds in cdsContracts:
#            unadjustedSpread = cds.parSpread(valuationDate,issuerCurves[m])
#            adjustedSpread = cds.parSpread(valuationDate,adjustedIssuerCurves[m])
#            testCases.print(m,str(cds._maturityDate),"%10.3f"%(unadjustedSpread*10000),"%10.3f" %(adjustedSpread*10000))

    cdsIndex = FinCDSIndexPortfolio()

    intrinsicSpd3Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              indexMaturityDates[0],
                                              adjustedIssuerCurves) * 10000.0

    intrinsicSpd5Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              indexMaturityDates[1],
                                              adjustedIssuerCurves) * 10000.0

    intrinsicSpd7Y = cdsIndex.intrinsicSpread(valuationDate,
                                              stepInDate,
                                              indexMaturityDates[2],
                                              adjustedIssuerCurves) * 10000.0

    intrinsicSpd10Y = cdsIndex.intrinsicSpread(valuationDate,
                                               stepInDate,
                                               indexMaturityDates[3],
                                               adjustedIssuerCurves) * 10000.0

    # If the adjustment works then this should equal the index spreads
    testCases.header("LABEL", "VALUE")
    testCases.print("ADJUSTED INTRINSIC SPD 3Y", intrinsicSpd3Y)
    testCases.print("ADJUSTED INTRINSIC SPD 5Y", intrinsicSpd5Y)
    testCases.print("ADJUSTED INTRINSIC SPD 7Y", intrinsicSpd7Y)
    testCases.print("ADJUSTED INTRINSIC SPD 10Y", intrinsicSpd10Y)

###############################################################################


test_performCDSIndexHazardRateAdjustment()
testCases.compareTestCases()
