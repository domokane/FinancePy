###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.credit.FinCDS import FinCDS
from financepy.finutils.FinMath import ONE_MILLION
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

##########################################################################


def buildIborCurve(tradeDate):

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
    maturityDate = FinDate(20, 6, 2010)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    recoveryRate = 0.40

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsMarketContracts,
                              liborCurve,
                              recoveryRate)
    return issuerCurve

##########################################################################


def test_valueCDSIndex():

    # We treat an index as a CDS contract with a flat CDS curve
    tradeDate = FinDate(7, 2, 2006)
    liborCurve = buildIborCurve(tradeDate)
    issuerCurve = buildIssuerCurve(tradeDate, liborCurve)
    stepInDate = tradeDate.addDays(1)
    valuationDate = stepInDate
    maturityDate = FinDate(20, 6, 2010)

    cdsRecovery = 0.40
    notional = 10.0 * ONE_MILLION
    longProtection = True
    indexCoupon = 0.004

    cdsIndexContract = FinCDS(stepInDate,
                              maturityDate,
                              indexCoupon,
                              notional,
                              longProtection)

#    cdsIndexContract.print(valuationDate)

    testCases.header("LABEL", "VALUE")

    spd = cdsIndexContract.parSpread(
        valuationDate, issuerCurve, cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cdsIndexContract.value(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("FULL VALUE", v['full_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    p = cdsIndexContract.cleanPrice(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accruedDays = cdsIndexContract.accruedDays()
    testCases.print("ACCRUED DAYS", accruedDays)

    accruedInterest = cdsIndexContract.accruedInterest()
    testCases.print("ACCRUED COUPON", accruedInterest)

    protPV = cdsIndexContract.protectionLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsIndexContract.premiumLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsIndexContract.riskyPV01(
        valuationDate, issuerCurve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsIndexContract.printFlows(issuerCurve)


test_valueCDSIndex()
testCases.compareTestCases()
