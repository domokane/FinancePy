# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.credit.FinCDSOption import FinCDSOption
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
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


def buildFullIssuerCurve(tradeDate):

    valuationDate = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360
    depos = []
    irBump = 0.0

    m = 1.0  # 0.00000000000

    spotDays = 2
    settlementDate = valuationDate.addDays(spotDays)

    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(settlementDate, maturityDate, m * 0.0016, dcType)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(settlementDate, maturityDate, m * 0.0020, dcType)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(settlementDate, maturityDate, m * 0.0024, dcType)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(settlementDate, maturityDate, m * 0.0033, dcType)

    maturityDate = settlementDate.addMonths(12)
    depo5 = FinLiborDeposit(settlementDate, maturityDate, m * 0.0056, dcType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    swaps = []
    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL

    maturityDate = settlementDate.addMonths(24)
    swap1 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0044 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0078 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0119 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0158 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0192 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0219 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0242 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0261 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.0276 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    liborCurve = FinLiborCurve(
        "USD_LIBOR", settlementDate, depos, fras, swaps)

    cdsMarketContracts = []

    cdsCoupon = 0.005743
    maturityDate = valuationDate.nextCDSDate(6)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.007497
    maturityDate = valuationDate.nextCDSDate(12)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.011132
    maturityDate = valuationDate.nextCDSDate(24)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.013932
    maturityDate = valuationDate.nextCDSDate(36)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.015764
    maturityDate = valuationDate.nextCDSDate(48)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.017366
    maturityDate = valuationDate.nextCDSDate(60)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.020928
    maturityDate = valuationDate.nextCDSDate(84)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.022835
    maturityDate = valuationDate.nextCDSDate(120)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    recoveryRate = 0.40

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsMarketContracts,
                              liborCurve,
                              recoveryRate)

    return liborCurve, issuerCurve

##########################################################################


def test_fullPriceCDSwaption():

    # This reproduces example on page 38 of Open Gamma note on CDS Option
    tradeDate = FinDate(2014, 2, 5)
    liborCurve, issuerCurve = buildFullIssuerCurve(tradeDate)
    stepInDate = tradeDate.addDays(1)
    valuationDate = stepInDate
    expiryDate = FinDate(2014, 3, 20)
    maturityDate = FinDate(2019, 6, 20)

    cdsRecovery = 0.40
    notional = 100.0
    longProtection = False
    cdsCoupon = 0.0  # NOT KNOWN

    cdsContract = FinCDS(stepInDate,
                         maturityDate,
                         cdsCoupon,
                         notional,
                         longProtection)

    testCases.banner(
        "=============================== CDS ===============================")
#    cdsContract.print(valuationDate)

    testCases.header("LABEL", "VALUE")
    spd = cdsContract.parSpread(
        valuationDate,
        issuerCurve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD:", spd)

    v = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("FULL VALUE", v[0])
    testCases.print("CLEAN VALUE", v[1])

    p = cdsContract.cleanPrice(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accruedDays = cdsContract.accruedDays()
    testCases.print("ACCRUED DAYS", accruedDays)

    accruedInterest = cdsContract.accruedInterest()
    testCases.print("ACCRUED COUPON", accruedInterest)

    protPV = cdsContract.protectionLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsContract.premiumLegPV(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuationDate, issuerCurve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsContract.printFlows(issuerCurve)

    testCases.banner(
        "=========================== FORWARD CDS ===========================")

    cdsContract = FinCDS(expiryDate,
                         maturityDate,
                         cdsCoupon,
                         notional,
                         longProtection)

#    cdsContract.print(valuationDate)

    spd = cdsContract.parSpread(
        valuationDate,
        issuerCurve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("FULL VALUE", v[0])
    testCases.print("CLEAN VALUE", v[1])

    protPV = cdsContract.protectionLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsContract.premiumLegPV(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuationDate, issuerCurve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsContract.printFlows(issuerCurve)

    testCases.banner(
        "========================== CDS OPTIONS ============================")

    cdsCoupon = 0.01
    volatility = 0.3
    testCases.print("Expiry Date:", str(expiryDate))
    testCases.print("Maturity Date:", str(maturityDate))
    testCases.print("CDS Coupon:", cdsCoupon)

    testCases.header("STRIKE", "FULL VALUE", "IMPLIED VOL")

    for strike in np.linspace(100, 300, 21):

        cdsOption = FinCDSOption(expiryDate,
                                 maturityDate,
                                 strike / 10000.0,
                                 notional)

        v = cdsOption.value(valuationDate,
                            issuerCurve,
                            volatility)

        vol = cdsOption.impliedVolatility(valuationDate,
                                          issuerCurve,
                                          v)

        testCases.print(strike, v, vol)

##########################################################################


test_fullPriceCDSwaption()
testCases.compareTestCases()
