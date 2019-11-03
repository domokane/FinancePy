# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinMath import ONE_MILLION

from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.credit.FinCDS import FinCDS

from financepy.market.curves.FinLiborOneCurve import FinLiborOneCurve
from financepy.market.curves.FinCDSCurve import FinCDSCurve


def buildFullIssuerCurve(mktSpreadBump, irBump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 8-AUG-2019 SNAP AT 1600

    tradeDate = FinDate(2019, 8, 9)
    valuationDate = tradeDate.addDays(1)

    dcType = FinDayCountTypes.ACT_360
    depos = []

    m = 1.0  # 0.00000000000

    spotDays = 2
    settlementDate = valuationDate.addDays(spotDays)

    maturityDate = settlementDate.addMonths(1)
    depo1 = FinLiborDeposit(settlementDate, maturityDate, m * 0.022009, dcType)

    maturityDate = settlementDate.addMonths(2)
    depo2 = FinLiborDeposit(settlementDate, maturityDate, m * 0.022138, dcType)

    maturityDate = settlementDate.addMonths(3)
    depo3 = FinLiborDeposit(settlementDate, maturityDate, m * 0.021810, dcType)

    maturityDate = settlementDate.addMonths(6)
    depo4 = FinLiborDeposit(settlementDate, maturityDate, m * 0.020503, dcType)

    maturityDate = settlementDate.addMonths(12)
    depo5 = FinLiborDeposit(settlementDate, maturityDate, m * 0.019930, dcType)

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
        m * 0.015910 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturityDate = settlementDate.addMonths(36)
    swap2 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.014990 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturityDate = settlementDate.addMonths(48)
    swap3 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.014725 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturityDate = settlementDate.addMonths(60)
    swap4 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.014640 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturityDate = settlementDate.addMonths(72)
    swap5 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.014800 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturityDate = settlementDate.addMonths(84)
    swap6 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.014995 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturityDate = settlementDate.addMonths(96)
    swap7 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.015180 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturityDate = settlementDate.addMonths(108)
    swap8 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.015610 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturityDate = settlementDate.addMonths(120)
    swap9 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.015880 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    maturityDate = settlementDate.addMonths(144)
    swap10 = FinLiborSwap(
        settlementDate,
        maturityDate,
        m * 0.016430 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap10)

    liborCurve = FinLiborOneCurve(
        "USD_LIBOR", settlementDate, depos, fras, swaps)

    cdsMarketContracts = []

    cdsCoupon = 0.04 + mktSpreadBump

#    maturityDate = valuationDate.nextCDSDate(6)
#    cds = FinCDS(valuationDate,maturityDate, cdsCoupon)
#    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(12)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(24)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(36)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(48)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(60)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(84)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(120)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturityDate = valuationDate.nextCDSDate(180)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsMarketContracts.append(cds)

#    for cds in cdsMarketContracts:
#        print("CDS Maturity Date",cds._maturityDate)

    recoveryRate = 0.40

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsMarketContracts,
                              liborCurve,
                              recoveryRate)

    return liborCurve, issuerCurve

##########################################################################


def test_fullPriceCDS():

    liborCurve, issuerCurve = buildFullIssuerCurve(0.0, 0.0)

    # This is the 10 year contract at an off market coupon
    maturityDate = FinDate(2029, 6, 20)
    cdsCoupon = 0.0150
    notional = ONE_MILLION
    longProtection = False
    tradeDate = FinDate(2019, 8, 9)
    valuationDate = tradeDate.addDays(1)
    effectiveDate = valuationDate

    cdsContract = FinCDS(effectiveDate,
                         maturityDate,
                         cdsCoupon,
                         notional,
                         longProtection)

    cdsRecovery = 0.40

    print("LABEL", "VALUE")
    spd = cdsContract.parSpread(
        valuationDate,
        issuerCurve,
        cdsRecovery) * 10000.0
    print("PAR_SPREAD", spd)

    v = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    print("FULL_VALUE", v[0])
    print("CLEAN_VALUE", v[1])

    p = cdsContract.cleanPrice(valuationDate, issuerCurve, cdsRecovery)
    print("CLEAN_PRICE", p)

    accruedDays = cdsContract.accruedDays()
    print("ACCRUED_DAYS", accruedDays)

    accruedInterest = cdsContract.accruedInterest()
    print("ACCRUED_COUPON", accruedInterest)

    protPV = cdsContract.protectionLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    print("PROTECTION_PV", protPV)

    premPV = cdsContract.premiumLegPV(valuationDate, issuerCurve, cdsRecovery)
    print("PREMIUM_PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuationDate, issuerCurve)
    print("FULL_RPV01", fullRPV01)
    print("CLEAN_RPV01", cleanRPV01)

    bump = 1.0 / 10000.0  # 1 bp

    liborCurve, issuerCurve = buildFullIssuerCurve(bump, 0)
    v_bump = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    dv = v_bump[0] - v[0]
    print("CREDIT_DV01", dv)

    # Interest Rate Bump
    liborCurve, issuerCurve = buildFullIssuerCurve(0, bump)
    v_bump = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    dv = v_bump[0] - v[0]
    print("INTEREST_DV01", dv)

##########################################################################


test_fullPriceCDS()
