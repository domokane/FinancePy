###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.credit.FinCDS import FinCDS
from financepy.finutils.FinMath import ONE_MILLION
from financepy.market.curves.FinInterpolate import FinInterpTypes
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.products.libor.FinLiborCurve import FinLiborCurve
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.finutils.FinGlobalVariables import gDaysInYear
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
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


def test_CDSFastApproximation():

    valueDate = FinDate(2018, 6, 20)
    # I build a discount curve that requires no bootstrap
    times = np.linspace(0, 10.0, 11)
    r = 0.05

    discountFactors = np.power((1.0 + r), -times)
    dates = valueDate.addYears(times)

    liborCurve = FinDiscountCurve(valueDate,
                                  dates,
                                  discountFactors,
                                  FinInterpTypes.FLAT_FORWARDS)

    ##########################################################################

    maturityDate = valueDate.nextCDSDate(120)
    t = (maturityDate - valueDate) / 365.242
    z = liborCurve.df(maturityDate)
    r = -np.log(z) / t

    recoveryRate = 0.40

    contractCoupon = 0.010

    testCases.header("MKT_SPD", "EXACT_VALUE", "APPROX_VALUE", "DIFF(%NOT)")

    for mktCoupon in np.linspace(0.000, 0.05, 21):

        cdsContracts = []

        cdsMkt = FinCDS(valueDate, maturityDate, mktCoupon, ONE_MILLION)

        cdsContracts.append(cdsMkt)

        issuerCurve = FinCDSCurve(valueDate,
                                  cdsContracts,
                                  liborCurve,
                                  recoveryRate)

        cdsContract = FinCDS(valueDate, maturityDate, contractCoupon)
        v_exact = cdsContract.value(
            valueDate, issuerCurve, recoveryRate)['full_pv']
        v_approx = cdsContract.valueFastApprox(
            valueDate, r, mktCoupon, recoveryRate)[0]
        pctdiff = (v_exact - v_approx) / ONE_MILLION * 100.0
        testCases.print(mktCoupon * 10000, v_exact, v_approx, pctdiff)

##########################################################################


def test_CDSCurveRepricing():

    valuationDate = FinDate(2018, 6, 20)
    recoveryRate = 0.40

    cdsContracts, issuerCurve = test_CurveBuild()
    testCases.header("CDS_MATURITY_DATE", "PAR_SPREAD")
    for cds in cdsContracts:
        spd = cds.parSpread(valuationDate, issuerCurve, recoveryRate)
        testCases.print(str(cds._maturityDate), spd * 10000.0)

##########################################################################


def test_CDSCurveBuildTiming():

    numCurves = 10

    start = time.time()
    for i in range(0, numCurves):
        test_CurveBuild()

    end = time.time()

    testCases.header("Label", "TIME")
    duration = (end - start) / numCurves
    testCases.print(str(numCurves) + " Libor curves", duration)

    start = time.time()
    for i in range(0, numCurves):
        buildFullIssuerCurve(0, 0)
    end = time.time()

    duration = (end - start) / numCurves
    testCases.print(str(numCurves) + " Libor + CDS curves:", duration)

##########################################################################


def test_CurveBuild():

    valuationDate = FinDate(2018, 6, 20)

    times = np.linspace(0.0, 10.0, 11)
    r = 0.05
    discountFactors = np.power((1.0 + r), -times)
    dates = valuationDate.addYears(times)
    liborCurve = FinDiscountCurve(valuationDate,
                                  dates,
                                  discountFactors,
                                  FinInterpTypes.FLAT_FORWARDS)
    recoveryRate = 0.40

    cdsContracts = []

    cdsCoupon = 0.005  # 50 bps
    maturityDate = valuationDate.addMonths(12)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    cdsCoupon = 0.0055
    maturityDate = valuationDate.addMonths(24)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    cdsCoupon = 0.0060
    maturityDate = valuationDate.addMonths(36)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    cdsCoupon = 0.0065
    maturityDate = valuationDate.addMonths(60)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    cdsCoupon = 0.0070
    maturityDate = valuationDate.addMonths(84)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    cdsCoupon = 0.0073
    maturityDate = valuationDate.addMonths(120)
    cds = FinCDS(valuationDate, maturityDate, cdsCoupon)
    cdsContracts.append(cds)

    issuerCurve = FinCDSCurve(valuationDate,
                              cdsContracts,
                              liborCurve,
                              recoveryRate)

    return cdsContracts, issuerCurve

##########################################################################


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

    liborCurve = FinLiborCurve(
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

    mktSpread = 0.040

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

    testCases.header("LABEL", "VALUE")
    spd = cdsContract.parSpread(
        valuationDate,
        issuerCurve,
        cdsRecovery) * 10000.0
    testCases.print("PAR_SPREAD", spd)

    v = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("FULL_VALUE", v['full_pv'])
    testCases.print("CLEAN_VALUE", v['clean_pv'])

    p = cdsContract.cleanPrice(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("CLEAN_PRICE", p)

    accruedDays = cdsContract.accruedDays()
    testCases.print("ACCRUED_DAYS", accruedDays)

    accruedInterest = cdsContract.accruedInterest()
    testCases.print("ACCRUED_COUPON", accruedInterest)

    protPV = cdsContract.protectionLegPV(
        valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PROTECTION_PV", protPV)

    premPV = cdsContract.premiumLegPV(valuationDate, issuerCurve, cdsRecovery)
    testCases.print("PREMIUM_PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuationDate, issuerCurve)
    testCases.print("FULL_RPV01", fullRPV01)
    testCases.print("CLEAN_RPV01", cleanRPV01)

    bump = 1.0 / 10000.0  # 1 bp

    liborCurve, issuerCurve = buildFullIssuerCurve(bump, 0)
    v_bump = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    dv = v_bump['full_pv'] - v['full_pv']
    testCases.print("CREDIT_DV01", dv)

    # Interest Rate Bump
    liborCurve, issuerCurve = buildFullIssuerCurve(0, bump)
    v_bump = cdsContract.value(valuationDate, issuerCurve, cdsRecovery)
    dv = v_bump['full_pv'] - v['full_pv']
    testCases.print("INTEREST_DV01", dv)

    t = (maturityDate - valuationDate) / gDaysInYear
    z = liborCurve.df(maturityDate)
    r = -np.log(z) / t

    v_approx = cdsContract.valueFastApprox(valuationDate,
                                           r,
                                           mktSpread,
                                           cdsRecovery)

    testCases.print("FULL APPROX VALUE", v_approx[0])
    testCases.print("CLEAN APPROX VALUE", v_approx[1])
    testCases.print("APPROX CREDIT DV01", v_approx[2])
    testCases.print("APPROX INTEREST DV01", v_approx[3])

##########################################################################


def test_fullPriceCDSConvergence():

    liborCurve, issuerCurve = buildFullIssuerCurve(0.0, 0.0)

    # This is the 10 year contract at an off market coupon
    maturityDate = FinDate(2029, 6, 20)
    cdsCoupon = 0.0150
    notional = ONE_MILLION
    longProtection = False
    tradeDate = FinDate(2019, 8, 9)
    valuationDate = tradeDate.addDays(1)

    cdsContract = FinCDS(valuationDate,
                         maturityDate,
                         cdsCoupon,
                         notional,
                         longProtection)

    cdsRecovery = 0.40

    testCases.header("NumSteps", "Value")
    for n in [10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000]:
        v_full = cdsContract.value(
            valuationDate, issuerCurve, cdsRecovery, 0, 1, n)['full_pv']
        testCases.print(n, v_full)

##########################################################################


def test_CDSDateGeneration():

    # This is the 10 year contract at an off market coupon
    maturityDate = FinDate(2024, 6, 20)
    cdsCoupon = 0.0100

    tradeDate = FinDate(2019, 8, 9)
    valuationDate = tradeDate.addDays(1)

    cdsContract = FinCDS(valuationDate,
                         maturityDate,
                         cdsCoupon,
                         ONE_MILLION,
                         True,
                         FinFrequencyTypes.QUARTERLY,
                         FinDayCountTypes.ACT_360,
                         FinCalendarTypes.WEEKEND,
                         FinBusDayAdjustTypes.FOLLOWING,
                         FinDateGenRuleTypes.BACKWARD)

    testCases.header("Flow Date", "AccrualFactor", "Flow")
    numFlows = len(cdsContract._adjustedDates)
    for n in range(0, numFlows):
        testCases.print(str(
            cdsContract._adjustedDates[n]), cdsContract._accrualFactors[n],
            cdsContract._flows[n])

##########################################################################


test_CDSDateGeneration()
test_fullPriceCDS()
test_fullPriceCDSConvergence()
test_CDSCurveBuildTiming()
test_CDSCurveRepricing()
test_CDSFastApproximation()
testCases.compareTestCases()
