###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import sys
sys.path.append("..")

from financepy.products.credit.FinCDSOption import FinCDSOption
from financepy.products.credit.cds import FinCDS
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborDeposit import FinIborDeposit
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


def buildFullIssuerCurve(valuation_date):

    dcType = DayCountTypes.ACT_360
    depos = []
    irBump = 0.0

    m = 1.0  # 0.00000000000

    spotDays = 0
    settlement_date = valuation_date.addDays(spotDays)

    maturity_date = settlement_date.addMonths(1)
    depo1 = FinIborDeposit(settlement_date, maturity_date, m * 0.0016, dcType)

    maturity_date = settlement_date.addMonths(2)
    depo2 = FinIborDeposit(settlement_date, maturity_date, m * 0.0020, dcType)

    maturity_date = settlement_date.addMonths(3)
    depo3 = FinIborDeposit(settlement_date, maturity_date, m * 0.0024, dcType)

    maturity_date = settlement_date.addMonths(6)
    depo4 = FinIborDeposit(settlement_date, maturity_date, m * 0.0033, dcType)

    maturity_date = settlement_date.addMonths(12)
    depo5 = FinIborDeposit(settlement_date, maturity_date, m * 0.0056, dcType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    spotDays = 2
    settlement_date = valuation_date.addDays(spotDays)

    swaps = []
    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL

    maturity_date = settlement_date.addMonths(24)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0044 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(36)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0078 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(48)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0119 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(60)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0158 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(72)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0192 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturity_date = settlement_date.addMonths(84)
    swap6 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0219 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturity_date = settlement_date.addMonths(96)
    swap7 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0242 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturity_date = settlement_date.addMonths(108)
    swap8 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0261 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturity_date = settlement_date.addMonths(120)
    swap9 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        m * 0.0276 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    cdsMarketContracts = []
    cdsCoupon = 0.005743
    maturity_date = valuation_date.nextCDSDate(6)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.007497
    maturity_date = valuation_date.nextCDSDate(12)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.011132
    maturity_date = valuation_date.nextCDSDate(24)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.013932
    maturity_date = valuation_date.nextCDSDate(36)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.015764
    maturity_date = valuation_date.nextCDSDate(48)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.017366
    maturity_date = valuation_date.nextCDSDate(60)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.020928
    maturity_date = valuation_date.nextCDSDate(84)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    cdsCoupon = 0.022835
    maturity_date = valuation_date.nextCDSDate(120)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = FinCDSCurve(valuation_date,
                              cdsMarketContracts,
                              libor_curve,
                              recovery_rate)

    return libor_curve, issuer_curve

##########################################################################


def test_full_priceCDSwaption():

    # This reproduces example on page 38 of Open Gamma note on CDS Option
    tradeDate = Date(5, 2, 2014)
    _, issuer_curve = buildFullIssuerCurve(tradeDate)
    step_in_date = tradeDate.addDays(1)
    valuation_date = step_in_date
    expiry_date = Date(20, 3, 2014)
    maturity_date = Date(20, 6, 2019)

    cdsRecovery = 0.40
    notional = 100.0
    long_protection = False
    cdsCoupon = 0.0  # NOT KNOWN

    cdsContract = FinCDS(step_in_date,
                         maturity_date,
                         cdsCoupon,
                         notional,
                         long_protection)

    testCases.banner(
        "=============================== CDS ===============================")
#    cdsContract.print(valuation_date)

    testCases.header("LABEL", "VALUE")
    spd = cdsContract.parSpread(
        valuation_date,
        issuer_curve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD:", spd)

    v = cdsContract.value(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("FULL VALUE", v['full_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    p = cdsContract.clean_price(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accrued_days = cdsContract.accrued_days()
    testCases.print("ACCRUED DAYS", accrued_days)

    accruedInterest = cdsContract.accruedInterest()
    testCases.print("ACCRUED COUPON", accruedInterest)

    protPV = cdsContract.protectionLegPV(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsContract.premiumLegPV(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuation_date, issuer_curve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsContract.printFlows(issuer_curve)

    testCases.banner(
        "=========================== FORWARD CDS ===========================")

    cdsContract = FinCDS(expiry_date,
                         maturity_date,
                         cdsCoupon,
                         notional,
                         long_protection)

#    cdsContract.print(valuation_date)

    spd = cdsContract.parSpread(
        valuation_date,
        issuer_curve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cdsContract.value(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("FULL VALUE", v['full_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    protPV = cdsContract.protectionLegPV(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsContract.premiumLegPV(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsContract.riskyPV01(valuation_date, issuer_curve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsContract.printFlows(issuer_curve)

    testCases.banner(
        "========================== CDS OPTIONS ============================")

    cdsCoupon = 0.01
    volatility = 0.3
    testCases.print("Expiry Date:", str(expiry_date))
    testCases.print("Maturity Date:", str(maturity_date))
    testCases.print("CDS Coupon:", cdsCoupon)

    testCases.header("STRIKE", "FULL VALUE", "IMPLIED VOL")

    for strike in np.linspace(100, 300, 41):

        cdsOption = FinCDSOption(expiry_date,
                                 maturity_date,
                                 strike / 10000.0,
                                 notional)

        v = cdsOption.value(valuation_date,
                            issuer_curve,
                            volatility)

        vol = cdsOption.impliedVolatility(valuation_date,
                                          issuer_curve,
                                          v)

        testCases.print(strike, v, vol)

##########################################################################


test_full_priceCDSwaption()
testCases.compareTestCases()
