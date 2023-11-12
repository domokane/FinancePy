###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_option import CDSOption


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################
##########################################################################


def buildFullIssuerCurve(value_date):

    dcType = DayCountTypes.ACT_360
    depos = []
    irBump = 0.0

    m = 1.0  # 0.00000000000

    spot_days = 0
    settle_date = value_date.add_days(spot_days)

    maturity_date = settle_date.add_months(1)
    depo1 = IborDeposit(settle_date, maturity_date, m * 0.0016, dcType)

    maturity_date = settle_date.add_months(2)
    depo2 = IborDeposit(settle_date, maturity_date, m * 0.0020, dcType)

    maturity_date = settle_date.add_months(3)
    depo3 = IborDeposit(settle_date, maturity_date, m * 0.0024, dcType)

    maturity_date = settle_date.add_months(6)
    depo4 = IborDeposit(settle_date, maturity_date, m * 0.0033, dcType)

    maturity_date = settle_date.add_months(12)
    depo5 = IborDeposit(settle_date, maturity_date, m * 0.0056, dcType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    spot_days = 2
    settle_date = value_date.add_days(spot_days)

    swaps = []
    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL

    maturity_date = settle_date.add_months(24)
    swap1 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0044 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settle_date.add_months(36)
    swap2 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0078 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settle_date.add_months(48)
    swap3 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0119 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settle_date.add_months(60)
    swap4 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0158 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settle_date.add_months(72)
    swap5 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0192 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturity_date = settle_date.add_months(84)
    swap6 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0219 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturity_date = settle_date.add_months(96)
    swap7 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0242 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturity_date = settle_date.add_months(108)
    swap8 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0261 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturity_date = settle_date.add_months(120)
    swap9 = IborSwap(
        settle_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.0276 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    libor_curve = IborSingleCurve(value_date, depos, fras, swaps)

    cdsMarketContracts = []
    cds_coupon = 0.005743
    maturity_date = value_date.next_cds_date(6)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.007497
    maturity_date = value_date.next_cds_date(12)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.011132
    maturity_date = value_date.next_cds_date(24)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.013932
    maturity_date = value_date.next_cds_date(36)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.015764
    maturity_date = value_date.next_cds_date(48)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.017366
    maturity_date = value_date.next_cds_date(60)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.020928
    maturity_date = value_date.next_cds_date(84)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    cds_coupon = 0.022835
    maturity_date = value_date.next_cds_date(120)
    cds = CDS(value_date, maturity_date, cds_coupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(value_date,
                            cdsMarketContracts,
                            libor_curve,
                            recovery_rate)

    return libor_curve, issuer_curve

##########################################################################


def test_dirty_priceCDSwaption():

    # This reproduces example on page 38 of Open Gamma note on CDS Option
    tradeDate = Date(5, 2, 2014)
    _, issuer_curve = buildFullIssuerCurve(tradeDate)
    step_in_date = tradeDate.add_days(1)
    value_date = step_in_date
    expiry_date = Date(20, 3, 2014)
    maturity_date = Date(20, 6, 2019)

    cdsRecovery = 0.40
    notional = 100.0
    long_protection = False
    cds_coupon = 0.0  # NOT KNOWN

    cds_contract = CDS(step_in_date,
                       maturity_date,
                       cds_coupon,
                       notional,
                       long_protection)

    testCases.banner(
        "=============================== CDS ===============================")
#    cds_contract.print(value_date)

    testCases.header("LABEL", "VALUE")
    spd = cds_contract.par_spread(
        value_date,
        issuer_curve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD:", spd)

    v = cds_contract.value(value_date, issuer_curve, cdsRecovery)
    testCases.print("DIRTY VALUE", v['dirty_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    p = cds_contract.clean_price(value_date, issuer_curve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accrued_days = cds_contract.accrued_days()
    testCases.print("ACCRUED DAYS", accrued_days)

    accrued_interest = cds_contract.accrued_interest()
    testCases.print("ACCRUED COUPON", accrued_interest)

    prot_pv = cds_contract.protection_leg_pv(
        value_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", prot_pv)

    premPV = cds_contract.premium_leg_pv(
        value_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cds_contract.risky_pv01(
        value_date, issuer_curve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cds_contract.print_payments(issuer_curve)

    testCases.banner(
        "=========================== FORWARD CDS ===========================")

    cds_contract = CDS(expiry_date,
                       maturity_date,
                       cds_coupon,
                       notional,
                       long_protection)

#    cds_contract.print(value_date)

    spd = cds_contract.par_spread(
        value_date,
        issuer_curve,
        cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cds_contract.value(value_date, issuer_curve, cdsRecovery)
    testCases.print("DIRTY VALUE", v['dirty_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    prot_pv = cds_contract.protection_leg_pv(
        value_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", prot_pv)

    premPV = cds_contract.premium_leg_pv(
        value_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    dirtyRPV01, cleanRPV01 = cds_contract.risky_pv01(
        value_date, issuer_curve)
    testCases.print("DIRTY RPV01", dirtyRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cds_contract.print_payments(issuer_curve)

    testCases.banner(
        "========================== CDS OPTIONS ============================")

    cds_coupon = 0.01
    volatility = 0.3
    testCases.print("Expiry Date:", str(expiry_date))
    testCases.print("Maturity Date:", str(maturity_date))
    testCases.print("CDS Coupon:", cds_coupon)

    testCases.header("STRIKE", "LONG PROTECTION", "DIRTY VALUE", "IMPLIED VOL")

    for strike in np.linspace(100, 300, 41):

        long_protection = True  # long protection

        cdsOption = CDSOption(expiry_date,
                              maturity_date,
                              strike / 10000.0,
                              notional,
                              long_protection)

        v = cdsOption.value(value_date,
                            issuer_curve,
                            volatility)

        vol = cdsOption.implied_volatility(value_date,
                                           issuer_curve,
                                           v)

        testCases.print(strike, long_protection, v, vol)

    for strike in np.linspace(100, 300, 41):

        long_protection = False  # long protection

        cdsOption = CDSOption(expiry_date,
                              maturity_date,
                              strike / 10000.0,
                              notional,
                              long_protection)

        v = cdsOption.value(value_date,
                            issuer_curve,
                            volatility)

        vol = cdsOption.implied_volatility(value_date,
                                           issuer_curve,
                                           v)

        testCases.print(strike, long_protection, v, vol)

##########################################################################


test_dirty_priceCDSwaption()
testCases.compareTestCases()
