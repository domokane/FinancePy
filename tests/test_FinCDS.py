###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_vars import gDaysInYear
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS
import time
import numpy as np


def buildFullIssuerCurve1(mktSpreadBump, irBump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 8-AUG-2019 SNAP AT 1600

    tradeDate = Date(9, 8, 2019)
    valuation_date = tradeDate.add_days(1)

    m = 1.0  # 0.00000000000

    dcType = DayCountTypes.ACT_360
    depos = []
    depo1 = IborDeposit(valuation_date, "1D", m * 0.0220, dcType)
    depos.append(depo1)

    spot_days = 2
    settlement_date = valuation_date.add_days(spot_days)

    maturity_date = settlement_date.add_months(1)
    depo1 = IborDeposit(settlement_date, maturity_date, m * 0.022009, dcType)

    maturity_date = settlement_date.add_months(2)
    depo2 = IborDeposit(settlement_date, maturity_date, m * 0.022138, dcType)

    maturity_date = settlement_date.add_months(3)
    depo3 = IborDeposit(settlement_date, maturity_date, m * 0.021810, dcType)

    maturity_date = settlement_date.add_months(6)
    depo4 = IborDeposit(settlement_date, maturity_date, m * 0.020503, dcType)

    maturity_date = settlement_date.add_months(12)
    depo5 = IborDeposit(settlement_date, maturity_date, m * 0.019930, dcType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    swaps = []
    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL

    maturity_date = settlement_date.add_months(24)
    swap1 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.015910 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.add_months(36)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.014990 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.add_months(48)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.014725 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.add_months(60)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.014640 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.add_months(72)
    swap5 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.014800 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    maturity_date = settlement_date.add_months(84)
    swap6 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.014995 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap6)

    maturity_date = settlement_date.add_months(96)
    swap7 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.015180 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap7)

    maturity_date = settlement_date.add_months(108)
    swap8 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.015610 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap8)

    maturity_date = settlement_date.add_months(120)
    swap9 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.015880 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap9)

    maturity_date = settlement_date.add_months(144)
    swap10 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.016430 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap10)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    cdsMarketContracts = []

    cdsCoupon = 0.04 + mktSpreadBump

    maturity_date = valuation_date.next_cds_date(6)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(12)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(24)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(36)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(48)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(60)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(84)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(120)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    maturity_date = valuation_date.next_cds_date(180)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(valuation_date,
                            cdsMarketContracts,
                            libor_curve,
                            recovery_rate)

    return libor_curve, issuer_curve


def buildFullIssuerCurve2(mktSpreadBump, irBump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 20 August 2020 SNAP AT 1600

    m = 1.0

    valuation_date = Date(24, 8, 2020)
    settlement_date = Date(24, 8, 2020)
    dcType = DayCountTypes.ACT_360
    depos = []

    maturity_date = settlement_date.add_months(1)
    depo1 = IborDeposit(settlement_date, maturity_date, m * 0.001709, dcType)

    maturity_date = settlement_date.add_months(2)
    depo2 = IborDeposit(settlement_date, maturity_date, m * 0.002123, dcType)

    maturity_date = settlement_date.add_months(3)
    depo3 = IborDeposit(settlement_date, maturity_date, m * 0.002469, dcType)

    maturity_date = settlement_date.add_months(6)
    depo4 = IborDeposit(settlement_date, maturity_date, m * 0.003045, dcType)

    maturity_date = settlement_date.add_months(12)
    depo5 = IborDeposit(settlement_date, maturity_date, m * 0.004449, dcType)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    swaps = []
    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL

    maturity_date = settlement_date.add_months(24)
    swap1 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.002155 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.add_months(36)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.002305 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.add_months(48)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.002665 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.add_months(60)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        m * 0.003290 + irBump,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    libor_curve = IborSingleCurve(valuation_date, depos, [], swaps)

    cdsCoupon = 0.01 + mktSpreadBump

    cdsMarketContracts = []
    effective_date = Date(21, 8, 2020)
    cds = CDS(effective_date, "6M", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "1Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "2Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "3Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "4Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "5Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "7Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_date, "10Y", cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(settlement_date,
                            cdsMarketContracts,
                            libor_curve,
                            recovery_rate)

    years = np.linspace(0.0, 10.0, 20)
    dates = settlement_date.add_years(years)
    for dt in dates:
        df = libor_curve.df(dt)
        q = issuer_curve.survival_prob(dt)

    return libor_curve, issuer_curve


cdsRecovery = 0.40

libor_curve, issuer_curve1 = buildFullIssuerCurve1(0.0, 0.0)

# This is the 10 year contract at an off market coupon
maturity_date = Date(20, 6, 2029)
cdsCoupon = 0.0150
notional = ONE_MILLION
long_protection = True
tradeDate = Date(9, 8, 2019)
valuation_date1 = tradeDate.add_days(1)
effective_date = valuation_date1

cds_contract1 = CDS(effective_date,
                    maturity_date,
                    cdsCoupon,
                    notional,
                    long_protection)
t = (maturity_date - valuation_date1) / gDaysInYear
z = libor_curve.df(maturity_date)
r1 = -np.log(z) / t
print(t, z, r1, maturity_date)
mktSpread1 = 0.040

libor_curve, issuer_curve2 = buildFullIssuerCurve2(0.0, 0.0)

# This is the 10 year contract at an off market coupon
maturity_date = Date(20, 6, 2025)
cdsCoupon = 0.050
notional = ONE_MILLION
long_protection = True
tradeDate = Date(20, 8, 2020)
effective_date = Date(21, 8, 2020)
valuation_date2 = tradeDate

cds_contract2 = CDS(effective_date,
                    maturity_date,
                    cdsCoupon,
                    notional,
                    long_protection)
t = (maturity_date - valuation_date2) / gDaysInYear
z = libor_curve.df(maturity_date)
r2 = -np.log(z) / t
mktSpread2 = 0.01


def test_par_spread():
    spd = cds_contract1.par_spread(
        valuation_date1,
        issuer_curve1,
        cdsRecovery) * 10000.0
    assert round(spd, 4) == 399.9996

    spd = cds_contract2.par_spread(
        valuation_date2,
        issuer_curve2,
        cdsRecovery) * 10000.0
    assert round(spd, 4) == 98.7148


def test_value():
    v = cds_contract1.value(valuation_date1, issuer_curve1, cdsRecovery)
    assert round(v['dirty_pv'], 4) == 168514.5961
    assert round(v['clean_pv'], 4) == 170639.5961

    v = cds_contract2.value(valuation_date2, issuer_curve2, cdsRecovery)
    assert round(v['dirty_pv'], 4) == -200111.1901
    assert round(v['clean_pv'], 4) == -191777.8568


def test_clean_price():
    p = cds_contract1.clean_price(valuation_date1, issuer_curve1, cdsRecovery)
    assert round(p, 4) == 82.936

    p = cds_contract2.clean_price(valuation_date2, issuer_curve2, cdsRecovery)
    assert round(p, 4) == 119.1926


def test_accrued_days():
    accrued_days = cds_contract1.accrued_days()
    assert accrued_days == 51.0

    accrued_days = cds_contract2.accrued_days()
    assert accrued_days == 60.0


def test_accrued_interest():
    accrued_interest = cds_contract1.accrued_interest()
    assert accrued_interest == -2125.0

    accrued_interest = cds_contract2.accrued_interest()
    assert round(accrued_interest, 4) == -8333.3333


def test_protection_leg_pv():
    prot_pv = cds_contract1.protection_leg_pv(
        valuation_date1, issuer_curve1, cdsRecovery)
    assert round(prot_pv, 4) == 273023.5197

    prot_pv = cds_contract2.protection_leg_pv(
        valuation_date2, issuer_curve2, cdsRecovery)
    assert round(prot_pv, 4) == 47213.1783


def test_premium_leg_pv():
    premPV = cds_contract1.premium_leg_pv(
        valuation_date1, issuer_curve1, cdsRecovery)
    assert round(premPV, 4) == 104508.9236

    premPV = cds_contract2.premium_leg_pv(
        valuation_date2, issuer_curve2, cdsRecovery)
    assert round(premPV, 4) == 247472.5265


def test_value_approx():

    v_approx = cds_contract1.value_fast_approx(valuation_date1,
                                               r1,
                                               mktSpread1,
                                               cdsRecovery)
    print(valuation_date1, r1, mktSpread1, cdsRecovery)
    assert round(v_approx[0], 4) == 165262.8062
    assert round(v_approx[1], 4) == 167387.8062
    assert round(v_approx[2], 4) == 555.5746
    assert round(v_approx[3], 4) == -71.4881

    v_approx = cds_contract2.value_fast_approx(valuation_date2,
                                               r2,
                                               mktSpread2,
                                               cdsRecovery)
    assert round(v_approx[0], 4) == -195853.3675
    assert round(v_approx[1], 4) == -187520.0342
    assert round(v_approx[2], 4) == 534.9973
    assert round(v_approx[3], 4) == 44.6327
