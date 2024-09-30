###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_vars import g_days_in_year
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS


def buildFullIssuerCurve1(mktSpreadBump, irBump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 8-AUG-2019 SNAP AT 1600

    tradeDate = Date(9, 8, 2019)
    value_dt = tradeDate.add_days(1)

    m = 1.0  # 0.00000000000

    dc_type = DayCountTypes.ACT_360
    depos = []
    depo1 = IborDeposit(value_dt, "1D", m * 0.0220, dc_type)
    depos.append(depo1)

    spot_days = 2
    settle_dt = value_dt.add_days(spot_days)

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, m * 0.022009, dc_type)

    maturity_dt = settle_dt.add_months(2)
    depo2 = IborDeposit(settle_dt, maturity_dt, m * 0.022138, dc_type)

    maturity_dt = settle_dt.add_months(3)
    depo3 = IborDeposit(settle_dt, maturity_dt, m * 0.021810, dc_type)

    maturity_dt = settle_dt.add_months(6)
    depo4 = IborDeposit(settle_dt, maturity_dt, m * 0.020503, dc_type)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, m * 0.019930, dc_type)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    swaps = []
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL

    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015910 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014990 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014725 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014640 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014800 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014995 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015180 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015610 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015880 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap9)

    maturity_dt = settle_dt.add_months(144)
    swap10 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.016430 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap10)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    cdsMarketContracts = []

    cds_cpn = 0.04 + mktSpreadBump

    maturity_dt = value_dt.next_cds_date(6)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(12)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(24)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(36)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(48)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(60)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(84)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(120)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    maturity_dt = value_dt.next_cds_date(180)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return libor_curve, issuer_curve


def buildFullIssuerCurve2(mktSpreadBump, irBump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 20 August 2020 SNAP AT 1600

    m = 1.0

    value_dt = Date(24, 8, 2020)
    settle_dt = Date(24, 8, 2020)
    dc_type = DayCountTypes.ACT_360
    depos = []

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, m * 0.001709, dc_type)

    maturity_dt = settle_dt.add_months(2)
    depo2 = IborDeposit(settle_dt, maturity_dt, m * 0.002123, dc_type)

    maturity_dt = settle_dt.add_months(3)
    depo3 = IborDeposit(settle_dt, maturity_dt, m * 0.002469, dc_type)

    maturity_dt = settle_dt.add_months(6)
    depo4 = IborDeposit(settle_dt, maturity_dt, m * 0.003045, dc_type)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, m * 0.004449, dc_type)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    swaps = []
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL

    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.002155 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.002305 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.002665 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.003290 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    libor_curve = IborSingleCurve(value_dt, depos, [], swaps)

    cds_cpn = 0.01 + mktSpreadBump

    cdsMarketContracts = []
    effective_dt = Date(21, 8, 2020)
    cds = CDS(effective_dt, "6M", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "1Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "2Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "3Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "4Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "5Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "7Y", cds_cpn)
    cdsMarketContracts.append(cds)

    cds = CDS(effective_dt, "10Y", cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        settle_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    years = np.linspace(0.0, 10.0, 20)
    dates = settle_dt.add_years(years)
    for dt in dates:
        df = libor_curve.df(dt)
        q = issuer_curve.survival_prob(dt)

    return libor_curve, issuer_curve


cdsRecovery = 0.40

libor_curve, issuer_curve1 = buildFullIssuerCurve1(0.0, 0.0)

# This is the 10 year contract at an off market cpn
maturity_dt = Date(20, 6, 2029)
cds_cpn = 0.0150
notional = ONE_MILLION
long_protection = True
tradeDate = Date(9, 8, 2019)
value_dt1 = tradeDate.add_days(1)
effective_dt = value_dt1

cds_contract1 = CDS(
    effective_dt, maturity_dt, cds_cpn, notional, long_protection
)
t = (maturity_dt - value_dt1) / g_days_in_year
z = libor_curve.df(maturity_dt)
r1 = -np.log(z) / t
print(t, z, r1, maturity_dt)
mktSpread1 = 0.040

libor_curve, issuer_curve2 = buildFullIssuerCurve2(0.0, 0.0)

# This is the 10 year contract at an off market cpn
maturity_dt = Date(20, 6, 2025)
cds_cpn = 0.050
notional = ONE_MILLION
long_protection = True
tradeDate = Date(20, 8, 2020)
effective_dt = Date(21, 8, 2020)
value_dt2 = tradeDate

cds_contract2 = CDS(
    effective_dt, maturity_dt, cds_cpn, notional, long_protection
)
t = (maturity_dt - value_dt2) / g_days_in_year
z = libor_curve.df(maturity_dt)
r2 = -np.log(z) / t
mktSpread2 = 0.01


def test_par_spread():
    spd = (
        cds_contract1.par_spread(value_dt1, issuer_curve1, cdsRecovery)
        * 10000.0
    )
    assert round(spd, 4) == 399.9996

    spd = (
        cds_contract2.par_spread(value_dt2, issuer_curve2, cdsRecovery)
        * 10000.0
    )
    assert round(spd, 4) == 99.5858


def test_value():
    v = cds_contract1.value(value_dt1, issuer_curve1, cdsRecovery)
    assert round(v["dirty_pv"], 4) == 168514.5956
    assert round(v["clean_pv"], 4) == 170639.5956

    v = cds_contract2.value(value_dt2, issuer_curve2, cdsRecovery)
    assert round(v["dirty_pv"], 4) == -199842.7922
    assert round(v["clean_pv"], 4) == -191509.4589


def test_clean_price():
    p = cds_contract1.clean_price(value_dt1, issuer_curve1, cdsRecovery)
    assert round(p, 4) == 82.936

    p = cds_contract2.clean_price(value_dt2, issuer_curve2, cdsRecovery)
    assert round(p, 4) == 119.1509


def testaccrued_days():
    accrued_days = cds_contract1.accrued_days()
    assert accrued_days == 51.0

    accrued_days = cds_contract2.accrued_days()
    assert accrued_days == 60.0


def test_accrued_interest():
    accrued_interest = cds_contract1.accrued_interest()
    assert accrued_interest == -2125.0

    accrued_interest = cds_contract2.accrued_interest()
    assert round(accrued_interest, 4) == -8333.3333


def test_prot_leg_pv():
    prot_pv = cds_contract1.prot_leg_pv(value_dt1, issuer_curve1, cdsRecovery)
    assert round(prot_pv, 4) == 273023.5221

    prot_pv = cds_contract2.prot_leg_pv(value_dt2, issuer_curve2, cdsRecovery)
    assert round(prot_pv, 4) == 47629.7343


def test_premium_leg_pv():
    premPV = cds_contract1.premium_leg_pv(
        value_dt1, issuer_curve1, cdsRecovery
    )
    assert round(premPV, 4) == 104508.9265

    premPV = cds_contract2.premium_leg_pv(
        value_dt2, issuer_curve2, cdsRecovery
    )
    assert round(premPV, 4) == 247472.5265


def test_value_approx():

    v_approx = cds_contract1.value_fast_approx(
        value_dt1, r1, mktSpread1, cdsRecovery
    )
    print(value_dt1, r1, mktSpread1, cdsRecovery)
    assert round(v_approx[0], 4) == 165262.8062
    assert round(v_approx[1], 4) == 167387.8062
    assert round(v_approx[2], 4) == 555.5746
    assert round(v_approx[3], 4) == -71.4881

    v_approx = cds_contract2.value_fast_approx(
        value_dt2, r2, mktSpread2, cdsRecovery
    )
    assert round(v_approx[0], 4) == -195853.3675
    assert round(v_approx[1], 4) == -187520.0342
    assert round(v_approx[2], 4) == 534.9973
    assert round(v_approx[3], 4) == 44.6327
