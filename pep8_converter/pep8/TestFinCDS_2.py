# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import numpy as np
import sys

sys.path.append("..")

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS


test_cases = FinTestCases(__file__, global_test_case_mode)

# TO DO


########################################################################################


def test_cds_fast_approximation():

    value_dt = Date(20, 6, 2018)
    # I build a discount curve that requires no bootstrap
    times = np.linspace(0, 10.0, 11)
    r = 0.05

    discount_factors = np.power((1.0 + r), -times)
    dates = value_dt.add_years(times)

    libor_curve = DiscountCurve(
        value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES
    )

    ##########################################################################

    maturity_dt = value_dt.next_cds_date(120)
    t = (maturity_dt - value_dt) / 365.242
    z = libor_curve.df(maturity_dt)
    r = -np.log(z) / t

    recovery_rate = 0.40

    contract_cpn = 0.010

    test_cases.header("MKT_SPD", "EXACT_VALUE", "APPROX_VALUE", "DIFF(%NOT)")

    for mkt_coupon in np.linspace(0.000, 0.05, 21):

        cds_contracts = []

        cds_mkt = CDS(value_dt, maturity_dt, mkt_coupon, ONE_MILLION)

        cds_contracts.append(cds_mkt)

        issuer_curve = CDSCurve(value_dt, cds_contracts, libor_curve, recovery_rate)

        cds_contract = CDS(value_dt, maturity_dt, contract_cpn)
        v_exact = cds_contract.value(value_dt, issuer_curve, recovery_rate)["dirty_pv"]
        v_approx = cds_contract.value_fast_approx(
            value_dt, r, mkt_coupon, recovery_rate
        )[0]
        pct_diff = (v_exact - v_approx) / ONE_MILLION * 100.0
        test_cases.print(mkt_coupon * 10000, v_exact, v_approx, pct_diff)


##########################################################################


########################################################################################


def test_cds_curve_repricing():

    value_dt = Date(20, 6, 2018)
    recovery_rate = 0.40

    cds_contracts, issuer_curve = test_IssuerCurveBuild()
    test_cases.header("CDS_MATURITY_dt", "PAR_SPREAD")
    for cds in cds_contracts:
        spd = cds.par_spread(value_dt, issuer_curve, recovery_rate)
        test_cases.print(str(cds.maturity_dt), spd * 10000.0)


##########################################################################


########################################################################################


def test_cds_curve_build_timing():

    num_curves = 1000

    start = time.time()
    for _ in range(0, num_curves):
        test_IssuerCurveBuild()

    end = time.time()

    test_cases.header("LABEL", "TIME")
    duration = (end - start) / num_curves
    test_cases.print(str(num_curves) + " Libor curves", duration)


##########################################################################


########################################################################################


def test_IssuerCurveBuild():
    """Test issuer curve build with simple libor curve to isolate cds
    curve building time cost."""

    value_dt = Date(20, 6, 2018)

    times = np.linspace(0.0, 10.0, 11)
    r = 0.05
    discount_factors = np.power((1.0 + r), -times)
    dates = value_dt.add_years(times)
    libor_curve = DiscountCurve(
        value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES
    )
    recovery_rate = 0.40

    cds_contracts = []

    cds_cpn = 0.005  # 50 bps
    maturity_dt = value_dt.add_months(12)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    cds_cpn = 0.0055
    maturity_dt = value_dt.add_months(24)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    cds_cpn = 0.0060
    maturity_dt = value_dt.add_months(36)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    cds_cpn = 0.0065
    maturity_dt = value_dt.add_months(60)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    cds_cpn = 0.0070
    maturity_dt = value_dt.add_months(84)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    cds_cpn = 0.0073
    maturity_dt = value_dt.add_months(120)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_contracts.append(cds)

    issuer_curve = CDSCurve(value_dt, cds_contracts, libor_curve, recovery_rate)

    return cds_contracts, issuer_curve


##########################################################################


########################################################################################


def build_full_issuer_curve1(mkt_spread_bump, ir_bump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 8-AUG-2019 SNAP AT 1600

    trade_dt = Date(9, 8, 2019)
    value_dt = trade_dt.add_days(1)

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
        m * 0.015910 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014990 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014725 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014640 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014800 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.014995 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015180 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015610 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.015880 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap9)

    maturity_dt = settle_dt.add_months(144)
    swap10 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.016430 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap10)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    cds_mkt_contracts = []

    cds_cpn = 0.04 + mkt_spread_bump

    maturity_dt = value_dt.next_cds_date(6)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(12)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(24)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(36)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(48)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(60)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(84)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(120)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    maturity_dt = value_dt.next_cds_date(180)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(value_dt, cds_mkt_contracts, libor_curve, recovery_rate)

    return libor_curve, issuer_curve


##########################################################################


########################################################################################


def test_dirty_price_cds():

    mkt_spread = 0.040

    test_cases.header("Example", "Markit 9 Aug 2019")

    libor_curve, issuer_curve = build_full_issuer_curve1(0.0, 0.0)

    # This is the 10 year contract at an off market cpn
    maturity_dt = Date(20, 6, 2029)
    cds_cpn = 0.0150
    notional = ONE_MILLION
    long_protection = True
    trade_dt = Date(9, 8, 2019)
    value_dt = trade_dt.add_days(1)
    effective_dt = value_dt

    cds_contract = CDS(effective_dt, maturity_dt, cds_cpn, notional, long_protection)

    cds_recovery = 0.40

    test_cases.header("LABEL", "VALUE")
    spd = cds_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
    test_cases.print("PAR_SPREAD", spd)

    v = cds_contract.value(value_dt, issuer_curve, cds_recovery)
    test_cases.print("DIRTY_VALUE", v["dirty_pv"])
    test_cases.print("CLEAN_VALUE", v["clean_pv"])

    p = cds_contract.clean_price(value_dt, issuer_curve, cds_recovery)
    test_cases.print("CLEAN_PRICE", p)

    # MARKIT PRICE IS 168517

    accrued_days = cds_contract.accrued_days()
    test_cases.print("ACCRUED_DAYS", accrued_days)

    accrued_interest = cds_contract.accrued_interest()
    test_cases.print("ACCRUED_COUPON", accrued_interest)

    prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("prot_PV", prot_pv)

    prem_pv = cds_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("PREMIUM_PV", prem_pv)

    dirty_rpv01, clean_rpv01 = cds_contract.risky_pv01(value_dt, issuer_curve)
    test_cases.print("DIRTY_RPV01", dirty_rpv01)
    test_cases.print("CLEAN_RPV01", clean_rpv01)

    # cds_contract.print_payments(issuer_curve)

    bump = 1.0 / 10000.0  # 1 bp

    libor_curve, issuer_curve = build_full_issuer_curve1(bump, 0)
    v_bump = cds_contract.value(value_dt, issuer_curve, cds_recovery)
    dv = v_bump["dirty_pv"] - v["dirty_pv"]
    test_cases.print("CREDIT_DV01", dv)

    # Interest Rate Bump
    libor_curve, issuer_curve = build_full_issuer_curve1(0, bump)
    v_bump = cds_contract.value(value_dt, issuer_curve, cds_recovery)
    dv = v_bump["dirty_pv"] - v["dirty_pv"]
    test_cases.print("INTEREST_DV01", dv)

    t = (maturity_dt - value_dt) / G_DAYS_IN_YEARS
    z = libor_curve.df(maturity_dt)
    r = -np.log(z) / t

    v_approx = cds_contract.value_fast_approx(value_dt, r, mkt_spread, cds_recovery)

    test_cases.print("DIRTY APPROX VALUE", v_approx[0])
    test_cases.print("CLEAN APPROX VALUE", v_approx[1])
    test_cases.print("APPROX CREDIT DV01", v_approx[2])
    test_cases.print("APPROX INTEREST DV01", v_approx[3])


##########################################################################


########################################################################################


def build_full_issuer_curve2(mkt_spread_bump, ir_bump):

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
        m * 0.002155 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.002305 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.002665 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.003290 + ir_bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    libor_curve = IborSingleCurve(value_dt, depos, [], swaps)

    cds_cpn = 0.01 + mkt_spread_bump

    cds_mkt_contracts = []
    effective_dt = Date(21, 8, 2020)
    cds = CDS(effective_dt, "6M", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "1Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "2Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "3Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "4Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "5Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "7Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    cds = CDS(effective_dt, "10Y", cds_cpn)
    cds_mkt_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(settle_dt, cds_mkt_contracts, libor_curve, recovery_rate)

    test_cases.header("DATE", "DISCOUNT_FACTOR", "SURV_PROB")
    years = np.linspace(0.0, 10.0, 20)
    dates = settle_dt.add_years(years)
    for dt in dates:
        df = libor_curve.df(dt)
        q = issuer_curve.survival_prob(dt)
        test_cases.print("%16s" % dt, "%12.8f" % df, "%12.8f" % q)

    return libor_curve, issuer_curve


##########################################################################


########################################################################################


def test_dirty_price_cds_model_check():

    test_cases.print("Example", "MARKIT CHECK 19 Aug 2020")

    libor_curve, issuer_curve = build_full_issuer_curve2(0.0, 0.0)

    # This is the 10 year contract at an off market cpn
    maturity_dt = Date(20, 6, 2025)
    cds_cpn = 0.050
    notional = ONE_MILLION
    long_protection = True
    trade_dt = Date(20, 8, 2020)
    effective_dt = Date(21, 8, 2020)
    value_dt = trade_dt

    cds_contract = CDS(effective_dt, maturity_dt, cds_cpn, notional, long_protection)

    cds_recovery = 0.40

    test_cases.header("LABEL", "VALUE")
    spd = cds_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
    test_cases.print("PAR_SPREAD", spd)

    v = cds_contract.value(value_dt, issuer_curve, cds_recovery)
    test_cases.print("FULL_VALUE", v["dirty_pv"])
    test_cases.print("CLEAN_VALUE", v["clean_pv"])

    p = cds_contract.clean_price(value_dt, issuer_curve, cds_recovery)
    test_cases.print("CLEAN_PRICE", p)

    accrued_days = cds_contract.accrued_days()
    test_cases.print("ACCRUED_DAYS", accrued_days)

    accrued_interest = cds_contract.accrued_interest()
    test_cases.print("ACCRUED_COUPON", accrued_interest)

    prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("prot_PV", prot_pv)

    prem_pv = cds_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("PREMIUM_PV", prem_pv)

    rpv01 = cds_contract.risky_pv01(value_dt, issuer_curve)
    test_cases.print("FULL_RPV01", rpv01["dirty_rpv01"])
    test_cases.print("CLEAN_RPV01", rpv01["clean_rpv01"])

    credit_dv01 = cds_contract.credit_dv01(value_dt, issuer_curve, cds_recovery)
    test_cases.print("CREDIT DV01", credit_dv01)

    interest_dv01 = cds_contract.interest_dv01(value_dt, issuer_curve, cds_recovery)
    test_cases.print("INTEREST DV01", interest_dv01)

    # Consider fast approximation
    t = (maturity_dt - value_dt) / G_DAYS_IN_YEARS
    z = libor_curve.df(maturity_dt)
    r = -np.log(z) / t

    mkt_spread = 0.01
    v_approx = cds_contract.value_fast_approx(value_dt, r, mkt_spread, cds_recovery)

    test_cases.header("FAST VALUATIONS", "VALUE")

    test_cases.print("DIRTY APPROX VALUE", v_approx[0])
    test_cases.print("CLEAN APPROX VALUE", v_approx[1])
    test_cases.print("APPROX CREDIT DV01", v_approx[2])
    test_cases.print("APPROX INTEREST DV01", v_approx[3])


##########################################################################


########################################################################################


def test_dirty_price_cds_convergence():

    _, issuer_curve = build_full_issuer_curve1(0.0, 0.0)

    # This is the 10 year contract at an off market cpn
    maturity_dt = Date(20, 6, 2029)
    cds_cpn = 0.0150
    notional = ONE_MILLION
    long_protection = False
    trade_dt = Date(9, 8, 2019)
    value_dt = trade_dt.add_days(1)

    cds_contract = CDS(value_dt, maturity_dt, cds_cpn, notional, long_protection)

    cds_recovery = 0.40

    test_cases.header("num_steps", "Value")
    for n in [10, 50, 100, 500, 1000]:
        v_dirty = cds_contract.value(value_dt, issuer_curve, cds_recovery, 0, 1, n)[
            "dirty_pv"
        ]
        test_cases.print(n, v_dirty)


##########################################################################


########################################################################################


def test_cds_date_generation():

    # This is the 10 year contract at an off market cpn
    maturity_dt = Date(20, 6, 2029)
    cds_cpn = 0.0100

    trade_dt = Date(9, 8, 2019)
    value_dt = trade_dt.add_days(1)

    cds_contract = CDS(
        value_dt,
        maturity_dt,
        cds_cpn,
        ONE_MILLION,
        True,
        FrequencyTypes.QUARTERLY,
        DayCountTypes.ACT_360,
        CalendarTypes.WEEKEND,
        BusDayAdjustTypes.FOLLOWING,
        DateGenRuleTypes.BACKWARD,
    )

    test_cases.header("Flow Date", "accrual_factor", "Flow")
    num_flows = len(cds_contract.payment_dts)
    for n in range(0, num_flows):
        test_cases.print(
            str(cds_contract.payment_dts[n]),
            cds_contract.accrual_factors[n],
            cds_contract.flows[n],
        )


##########################################################################


test_cds_curve_build_timing()
test_dirty_price_cds_model_check()
test_cds_date_generation()
test_dirty_price_cds()
test_dirty_price_cds_convergence()
test_cds_curve_repricing()
test_cds_fast_approximation()

test_cases.compare_test_cases()

########################################################################################
