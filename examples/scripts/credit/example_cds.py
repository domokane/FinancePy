# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time
import numpy as np


from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.ibor_single_curve import IborSingleCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS

# ============================================================================
# FINANCEPY EXAMPLES - Cds
# ============================================================================

DIRTY = 0
CLEAN = 1




# TO DO

########################################################################################




########################################################################################




########################################################################################




########################################################################################


def test_issuer_curve_build():
    """Test issuer curve build with simple libor curve to isolate cds
    curve building time cost."""

    value_dt = Date(20, 6, 2018)

    times = np.linspace(0.0, 10.0, 11)
    r = 0.05
    discount_factors = np.power((1.0 + r), -times)
    dates = value_dt.add_years(times)
    libor_curve = DiscountCurve(value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES)
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


########################################################################################


def build_full_issuer_curve1(mkt_spd_bump, ir_bump):

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

    cds_cpn = 0.04 + mkt_spd_bump

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


########################################################################################




########################################################################################


def build_full_issuer_curve2(mkt_spd_bump, ir_bump):

    # https://www.markit.com/markit.jsp?jsppage=pv.jsp
    # YIELD CURVE 20 August 2020 SNAP AT 1600

    m = 1.0

    value_dt = Date(20, 8, 2020)
    settle_dt = Date(20, 8, 2020)
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

    cds_cpn = 0.01 + mkt_spd_bump

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

    print("DATE", "DISCOUNT_FACTOR", "SURV_PROB")
    years = np.linspace(0.0, 10.0, 20)
    dates = settle_dt.add_years(years)
    for dt in dates:
        df = libor_curve.df(dt)
        q = issuer_curve.survival_prob(dt)
        print("%16s" % dt, "%12.8f" % df, "%12.8f" % q)

    return libor_curve, issuer_curve


########################################################################################




########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. CDS CURVE BUILD TIMING
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. CDS CURVE BUILD TIMING")
print("=" * 78)

num_curves = 1000

start = time.time()
for _ in range(0, num_curves):
    test_issuer_curve_build()

end = time.time()

print("LABEL", "TIME")
duration = (end - start) / num_curves
print(str(num_curves) + " Libor curves", duration)

# ============================================================================
# 2. DIRTY PRICE CDS MODEL CHECK
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# Calculates coupon interest earned since the previous coupon date and illustrates the clean/dirty price adjustment.

print("\n" + "=" * 78)
print("2. DIRTY PRICE CDS MODEL CHECK")
print("=" * 78)

print("Example", "MARKIT CHECK 19 Aug 2020")

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

print("LABEL", "VALUE")
spd = cds_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
print("PAR_SPREAD", spd)

v = cds_contract.value(value_dt, issuer_curve, cds_recovery)
print("DIRTY_VALUE", v[DIRTY])
print("CLEAN_VALUE", v[CLEAN])

p = cds_contract.clean_price(value_dt, issuer_curve, cds_recovery)
print("CLEAN_PRICE", p)

accrued_days = cds_contract.accrued_days(value_dt)
print("ACCRUED_DAYS", accrued_days)

accrued_interest = cds_contract.accrued_interest(value_dt)
print("ACCRUED_COUPON", accrued_interest)

prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
print("prot_PV", prot_pv)

prem_pv = cds_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
print("PREMIUM_PV", prem_pv)

rpv01 = cds_contract.rpv01(value_dt, issuer_curve)
print("FULL_RPV01", rpv01[DIRTY])
print("CLEAN_RPV01", rpv01[CLEAN])

credit_dv01 = cds_contract.spread_dv01(value_dt, issuer_curve, cds_recovery)
print("CREDIT DV01", credit_dv01)

interest_dv01 = cds_contract.ir_dv01(value_dt, issuer_curve, cds_recovery)
print("INTEREST DV01", interest_dv01)

recovery_dv01 = cds_contract.recovery_dv01(value_dt, issuer_curve, cds_recovery)
print("RECOVERY DV01", recovery_dv01)

#    csa = cds_contract.cash_settlement_amount(value_dt, value_dt, issuer_curve, cds_recovery)
#    print("CSA", csa)

# Consider fast approximation
t = (maturity_dt - value_dt) / G_DAYS_IN_YEAR
z = libor_curve.df(maturity_dt)
r = -np.log(z) / t

mkt_spd = 0.01
v_approx = cds_contract.value_fast_approx(value_dt, r, mkt_spd, cds_recovery)

print("FAST VALUATIONS", "VALUE")

print("DIRTY APPROX VALUE", v_approx[0])
print("CLEAN APPROX VALUE", v_approx[1])
print("APPROX CREDIT DV01", v_approx[2])
print("APPROX INTEREST DV01", v_approx[3])

# ============================================================================
# 3. CDS DATE GENERATION
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("3. CDS DATE GENERATION")
print("=" * 78)

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

print("Flow Date", "AccrualFactor", "Flow")
num_flows = len(cds_contract.payment_dts)
for n in range(0, num_flows):
    print(
        str(cds_contract.payment_dts[n]),
        cds_contract.accrual_factors[n],
        cds_contract.flows[n],
    )

# ============================================================================
# 4. DIRTY PRICE CDS
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# Calculates coupon interest earned since the previous coupon date and illustrates the clean/dirty price adjustment.

print("\n" + "=" * 78)
print("4. DIRTY PRICE CDS")
print("=" * 78)

mkt_spd = 0.040

print("Example", "Markit 9 Aug 2019")

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

print("LABEL", "VALUE")
spd = cds_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
print("PAR_SPREAD", spd)

v = cds_contract.value(value_dt, issuer_curve, cds_recovery)
print("DIRTY_VALUE", v[DIRTY])
print("CLEAN_VALUE", v[CLEAN])

p = cds_contract.clean_price(value_dt, issuer_curve, cds_recovery)
print("CLEAN_PRICE", p)

# MARKIT PRICE IS 168517

accrued_days = cds_contract.accrued_days(value_dt)
print("ACCRUED_DAYS", accrued_days)

accrued_interest = cds_contract.accrued_interest(value_dt)
print("ACCRUED_COUPON", accrued_interest)

prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
print("prot_PV", prot_pv)

prem_pv = cds_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
print("PREMIUM_PV", prem_pv)

dirty_rpv01, clean_rpv01 = cds_contract.rpv01(value_dt, issuer_curve)
print("DIRTY_RPV01", dirty_rpv01)
print("CLEAN_RPV01", clean_rpv01)

# cds_contract.print_payments(issuer_curve)

bump = 1.0 / 10000.0  # 1 bp

libor_curve, issuer_curve = build_full_issuer_curve1(bump, 0)
v_bump = cds_contract.value(value_dt, issuer_curve, cds_recovery)
dv = v_bump[DIRTY] - v[DIRTY]
print("CREDIT_DV01", dv)

# Interest Rate Bump
libor_curve, issuer_curve = build_full_issuer_curve1(0, bump)
v_bump = cds_contract.value(value_dt, issuer_curve, cds_recovery)
dv = v_bump[DIRTY] - v[DIRTY]
print("INTEREST_DV01", dv)

t = (maturity_dt - value_dt) / G_DAYS_IN_YEAR
z = libor_curve.df(maturity_dt)
r = -np.log(z) / t

v_approx = cds_contract.value_fast_approx(value_dt, r, mkt_spd, cds_recovery)

print("DIRTY APPROX VALUE", v_approx[0])
print("CLEAN APPROX VALUE", v_approx[1])
print("DIRTY RPV01 VALUE", v_approx[2])
print("CLEAN RPV01 VALUE", v_approx[3])
print("APPROX SPREAD DV01", v_approx[4])
print("APPROX INTEREST DV01", v_approx[5])
print("APPROX RECOVERY DV01", v_approx[6])

# ============================================================================
# 5. DIRTY PRICE CDS CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("5. DIRTY PRICE CDS CONVERGENCE")
print("=" * 78)

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

print("NumSteps", "Value")
for n in [10, 50, 100, 500, 1000]:
    v_dirty = cds_contract.value(value_dt, issuer_curve, cds_recovery, 0, 1, n)[DIRTY]
    print(n, v_dirty)

# ============================================================================
# 6. CDS CURVE REPRICING
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("6. CDS CURVE REPRICING")
print("=" * 78)

value_dt = Date(20, 6, 2018)
recovery_rate = 0.40

cds_contracts, issuer_curve = test_issuer_curve_build()
print("CDS_MATURITY_dt", "PAR_SPREAD")
for cds in cds_contracts:
    spd = cds.par_spread(value_dt, issuer_curve, recovery_rate)
    print(str(cds.maturity_dt), spd * 10000.0)

# ============================================================================
# 7. CDS FAST APPROXIMATION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("7. CDS FAST APPROXIMATION")
print("=" * 78)

value_dt = Date(20, 6, 2018)
# I build a discount curve that requires no bootstrap
times = np.linspace(0, 10.0, 11)
r = 0.05

discount_factors = np.power((1.0 + r), -times)
dates = value_dt.add_years(times)

libor_curve = DiscountCurve(value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES)

maturity_dt = value_dt.next_cds_date(120)
t = (maturity_dt - value_dt) / 365.242
z = libor_curve.df(maturity_dt)
r = -np.log(z) / t

recovery_rate = 0.40

contract_cpn = 0.010

print("MKT_SPD", "EXACT_VALUE", "APPROX_VALUE", "DIFF(%NOT)")

for mkt_cpn in np.linspace(0.000, 0.05, 21):

    cds_contracts = []

    cds_mkt = CDS(value_dt, maturity_dt, mkt_cpn, ONE_MILLION)

    cds_contracts.append(cds_mkt)

    issuer_curve = CDSCurve(value_dt, cds_contracts, libor_curve, recovery_rate)

    cds_contract = CDS(value_dt, maturity_dt, contract_cpn)
    v_exact = cds_contract.value(value_dt, issuer_curve, recovery_rate)[DIRTY]
    v_approx = cds_contract.value_fast_approx(value_dt, r, mkt_cpn, recovery_rate)[0]
    pct_diff = (v_exact - v_approx) / ONE_MILLION * 100.0
    print(mkt_cpn * 10000, v_exact, v_approx, pct_diff)

