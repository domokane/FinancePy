# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS

test_cases = FinTestCases(__file__, global_test_case_mode)

# TO DO

################################################################################


def build_ibor_curve(trade_dt):

    value_dt = trade_dt.add_days(1)
    dc_type = DayCountTypes.ACT_360
    depos = []

    depos = []
    fras = []
    swaps = []

    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    settle_dt = value_dt

    maturity_dt = settle_dt.add_months(12)
    swap1 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(24)
    swap2 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(36)
    swap3 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type)
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(48)
    swap4 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type)
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(60)
    swap5 = IborSwap(settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type)
    swaps.append(swap5)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve

################################################################################


def build_issuer_curve(trade_dt, libor_curve):

    value_dt = trade_dt.add_days(1)

    cds_mkt_contracts = []

    cds_cpn = 0.0048375
    maturity_dt = Date(20, 6, 2010)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(value_dt, cds_mkt_contracts, libor_curve, recovery_rate)
    return issuer_curve

################################################################################


def test_value_cds_index():

    # We treat an index as a CDS contract with a flat CDS curve
    trade_dt = Date(7, 2, 2006)
    libor_curve = build_ibor_curve(trade_dt)
    issuer_curve = build_issuer_curve(trade_dt, libor_curve)
    step_in_dt = trade_dt.add_days(1)
    value_dt = step_in_dt
    maturity_dt = Date(20, 6, 2010)

    cds_recovery = 0.40
    notional = 10.0 * ONE_MILLION
    long_protection = True
    index_cpn = 0.004

    cds_index_contract = CDS(
        step_in_dt, maturity_dt, index_cpn, notional, long_protection
    )

    #    cds_index_contract.print(value_dt)

    test_cases.header("LABEL", "VALUE")

    spd = cds_index_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
    test_cases.print("PAR SPREAD", spd)

    v = cds_index_contract.value(value_dt, issuer_curve, cds_recovery)
    test_cases.print("DIRTY VALUE", v["dirty_pv"])
    test_cases.print("CLEAN VALUE", v["clean_pv"])

    p = cds_index_contract.clean_price(value_dt, issuer_curve, cds_recovery)
    test_cases.print("CLEAN PRICE", p)

    accrued_days = cds_index_contract.accrued_days()
    test_cases.print("ACCRUED DAYS", accrued_days)

    accrued_interest = cds_index_contract.accrued_interest()
    test_cases.print("ACCRUED COUPON", accrued_interest)

    prot_pv = cds_index_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("PROTECTION LEG PV", prot_pv)

    prem_pv = cds_index_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
    test_cases.print("PREMIUM LEG PV", prem_pv)

    dirty_rpv01, clean_rpv01 = cds_index_contract.risky_pv01(value_dt, issuer_curve)
    test_cases.print("DIRTY RPV01", dirty_rpv01)
    test_cases.print("CLEAN RPV01", clean_rpv01)


################################################################################

#    cds_index_contract.print_payments(issuer_curve)


test_value_cds_index()
test_cases.compare_test_cases()
