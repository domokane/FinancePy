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


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################
##########################################################################


def buildFullIssuerCurve(value_dt):

    dc_type = DayCountTypes.ACT_360
    depos = []
    irBump = 0.0

    m = 1.0  # 0.00000000000

    spot_days = 0
    settle_dt = value_dt.add_days(spot_days)

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, m * 0.0016, dc_type)

    maturity_dt = settle_dt.add_months(2)
    depo2 = IborDeposit(settle_dt, maturity_dt, m * 0.0020, dc_type)

    maturity_dt = settle_dt.add_months(3)
    depo3 = IborDeposit(settle_dt, maturity_dt, m * 0.0024, dc_type)

    maturity_dt = settle_dt.add_months(6)
    depo4 = IborDeposit(settle_dt, maturity_dt, m * 0.0033, dc_type)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, m * 0.0056, dc_type)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    spot_days = 2
    settle_dt = value_dt.add_days(spot_days)

    swaps = []
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL

    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0044 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0078 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0119 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0158 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0192 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0219 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0242 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0261 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0276 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap9)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    cdsMarketContracts = []
    cds_cpn = 0.005743
    maturity_dt = value_dt.next_cds_date(6)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.007497
    maturity_dt = value_dt.next_cds_date(12)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.011132
    maturity_dt = value_dt.next_cds_date(24)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.013932
    maturity_dt = value_dt.next_cds_date(36)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.015764
    maturity_dt = value_dt.next_cds_date(48)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.017366
    maturity_dt = value_dt.next_cds_date(60)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.020928
    maturity_dt = value_dt.next_cds_date(84)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.022835
    maturity_dt = value_dt.next_cds_date(120)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return libor_curve, issuer_curve


##########################################################################


def test_dirty_priceCDSwaption():

    # This reproduces example on page 38 of Open Gamma note on CDS Option
    tradeDate = Date(5, 2, 2014)
    _, issuer_curve = buildFullIssuerCurve(tradeDate)
    step_in_dt = tradeDate.add_days(1)
    value_dt = step_in_dt
    expiry_dt = Date(20, 3, 2014)
    maturity_dt = Date(20, 6, 2019)

    cdsRecovery = 0.40
    notional = 100.0
    long_protection = False
    cds_cpn = 0.0  # NOT KNOWN

    cds_contract = CDS(
        step_in_dt, maturity_dt, cds_cpn, notional, long_protection
    )

    test_cases.banner(
        "=============================== CDS ==============================="
    )
    #    cds_contract.print(value_dt)

    test_cases.header("LABEL", "VALUE")
    spd = (
        cds_contract.par_spread(value_dt, issuer_curve, cdsRecovery) * 10000.0
    )
    test_cases.print("PAR SPREAD:", spd)

    v = cds_contract.value(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("DIRTY VALUE", v["dirty_pv"])
    test_cases.print("CLEAN VALUE", v["clean_pv"])

    p = cds_contract.clean_price(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("CLEAN PRICE", p)

    accrued_days = cds_contract.accrued_days()
    test_cases.print("ACCRUED DAYS", accrued_days)

    accrued_interest = cds_contract.accrued_interest()
    test_cases.print("ACCRUED COUPON", accrued_interest)

    prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("PROTECTION LEG PV", prot_pv)

    premPV = cds_contract.premium_leg_pv(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("PREMIUM LEG PV", premPV)

    fullRPV01, clean_rpv01 = cds_contract.risky_pv01(value_dt, issuer_curve)
    test_cases.print("FULL  RPV01", fullRPV01)
    test_cases.print("CLEAN RPV01", clean_rpv01)

    #    cds_contract.print_payments(issuer_curve)

    test_cases.banner(
        "=========================== FORWARD CDS ==========================="
    )

    cds_contract = CDS(
        expiry_dt, maturity_dt, cds_cpn, notional, long_protection
    )

    #    cds_contract.print(value_dt)

    spd = (
        cds_contract.par_spread(value_dt, issuer_curve, cdsRecovery) * 10000.0
    )
    test_cases.print("PAR SPREAD", spd)

    v = cds_contract.value(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("DIRTY VALUE", v["dirty_pv"])
    test_cases.print("CLEAN VALUE", v["clean_pv"])

    prot_pv = cds_contract.prot_leg_pv(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("PROTECTION LEG PV", prot_pv)

    premPV = cds_contract.premium_leg_pv(value_dt, issuer_curve, cdsRecovery)
    test_cases.print("PREMIUM LEG PV", premPV)

    dirty_rpv01, clean_rpv01 = cds_contract.risky_pv01(value_dt, issuer_curve)
    test_cases.print("DIRTY RPV01", dirty_rpv01)
    test_cases.print("CLEAN RPV01", clean_rpv01)

    #    cds_contract.print_payments(issuer_curve)

    test_cases.banner(
        "========================== CDS OPTIONS ============================"
    )

    cds_cpn = 0.01
    volatility = 0.3
    test_cases.print("Expiry Date:", str(expiry_dt))
    test_cases.print("Maturity Date:", str(maturity_dt))
    test_cases.print("CDS Coupon:", cds_cpn)

    test_cases.header(
        "STRIKE", "LONG PROTECTION", "DIRTY VALUE", "IMPLIED VOL"
    )

    for strike in np.linspace(100, 300, 41):

        long_protection = True  # long protection

        cdsOption = CDSOption(
            expiry_dt, maturity_dt, strike / 10000.0, notional, long_protection
        )

        v = cdsOption.value(value_dt, issuer_curve, volatility)

        vol = cdsOption.implied_volatility(value_dt, issuer_curve, v)

        test_cases.print(strike, long_protection, v, vol)

    for strike in np.linspace(100, 300, 41):

        long_protection = False  # long protection

        cdsOption = CDSOption(
            expiry_dt, maturity_dt, strike / 10000.0, notional, long_protection
        )

        v = cdsOption.value(value_dt, issuer_curve, volatility)

        vol = cdsOption.implied_volatility(value_dt, issuer_curve, v)

        test_cases.print(strike, long_protection, v, vol)


##########################################################################


test_dirty_priceCDSwaption()
test_cases.compareTestCases()
