###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################

##########################################################################


def build_Ibor_Curve(tradeDate):

    valuation_date = tradeDate.add_days(1)
    dcType = DayCountTypes.ACT_360
    depos = []

    depos = []
    fras = []
    swaps = []

    dcType = DayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    settlement_date = valuation_date

    maturity_date = settlement_date.add_months(12)
    swap1 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.add_months(24)
    swap2 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.add_months(36)
    swap3 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.add_months(48)
    swap4 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.add_months(60)
    swap5 = IborSwap(
        settlement_date,
        maturity_date,
        SwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    libor_curve = IborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve

##########################################################################


def buildIssuerCurve(tradeDate, libor_curve):

    valuation_date = tradeDate.add_days(1)

    cdsMarketContracts = []

    cdsCoupon = 0.0048375
    maturity_date = Date(20, 6, 2010)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(valuation_date,
                            cdsMarketContracts,
                            libor_curve,
                            recovery_rate)
    return issuer_curve

##########################################################################


def test_valueCDSIndex():

    # We treat an index as a CDS contract with a flat CDS curve
    tradeDate = Date(7, 2, 2006)
    libor_curve = build_Ibor_Curve(tradeDate)
    issuer_curve = buildIssuerCurve(tradeDate, libor_curve)
    step_in_date = tradeDate.add_days(1)
    valuation_date = step_in_date
    maturity_date = Date(20, 6, 2010)

    cdsRecovery = 0.40
    notional = 10.0 * ONE_MILLION
    long_protection = True
    index_coupon = 0.004

    cdsIndexContract = CDS(step_in_date,
                           maturity_date,
                           index_coupon,
                           notional,
                           long_protection)

#    cdsIndexContract.print(valuation_date)

    testCases.header("LABEL", "VALUE")

    spd = cdsIndexContract.par_spread(
        valuation_date, issuer_curve, cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cdsIndexContract.value(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("FULL VALUE", v['full_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    p = cdsIndexContract.clean_price(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accrued_days = cdsIndexContract.accrued_days()
    testCases.print("ACCRUED DAYS", accrued_days)

    accrued_interest = cdsIndexContract.accrued_interest()
    testCases.print("ACCRUED COUPON", accrued_interest)

    prot_pv = cdsIndexContract.protection_leg_pv(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", prot_pv)

    premPV = cdsIndexContract.premium_leg_pv(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsIndexContract.risky_pv01(
        valuation_date, issuer_curve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsIndexContract.print_flows(issuer_curve)


test_valueCDSIndex()
testCases.compareTestCases()
