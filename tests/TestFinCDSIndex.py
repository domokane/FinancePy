###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.credit.FinCDS import FinCDS
from financepy.utils.Math import ONE_MILLION
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.Date import Date
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
# TO DO
##########################################################################

##########################################################################


def buildIborCurve(tradeDate):

    valuation_date = tradeDate.addDays(1)
    dcType = FinDayCountTypes.ACT_360
    depos = []

    depos = []
    fras = []
    swaps = []

    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    settlement_date = valuation_date

    maturity_date = settlement_date.addMonths(12)
    swap1 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap1)

    maturity_date = settlement_date.addMonths(24)
    swap2 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap2)

    maturity_date = settlement_date.addMonths(36)
    swap3 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap3)

    maturity_date = settlement_date.addMonths(48)
    swap4 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0502,
        fixedFreq,
        dcType)
    swaps.append(swap4)

    maturity_date = settlement_date.addMonths(60)
    swap5 = FinIborSwap(
        settlement_date,
        maturity_date,
        FinSwapTypes.PAY,
        0.0501,
        fixedFreq,
        dcType)
    swaps.append(swap5)

    libor_curve = FinIborSingleCurve(valuation_date, depos, fras, swaps)

    return libor_curve

##########################################################################


def buildIssuerCurve(tradeDate, libor_curve):

    valuation_date = tradeDate.addDays(1)

    cdsMarketContracts = []

    cdsCoupon = 0.0048375
    maturity_date = Date(20, 6, 2010)
    cds = FinCDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = FinCDSCurve(valuation_date,
                              cdsMarketContracts,
                              libor_curve,
                              recovery_rate)
    return issuer_curve

##########################################################################


def test_valueCDSIndex():

    # We treat an index as a CDS contract with a flat CDS curve
    tradeDate = Date(7, 2, 2006)
    libor_curve = buildIborCurve(tradeDate)
    issuer_curve = buildIssuerCurve(tradeDate, libor_curve)
    step_in_date = tradeDate.addDays(1)
    valuation_date = step_in_date
    maturity_date = Date(20, 6, 2010)

    cdsRecovery = 0.40
    notional = 10.0 * ONE_MILLION
    long_protection = True
    indexCoupon = 0.004

    cdsIndexContract = FinCDS(step_in_date,
                              maturity_date,
                              indexCoupon,
                              notional,
                              long_protection)

#    cdsIndexContract.print(valuation_date)

    testCases.header("LABEL", "VALUE")

    spd = cdsIndexContract.parSpread(
        valuation_date, issuer_curve, cdsRecovery) * 10000.0
    testCases.print("PAR SPREAD", spd)

    v = cdsIndexContract.value(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("FULL VALUE", v['full_pv'])
    testCases.print("CLEAN VALUE", v['clean_pv'])

    p = cdsIndexContract.cleanPrice(valuation_date, issuer_curve, cdsRecovery)
    testCases.print("CLEAN PRICE", p)

    accrued_days = cdsIndexContract.accrued_days()
    testCases.print("ACCRUED DAYS", accrued_days)

    accruedInterest = cdsIndexContract.accruedInterest()
    testCases.print("ACCRUED COUPON", accruedInterest)

    protPV = cdsIndexContract.protectionLegPV(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PROTECTION LEG PV", protPV)

    premPV = cdsIndexContract.premiumLegPV(
        valuation_date, issuer_curve, cdsRecovery)
    testCases.print("PREMIUM LEG PV", premPV)

    fullRPV01, cleanRPV01 = cdsIndexContract.riskyPV01(
        valuation_date, issuer_curve)
    testCases.print("FULL  RPV01", fullRPV01)
    testCases.print("CLEAN RPV01", cleanRPV01)

#    cdsIndexContract.printFlows(issuer_curve)


test_valueCDSIndex()
testCases.compareTestCases()
