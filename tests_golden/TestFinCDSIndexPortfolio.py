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
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
import os

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

##############################################################################


def buildIssuerCurve(tradeDate, libor_curve):

    valuation_date = tradeDate.add_days(1)

    cdsMarketContracts = []

    cdsCoupon = 0.0048375
    maturity_date = Date(29, 6, 2010)
    cds = CDS(valuation_date, maturity_date, cdsCoupon)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(valuation_date,
                            cdsMarketContracts,
                            libor_curve,
                            recovery_rate)

    return issuer_curve

##########################################################################


def test_CDSIndexPortfolio():

    tradeDate = Date(1, 8, 2007)
    step_in_date = tradeDate.add_days(1)
    valuation_date = step_in_date

    libor_curve = build_Ibor_Curve(tradeDate)

    maturity3Y = tradeDate.next_cds_date(36)
    maturity5Y = tradeDate.next_cds_date(60)
    maturity7Y = tradeDate.next_cds_date(84)
    maturity10Y = tradeDate.next_cds_date(120)

    path = os.path.join(os.path.dirname(__file__),
                        './/data//CDX_NA_IG_S7_SPREADS.csv')
    f = open(path, 'r')
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = CDS(step_in_date, maturity3Y, spd3Y)
        cds5Y = CDS(step_in_date, maturity5Y, spd5Y)
        cds7Y = CDS(step_in_date, maturity7Y, spd7Y)
        cds10Y = CDS(step_in_date, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(valuation_date,
                                cds_contracts,
                                libor_curve,
                                recovery_rate)

        issuer_curves.append(issuer_curve)

    ##########################################################################
    # Now determine the average spread of the index
    ##########################################################################

    cdsIndex = CDSIndexPortfolio()

    averageSpd3Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity3Y,
                                           issuer_curves) * 10000.0

    averageSpd5Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity5Y,
                                           issuer_curves) * 10000.0

    averageSpd7Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity7Y,
                                           issuer_curves) * 10000.0

    averageSpd10Y = cdsIndex.average_spread(valuation_date,
                                            step_in_date,
                                            maturity10Y,
                                            issuer_curves) * 10000.0

    testCases.header("LABEL", "VALUE")
    testCases.print("AVERAGE SPD 3Y", averageSpd3Y)
    testCases.print("AVERAGE SPD 5Y", averageSpd5Y)
    testCases.print("AVERAGE SPD 7Y", averageSpd7Y)
    testCases.print("AVERAGE SPD 10Y", averageSpd10Y)

    ##########################################################################
    # Now determine the intrinsic spread of the index to the same maturity
    # dates. As the single name CDS contracts
    ##########################################################################

    cdsIndex = CDSIndexPortfolio()

    intrinsicSpd3Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity3Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd5Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity5Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd7Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity7Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd10Y = cdsIndex.intrinsic_spread(valuation_date,
                                                step_in_date,
                                                maturity10Y,
                                                issuer_curves) * 10000.0

    ##########################################################################
    ##########################################################################

    testCases.header("LABEL", "VALUE")
    testCases.print("INTRINSIC SPD 3Y", intrinsicSpd3Y)
    testCases.print("INTRINSIC SPD 5Y", intrinsicSpd5Y)
    testCases.print("INTRINSIC SPD 7Y", intrinsicSpd7Y)
    testCases.print("INTRINSIC SPD 10Y", intrinsicSpd10Y)


test_CDSIndexPortfolio()
testCases.compareTestCases()
