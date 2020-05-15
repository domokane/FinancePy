# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""
# TODO Set up test cases correctly

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.FinBondAnnuity import FinBondAnnuity
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinCalendar import FinDayAdjustTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinLiborCurve import FinLiborCurve

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinBondAnnuity():

    settlementDate = FinDate(20, 6, 2018)

    depos = []
    dcType = FinDayCountTypes.THIRTY_E_360_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    swap1 = FinLiborSwap(settlementDate, "1Y", 0.0502, fixedFreq, dcType)
    swap2 = FinLiborSwap(settlementDate, "2Y", 0.0502, fixedFreq, dcType)
    swap3 = FinLiborSwap(settlementDate, "3Y", 0.0501, fixedFreq, dcType)
    swap4 = FinLiborSwap(settlementDate, "4Y", 0.0502, fixedFreq, dcType)
    swap5 = FinLiborSwap(settlementDate, "5Y", 0.0501, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3, swap4, swap5]

    liborCurve = FinLiborCurve("USD_LIBOR", settlementDate, depos, [], swaps)

    #   print("==============================================================")
    #   print("SEMI-ANNUAL FREQUENCY")
    #   print("==============================================================")

    maturityDate = FinDate(20, 6, 2019)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360
    face = 1000000

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

#    print("===============================================================")
#    print("QUARTERLY FREQUENCY")
#    print("===============================================================")

    maturityDate = FinDate(20, 6, 2028)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

#    print("==================================================================")
#    print("MONTHLY FREQUENCY")
#    print("==================================================================")

    maturityDate = FinDate(20, 6, 2028)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

#    print("==================================================================")
#    print("FORWARD GEN")
#    print("==================================================================")

    maturityDate = FinDate(20, 6, 2028)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

#    print("==================================================================")
#    print("BACKWARD GEN WITH SHORT END STUB")
#    print("==================================================================")

    maturityDate = FinDate(20, 6, 2028)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

#    print("==================================================================")
#    print("FORWARD GEN WITH LONG END STUB")
#    print("==================================================================")

    maturityDate = FinDate(20, 6, 2028)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinBondAnnuity(
        maturityDate,
        coupon,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType,
        face)

    annuity.print()
    annuity.printFlows(settlementDate)
    v = annuity.cleanPriceFromDiscountCurve(settlementDate, liborCurve)
    print("VALUE:", v)

##########################################################################


test_FinBondAnnuity()
testCases.compareTestCases()
