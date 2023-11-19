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
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes

from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCDSCurve():

    curve_date = Date(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = DayCountTypes.ACT_365F
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    fixed_coupon = 0.05

    for i in range(1, 11):

        maturity_date = curve_date.add_months(12 * i)
        swap = IborSwap(curve_date,
                        maturity_date,
                        SwapTypes.PAY,
                        fixed_coupon,
                        fixedFreq,
                        fixedDCC)
        swaps.append(swap)

    libor_curve = IborSingleCurve(curve_date, depos, fras, swaps)

    cds_contracts = []

    for i in range(1, 11):
        maturity_date = curve_date.add_months(12 * i)
        cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
        cds_contracts.append(cds)

    issuer_curve = CDSCurve(curve_date,
                            cds_contracts,
                            libor_curve,
                            recovery_rate=0.40,
                            use_cache=False)

    testCases.header("T", "Q")
    n = len(issuer_curve._times)
    for i in range(0, n):
        testCases.print(issuer_curve._times[i], issuer_curve._values[i])

    recovery_rate = 0.40

    testCases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturity_date = curve_date.add_months(12 * i)
        cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
        v = cds.value(curve_date, issuer_curve, recovery_rate)
        testCases.print(i, v)

    if 1 == 0:
        x = [0.0, 1.2, 1.6, 1.7, 10.0]
        qs = issuer_curve.survival_prob(x)
        print("===>", qs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        xx = np.array(x)
        qs = issuer_curve.survival_prob(xx)
        print("===>", qs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        dfs = issuer_curve.df(x)
        print("===>", dfs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        xx = np.array(x)
        dfs = issuer_curve.df(xx)
        print("===>", dfs)

###############################################################################


def test_CDS_recovery_rate():

    value_date = Date(15, 8, 2022)
    settle_date = value_date

    swapType = SwapTypes.PAY
    dcType = DayCountTypes.ACT_360
    fixedFreq = FrequencyTypes.MONTHLY
    swap1 = IborSwap(settle_date,"2Y",swapType,0.03512100,fixedFreq,dcType)
    swap2 = IborSwap(settle_date,"3Y",swapType,0.03259300,fixedFreq,dcType)
    swap3 = IborSwap(settle_date,"4Y",swapType,0.03069300,fixedFreq,dcType)
    swap4 = IborSwap(settle_date,"5Y",swapType,0.02952200,fixedFreq,dcType)
    swap5 = IborSwap(settle_date,"6Y",swapType,0.02889300,fixedFreq,dcType)
    swap6 = IborSwap(settle_date,"7Y",swapType,0.02850200,fixedFreq,dcType)
    swap7 = IborSwap(settle_date,"8Y",swapType,0.02827400,fixedFreq,dcType)
    swap8 = IborSwap(settle_date,"9Y",swapType,0.02818500,fixedFreq,dcType)
    swap9 = IborSwap(settle_date,"10Y",swapType,0.02823000,fixedFreq,dcType)
    swaps = [swap1,swap2,swap3,swap4,swap5,swap6,swap7,swap8,swap9]

    libor_curve = IborSingleCurve(value_date, [], [], swaps)

    spreads = [0.000881720, 0.002246440, 0.004283100, 0.005730380, 0.006982450]
    tenors = ['1Y', '3Y', '5Y', '7Y', '10Y']

    cdss = []
    for i in range(len(spreads)):
        freq_type = FrequencyTypes.MONTHLY
        day_count_type = DayCountTypes.ACT_360
        cal_type = CalendarTypes.WEEKEND
        bd_type = BusDayAdjustTypes.FOLLOWING
        dg_type = DateGenRuleTypes.FORWARD
        cds1 = CDS(settle_date,
                   tenors[i],
                   spreads[i],
                   notional=100,
                   freq_type=freq_type,
                   dc_type=day_count_type,
                   cal_type=cal_type,
                   bd_type=bd_type,
                   dg_type=dg_type)

        cdss.append(cds1)

    recovery_rate = 0.575
    issuer_curve = CDSCurve(value_date, cdss, libor_curve, recovery_rate)
#    print(issuer_curve)

    recovery_rate = 0.1
    issuer_curve = CDSCurve(value_date, cdss, libor_curve, recovery_rate)
#    print(issuer_curve)

    recovery_rate = 0.9
    issuer_curve = CDSCurve(value_date, cdss, libor_curve, recovery_rate)
#    print(issuer_curve)


###############################################################################

test_FinCDSCurve()
test_CDS_recovery_rate()

testCases.compareTestCases()
