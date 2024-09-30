###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import numpy as np

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


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCDSCurve():

    curve_dt = Date(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = DayCountTypes.ACT_365F
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    fixed_cpn = 0.05

    for i in range(1, 11):

        maturity_dt = curve_dt.add_months(12 * i)
        swap = IborSwap(
            curve_dt,
            maturity_dt,
            SwapTypes.PAY,
            fixed_cpn,
            fixed_freq,
            fixedDCC,
        )
        swaps.append(swap)

    libor_curve = IborSingleCurve(curve_dt, depos, fras, swaps)

    cds_contracts = []

    for i in range(1, 11):
        maturity_dt = curve_dt.add_months(12 * i)
        cds = CDS(curve_dt, maturity_dt, 0.005 + 0.001 * (i - 1))
        cds_contracts.append(cds)

    issuer_curve = CDSCurve(
        curve_dt,
        cds_contracts,
        libor_curve,
        recovery_rate=0.40,
        use_cache=False,
    )

    test_cases.header("T", "Q")
    n = len(issuer_curve._times)
    for i in range(0, n):
        test_cases.print(issuer_curve._times[i], issuer_curve._values[i])

    recovery_rate = 0.40

    test_cases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturity_dt = curve_dt.add_months(12 * i)
        cds = CDS(curve_dt, maturity_dt, 0.005 + 0.001 * (i - 1))
        v = cds.value(curve_dt, issuer_curve, recovery_rate)
        test_cases.print(i, v)

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

    value_dt = Date(15, 8, 2022)
    settle_dt = value_dt

    swap_type = SwapTypes.PAY
    dc_type = DayCountTypes.ACT_360
    fixed_freq = FrequencyTypes.MONTHLY
    swap1 = IborSwap(
        settle_dt, "2Y", swap_type, 0.03512100, fixed_freq, dc_type
    )
    swap2 = IborSwap(
        settle_dt, "3Y", swap_type, 0.03259300, fixed_freq, dc_type
    )
    swap3 = IborSwap(
        settle_dt, "4Y", swap_type, 0.03069300, fixed_freq, dc_type
    )
    swap4 = IborSwap(
        settle_dt, "5Y", swap_type, 0.02952200, fixed_freq, dc_type
    )
    swap5 = IborSwap(
        settle_dt, "6Y", swap_type, 0.02889300, fixed_freq, dc_type
    )
    swap6 = IborSwap(
        settle_dt, "7Y", swap_type, 0.02850200, fixed_freq, dc_type
    )
    swap7 = IborSwap(
        settle_dt, "8Y", swap_type, 0.02827400, fixed_freq, dc_type
    )
    swap8 = IborSwap(
        settle_dt, "9Y", swap_type, 0.02818500, fixed_freq, dc_type
    )
    swap9 = IborSwap(
        settle_dt, "10Y", swap_type, 0.02823000, fixed_freq, dc_type
    )
    swaps = [swap1, swap2, swap3, swap4, swap5, swap6, swap7, swap8, swap9]

    libor_curve = IborSingleCurve(value_dt, [], [], swaps)

    spreads = [0.000881720, 0.002246440, 0.004283100, 0.005730380, 0.006982450]
    tenors = ["1Y", "3Y", "5Y", "7Y", "10Y"]

    cdss = []
    for i in range(len(spreads)):
        freq_type = FrequencyTypes.MONTHLY
        day_count_type = DayCountTypes.ACT_360
        cal_type = CalendarTypes.WEEKEND
        bd_type = BusDayAdjustTypes.FOLLOWING
        dg_type = DateGenRuleTypes.FORWARD
        cds1 = CDS(
            settle_dt,
            tenors[i],
            spreads[i],
            notional=100,
            freq_type=freq_type,
            dc_type=day_count_type,
            cal_type=cal_type,
            bd_type=bd_type,
            dg_type=dg_type,
        )

        cdss.append(cds1)

    recovery_rate = 0.575
    issuer_curve = CDSCurve(value_dt, cdss, libor_curve, recovery_rate)
    #    print(issuer_curve)

    recovery_rate = 0.1
    issuer_curve = CDSCurve(value_dt, cdss, libor_curve, recovery_rate)
    #    print(issuer_curve)

    recovery_rate = 0.9
    issuer_curve = CDSCurve(value_dt, cdss, libor_curve, recovery_rate)


#    print(issuer_curve)


###############################################################################

test_FinCDSCurve()
test_CDS_recovery_rate()

test_cases.compareTestCases()
