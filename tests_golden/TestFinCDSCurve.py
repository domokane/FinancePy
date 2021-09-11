###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS
import numpy as np

import sys
sys.path.append("..")


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

    testCases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturity_date = curve_date.add_months(12 * i)
        cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
        v = cds.value(curve_date, issuer_curve)
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


test_FinCDSCurve()
testCases.compareTestCases()
