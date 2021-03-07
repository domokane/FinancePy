###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.products.credit.cds import FinCDS
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.credit.cds_curve import FinCDSCurve
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCDSCurve():

    curveDate = Date(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = DayCountTypes.ACT_365F
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    fixedCoupon = 0.05

    for i in range(1, 11):

        maturity_date = curveDate.addMonths(12 * i)
        swap = FinIborSwap(curveDate,
                           maturity_date,
                           FinSwapTypes.PAY,
                           fixedCoupon,
                           fixedFreq,
                           fixedDCC)
        swaps.append(swap)

    libor_curve = IborSingleCurve(curveDate, depos, fras, swaps)

    cds_contracts = []

    for i in range(1, 11):
        maturity_date = curveDate.addMonths(12 * i)
        cds = FinCDS(curveDate, maturity_date, 0.005 + 0.001 * (i - 1))
        cds_contracts.append(cds)

    issuer_curve = FinCDSCurve(curveDate,
                              cds_contracts,
                              libor_curve,
                              recovery_rate=0.40,
                              useCache=False)

    testCases.header("T", "Q")
    n = len(issuer_curve._times)
    for i in range(0, n):
        testCases.print(issuer_curve._times[i], issuer_curve._values[i])

    testCases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturity_date = curveDate.addMonths(12 * i)
        cds = FinCDS(curveDate, maturity_date, 0.005 + 0.001 * (i - 1))
        v = cds.value(curveDate, issuer_curve)
        testCases.print(i, v)

    if 1 == 0:
        x = [0.0, 1.2, 1.6, 1.7, 10.0]
        qs = issuer_curve.survProb(x)
        print("===>", qs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        xx = np.array(x)
        qs = issuer_curve.survProb(xx)
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
