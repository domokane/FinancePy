###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.products.credit.cds import CDS
from financepy.products.rates.IborSwap import FinIborSwap
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCDSCurve():

    curve_date = FinDate(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = FinDayCountTypes.ACT_365F
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    fixedCoupon = 0.05

    for i in range(1, 11):

        maturityDate = curve_date.addMonths(12 * i)
        swap = FinIborSwap(curve_date,
                           maturityDate,
                           FinSwapTypes.PAY,
                           fixedCoupon,
                           fixedFreq,
                           fixedDCC)
        swaps.append(swap)

    libor_curve = FinIborSingleCurve(curve_date, depos, fras, swaps)

    cdsContracts = []

    for i in range(1, 11):
        maturityDate = curve_date.addMonths(12 * i)
        cds = CDS(curve_date, maturityDate, 0.005 + 0.001 * (i - 1))
        cdsContracts.append(cds)

    issuerCurve = CDSCurve(curve_date,
                           cdsContracts,
                           libor_curve,
                           recoveryRate=0.40,
                           use_cache=False)

    testCases.header("T", "Q")
    n = len(issuerCurve._times)
    for i in range(0, n):
        testCases.print(issuerCurve._times[i], issuerCurve._values[i])

    testCases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturityDate = curve_date.addMonths(12 * i)
        cds = CDS(curve_date, maturityDate, 0.005 + 0.001 * (i - 1))
        v = cds.value(curve_date, issuerCurve)
        testCases.print(i, v)

    if 1 == 0:
        x = [0.0, 1.2, 1.6, 1.7, 10.0]
        qs = issuerCurve.survProb(x)
        print("===>", qs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        xx = np.array(x)
        qs = issuerCurve.survProb(xx)
        print("===>", qs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        dfs = issuerCurve.df(x)
        print("===>", dfs)

        x = [0.3, 1.2, 1.6, 1.7, 10.0]
        xx = np.array(x)
        dfs = issuerCurve.df(xx)
        print("===>", dfs)

###############################################################################


test_FinCDSCurve()
testCases.compareTestCases()
