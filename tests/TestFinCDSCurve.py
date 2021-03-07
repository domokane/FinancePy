###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.products.credit.FinCDS import FinCDS
from financepy.products.rates.FinIborSwap import FinIborSwap
from financepy.products.credit.FinCDSCurve import FinCDSCurve
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCDSCurve():

    curveDate = FinDate(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = FinDayCountTypes.ACT_365F
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    fixedCoupon = 0.05

    for i in range(1, 11):

        maturityDate = curveDate.addMonths(12 * i)
        swap = FinIborSwap(curveDate,
                           maturityDate,
                           FinSwapTypes.PAY,
                           fixedCoupon,
                           fixedFreq,
                           fixedDCC)
        swaps.append(swap)

    libor_curve = FinIborSingleCurve(curveDate, depos, fras, swaps)

    cdsContracts = []

    for i in range(1, 11):
        maturityDate = curveDate.addMonths(12 * i)
        cds = FinCDS(curveDate, maturityDate, 0.005 + 0.001 * (i - 1))
        cdsContracts.append(cds)

    issuerCurve = FinCDSCurve(curveDate,
                              cdsContracts,
                              libor_curve,
                              recoveryRate=0.40,
                              useCache=False)

    testCases.header("T", "Q")
    n = len(issuerCurve._times)
    for i in range(0, n):
        testCases.print(issuerCurve._times[i], issuerCurve._values[i])

    testCases.header("CONTRACT", "VALUE")
    for i in range(1, 11):
        maturityDate = curveDate.addMonths(12 * i)
        cds = FinCDS(curveDate, maturityDate, 0.005 + 0.001 * (i - 1))
        v = cds.value(curveDate, issuerCurve)
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
