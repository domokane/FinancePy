# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.credit.FinCDS import FinCDS
from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinCDSCurve import FinCDSCurve
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinCDSCurve():

    curveDate = FinDate(2018, 12, 20)

    swaps = []
    depos = []
    fras = []

    fixedDCC = FinDayCountTypes.ACT_365_ISDA
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    fixedCoupon = 0.05

    for i in range(1, 11):

        maturityDate = curveDate.addMonths(12 * i)
        swap = FinLiborSwap(
            curveDate,
            maturityDate,
            fixedCoupon,
            fixedFreq,
            fixedDCC)
        swaps.append(swap)

    libor_curve = FinLiborCurve("USD_LIBOR", curveDate, depos, fras, swaps)

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


test_FinCDSCurve()
testCases.compareTestCases()
