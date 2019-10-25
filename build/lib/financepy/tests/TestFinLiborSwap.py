# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

import numpy as np

from financepy.finutils.FinDate import FinDate

from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayConventionTypes, FinDateGenRuleTypes
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve

from financepy.products.libor.FinLiborSwap import FinLiborSwap

from financepy.finutils.FinInterpolate import FinInterpMethods
from financepy.finutils.FinMath import ONE_MILLION

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)


def test_LiborSwap(): 
    
    startDate = FinDate(2018,6,20)
    endDate = FinDate(2028,6,20)

    fixedCoupon = 0.05
    fixedFreqType = FinFrequencyTypes.ANNUAL
    fixedDayCountType = FinDayCountTypes.ACT_360

    floatSpread = 0.0
    floatFreqType = FinFrequencyTypes.QUARTERLY
    floatDayCountType = FinDayCountTypes.THIRTY_360

    swapCalendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayConventionTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    fixedCoupon = 0.05
    firstFixing = None
    payFixedFlag = True
    notional = ONE_MILLION
    
    swap = FinLiborSwap(startDate, 
                        endDate,
                        fixedCoupon, 
                        fixedFreqType, 
                        fixedDayCountType,
                        notional,
                        floatSpread, 
                        floatFreqType, 
                        floatDayCountType,
                        firstFixing,
                        payFixedFlag,
                        swapCalendarType,
                        busDayAdjustType,
                        dateGenRuleType)
    
    times = np.linspace(0,10.0,10)
    rate = 0.05
    values = np.power((1+rate),-times)

    discountCurve = FinDiscountCurve(startDate,
                                     times,
                                     values,
                                     FinInterpMethods.FLAT_FORWARDS)

    v = swap.value(startDate,discountCurve)
    testCases.header("LABEL","VALUE")
    testCases.print("SWAP_VALUE",v)

test_LiborSwap()
testCases.compareTestCases()
