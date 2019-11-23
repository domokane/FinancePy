# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""
import numpy as np
import sys
#sys.path.append("..//..")
#
#from financepy import finutils
#from financepy.products import libor
#
#def test_Imports():
#
#    print("Testing Imports")
#    startDate = FinDate(2018, 11, 30)
#    endDate = FinDate(2028, 6, 20)
#
#    endDate = startDate.addMonths(60)
#    oisRate = 0.04
#    isPayer = True
#    fixedFreq = FinFrequencyTypes.ANNUAL
#    fixedDayCount = FinDayCountTypes.ACT_ACT_ISDA
#    floatFreq = FinFrequencyTypes.ANNUAL
#    floatDayCount = FinDayCountTypes.ACT_ACT_ISDA
#    notional = finutils.ONE_MILLION
#
#    ois = FinOIS(startDate,
#                 endDate,
#                 oisRate,
#                 fixedFreq,
#                 fixedDayCount,
#                 floatFreq,
#                 floatDayCount,
#                 isPayer,
#                 notional)
#
#    valueDate = FinDate(2018, 11, 30)
#    marketRate = 0.05
#    indexCurve = FinFlatCurve(
#        valueDate,
#        marketRate,
#        FinCompoundingMethods.ANNUAL)
#    ois.print(valueDate, indexCurve)
#
#    v = ois.value(startDate, indexCurve)
#    print("SWAP_VALUE", v)
#
#
#test_Imports()
