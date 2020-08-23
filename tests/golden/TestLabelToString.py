# -*- coding: utf-8 -*-

from financepy.finutils.FinHelperFunctions import labelToString

from financepy.finutils import *

from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.products.libor.FinLiborCurve import FinLiborCurve


s = labelToString("FRED", "WILMA")
print(s)

s = labelToString(2.3, 5.3)
print(s)

v = [0.5, 0.2, 0.8]
s = labelToString("Prices", v)


###############################################################################

settlementDate = FinDate(1, 1, 2007)
dcType = FinDayCountTypes.ACT_360
fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
swap1 = FinLiborSwap(settlementDate, FinDate(1,1,2008), 0.05, fixedFreq, dcType)
swap2 = FinLiborSwap(settlementDate, FinDate(1,1,2009), 0.05, fixedFreq, dcType)
swap3 = FinLiborSwap(settlementDate, FinDate(1,1,2010), 0.05, fixedFreq, dcType)
swaps = [swap1, swap2, swap3]
discountCurve = FinLiborCurve(settlementDate, [], [], swaps)

print(discountCurve)
