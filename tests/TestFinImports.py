###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.finutils import *
from financepy.products.libor import *
from financepy.products.bonds import *
from financepy.market.curves import *


def test_Imports():

    settlementDate = FinDate(1, 1, 2007)
    curve = FinDiscountCurveFlat(settlementDate, 0.05, FinFrequencyTypes.ANNUAL)

    dcType = FinDayCountTypes.ACT_360
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    swapType = FinSwapTypes.PAYER
    swap1 = FinLiborSwap(settlementDate, FinDate(1,1,2008), swapType, 0.05, fixedFreq, dcType)
    swap2 = FinLiborSwap(settlementDate, FinDate(1,1,2009), swapType, 0.05, fixedFreq, dcType)
    swap3 = FinLiborSwap(settlementDate, FinDate(1,1,2010), swapType, 0.05, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3]
    discountCurve = FinLiborCurve(settlementDate, [], [], swaps)

    issueDate = FinDate(1, 1, 2000)
    maturityDate = FinDate(1, 1, 2010)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(issueDate, maturityDate, coupon, frequencyType, accrualType)

###############################################################################

test_Imports()
