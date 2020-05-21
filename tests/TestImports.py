###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


# from financepy.finutils.FinDate import FinDate # Works
# from ..financepy import *  # fails


from financepy.finutils import *
from financepy.products.libor import *
from financepy.products.bonds import *
from financepy.market.curves import *
# from financepy.products.bonds.FinBond import *


def test_Imports():

    print(dir())

    if 1==1:
        settlementDate = FinDate(1, 1, 2007)
        curve = FinFlatCurve(settlementDate, 0.05, 1)
    
        dcType = FinDayCountTypes.ACT_360
        fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
        swap1 = FinLiborSwap(settlementDate, FinDate(1,1,2008), 0.05, fixedFreq, dcType)
        swap2 = FinLiborSwap(settlementDate, FinDate(1,1,2009), 0.05, fixedFreq, dcType)
        swap3 = FinLiborSwap(settlementDate, FinDate(1,1,2010), 0.05, fixedFreq, dcType)
        swaps = [swap1, swap2, swap3]
        discountCurve = FinLiborCurve("USD_LIBOR", settlementDate, [], [], swaps)
    
        print(discountCurve)

    maturityDate = FinDate(1, 1, 2010)
    coupon = 0.05
    frequencyType = FinFrequencyTypes.ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA

    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)

    print(bond)

###############################################################################

test_Imports()
