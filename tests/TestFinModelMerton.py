###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.models.FinModelMertonCredit import FinModelMertonCredit
from financepy.models.FinModelMertonCreditMkt import FinModelMertonCreditMkt

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelMertonCredit():

    # Input Equity values and equity vols
    equityValue = [2.6406, 2.6817, 3.977, 2.947, 2.528]
    equityVol = [0.7103, 0.3929, 0.3121, 0.4595, 0.6181]
    bondFace = [4.0, 3.5, 3.5, 3.2, 4.0]
    riskFreeRate = [0.05, 0.05, 0.05, 0.05, 0.05]
    assetGrowthRate = [0.0306, 0.03, 0.031, 0.0302, 0.0305]
    timeToMaturity = 1.0 # np.linspace(0.1, 10, 100)

    model = FinModelMertonCreditMkt(equityValue,
                                    bondFace,
                                    timeToMaturity,
                                    riskFreeRate,
                                    assetGrowthRate,
                                    equityVol)

    testCases.header("MERTON MARKET MODEL", "VALUE")    
    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debtValue())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.creditSpread())
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.probDefault())

    assetValue = model._A
    assetVol = model._vA
    
    model = FinModelMertonCredit(assetValue,
                                 bondFace,
                                 timeToMaturity,
                                 riskFreeRate,
                                 assetGrowthRate,
                                 assetVol)

    testCases.header("BASIC MERTON MODEL", "VALUE")    

    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debtValue())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.creditSpread())
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.probDefault())
    testCases.print("DISTANCE DEFAULT", model.distDefault())

    assetValue = 140.0
    bondFace = 100.0
    timeToMaturity = 1.0
    riskFreeRate = 0.05
    assetGrowthRate = 0.05
    assetVol = 0.20

    model = FinModelMertonCredit(assetValue,
                                 bondFace,
                                 timeToMaturity,
                                 riskFreeRate,
                                 assetGrowthRate,
                                 assetVol)

    testCases.header("BASIC MERTON MODEL", "VALUE")    

    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debtValue())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.creditSpread()*10000)
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.probDefault())
    testCases.print("DISTANCE DEFAULT", model.distDefault())
    
###############################################################################


test_FinModelMertonCredit()
testCases.compareTestCases()
