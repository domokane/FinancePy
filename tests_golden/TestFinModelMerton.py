###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.merton_firm import MertonFirm
from financepy.models.merton_firm_mkt import MertonFirmMkt
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelMertonCredit():

    # Input Equity values and equity vols
    equity_value = [2.6406, 2.6817, 3.977, 2.947, 2.528]
    equity_vol = [0.7103, 0.3929, 0.3121, 0.4595, 0.6181]
    bondFace = [4.0, 3.5, 3.5, 3.2, 4.0]
    risk_free_rate = [0.05, 0.05, 0.05, 0.05, 0.05]
    assetGrowthRate = [0.0306, 0.03, 0.031, 0.0302, 0.0305]
    timeToMaturity = 1.0  # np.linspace(0.1, 10, 100)

    model = MertonFirmMkt(equity_value,
                          bondFace,
                          timeToMaturity,
                          risk_free_rate,
                          assetGrowthRate,
                          equity_vol)

    testCases.header("MERTON MARKET MODEL", "VALUE")
    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debt_value())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.credit_spread())
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.prob_default())

    assetValue = model._A
    assetVol = model._vA

    model = MertonFirm(assetValue,
                       bondFace,
                       timeToMaturity,
                       risk_free_rate,
                       assetGrowthRate,
                       assetVol)

    testCases.header("BASIC MERTON MODEL", "VALUE")

    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debt_value())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.credit_spread())
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.prob_default())
    testCases.print("DISTANCE DEFAULT", model.dist_default())

    assetValue = 140.0
    bondFace = 100.0
    timeToMaturity = 1.0
    risk_free_rate = 0.05
    assetGrowthRate = 0.05
    assetVol = 0.20

    model = MertonFirm(assetValue,
                       bondFace,
                       timeToMaturity,
                       risk_free_rate,
                       assetGrowthRate,
                       assetVol)

    testCases.header("BASIC MERTON MODEL", "VALUE")

    testCases.print("ASSET VALUE", model._A)
    testCases.print("EQUITY VALUE", model._E)
    testCases.print("DEBT VALUE", model.debt_value())

    testCases.print("ASSET VOLATILITY", model._vA)
    testCases.print("EQUITY VOL", model._vE)

    testCases.print("CREDIT SPREAD", model.credit_spread()*10000)
    testCases.print("LEVERAGE", model.leverage())
    testCases.print("PROD DEFAULT", model.prob_default())
    testCases.print("DISTANCE DEFAULT", model.dist_default())

###############################################################################


test_FinModelMertonCredit()
testCases.compareTestCases()
