###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.merton_firm_mkt import MertonFirmMkt
from financepy.models.merton_firm import MertonFirm
import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelMertonCredit():

    # Input Equity values and equity vols
    equity_value = [2.6406, 2.6817, 3.977, 2.947, 2.528]
    equity_vol = [0.7103, 0.3929, 0.3121, 0.4595, 0.6181]
    bond_face = [4.0, 3.5, 3.5, 3.2, 4.0]
    risk_free_rate = [0.05, 0.05, 0.05, 0.05, 0.05]
    asset_growth_rate = [0.0306, 0.03, 0.031, 0.0302, 0.0305]
    years_to_maturity = 1.0  # np.linspace(0.1, 10, 100)

    model = MertonFirmMkt(equity_value,
                          bond_face,
                          years_to_maturity,
                          risk_free_rate,
                          asset_growth_rate,
                          equity_vol)

    test_cases.header("MERTON MARKET MODEL", "VALUE")
    test_cases.print("ASSET VALUE", model._A)
    test_cases.print("EQUITY VALUE", model._E)
    test_cases.print("DEBT VALUE", model.debt_value())

    test_cases.print("ASSET VOLATILITY", model._vA)
    test_cases.print("EQUITY VOL", model._vE)

    test_cases.print("CREDIT SPREAD", model.credit_spread())
    test_cases.print("LEVERAGE", model.leverage())
    test_cases.print("PROD DEFAULT", model.prob_default())

    asset_value = model._A
    asset_vol = model._vA

    model = MertonFirm(asset_value,
                       bond_face,
                       years_to_maturity,
                       risk_free_rate,
                       asset_growth_rate,
                       asset_vol)

    test_cases.header("BASIC MERTON MODEL", "VALUE")

    test_cases.print("ASSET VALUE", model._A)
    test_cases.print("EQUITY VALUE", model._E)
    test_cases.print("DEBT VALUE", model.debt_value())

    test_cases.print("ASSET VOLATILITY", model._vA)
    test_cases.print("EQUITY VOL", model._vE)

    test_cases.print("CREDIT SPREAD", model.credit_spread())
    test_cases.print("LEVERAGE", model.leverage())
    test_cases.print("PROD DEFAULT", model.prob_default())
    test_cases.print("DISTANCE DEFAULT", model.dist_default())

    asset_value = 140.0
    bond_face = 100.0
    years_to_maturity = 1.0
    risk_free_rate = 0.05
    asset_growth_rate = 0.05
    asset_vol = 0.20

    model = MertonFirm(asset_value,
                       bond_face,
                       years_to_maturity,
                       risk_free_rate,
                       asset_growth_rate,
                       asset_vol)

    test_cases.header("BASIC MERTON MODEL", "VALUE")

    test_cases.print("ASSET VALUE", model._A)
    test_cases.print("EQUITY VALUE", model._E)
    test_cases.print("DEBT VALUE", model.debt_value())

    test_cases.print("ASSET VOLATILITY", model._vA)
    test_cases.print("EQUITY VOL", model._vE)

    test_cases.print("CREDIT SPREAD", model.credit_spread()*10000)
    test_cases.print("LEVERAGE", model.leverage())
    test_cases.print("PROD DEFAULT", model.prob_default())
    test_cases.print("DISTANCE DEFAULT", model.dist_default())

###############################################################################


test_FinModelMertonCredit()
test_cases.compareTestCases()
