# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path
import numpy as np

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.models.merton_firm_mkt import MertonFirmMkt
from financepy.models.merton_firm import MertonFirm

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_model_merton_credit():

    # -------------------------------------------------------------------------
    # Market Merton model: infer asset value and asset volatility from equity.
    # -------------------------------------------------------------------------

    equity_value = np.array([2.6406, 2.6817, 3.9770, 2.9470, 2.5280])
    equity_vol = np.array([0.7103, 0.3929, 0.3121, 0.4595, 0.6181])
    bond_face = np.array([4.0, 3.5, 3.5, 3.2, 4.0])
    risk_free_rate = 0.05
    asset_growth_rate = np.array([0.0306, 0.0300, 0.0310, 0.0302, 0.0305])
    years_to_maturity = 1.0

    model_mkt = MertonFirmMkt(
        equity_value,
        bond_face,
        years_to_maturity,
        risk_free_rate,
        asset_growth_rate,
        equity_vol,
    )

    test_cases.header("MERTON MARKET MODEL", "VALUE")

    test_cases.print("ASSET VALUE", model_mkt.asset_value())
    test_cases.print("EQUITY VALUE", model_mkt.equity_value())
    test_cases.print("DEBT VALUE", model_mkt.debt_value())

    test_cases.print("ASSET VOLATILITY", model_mkt.asset_volatility())
    test_cases.print("EQUITY VOL", model_mkt.equity_volatility())

    test_cases.print("CREDIT SPREAD", model_mkt.credit_spread())
    test_cases.print("ASSET TO DEBT", model_mkt.asset_to_debt_ratio())
    test_cases.print("RISK NEUTRAL PROB DEFAULT", model_mkt.risk_neutral_default_probability())
    test_cases.print("PHYSICAL PROB DEFAULT", model_mkt.physical_default_probability())
    test_cases.print("DISTANCE DEFAULT", model_mkt.distance_to_default())

    # -------------------------------------------------------------------------
    # Check that the inferred asset quantities reproduce the market inputs.
    # -------------------------------------------------------------------------

    assert np.allclose(
        model_mkt.equity_value(),
        equity_value,
        rtol=1.0e-8,
        atol=1.0e-10,
    )

    assert np.allclose(
        model_mkt.equity_volatility(),
        equity_vol,
        rtol=1.0e-8,
        atol=1.0e-10,
    )

    # -------------------------------------------------------------------------
    # Pass inferred A and sigma_A into the basic Merton model.
    # -------------------------------------------------------------------------

    asset_value = model_mkt.asset_value()
    asset_vol = model_mkt.asset_volatility()

    model = MertonFirm(
        asset_value,
        bond_face,
        years_to_maturity,
        risk_free_rate,
        asset_growth_rate,
        asset_vol,
    )

    test_cases.header("BASIC MERTON MODEL", "VALUE")

    test_cases.print("ASSET VALUE", model.asset_value())
    test_cases.print("EQUITY VALUE", model.equity_value())
    test_cases.print("DEBT VALUE", model.debt_value())

    test_cases.print("ASSET VOLATILITY", model.asset_volatility())
    test_cases.print("EQUITY VOL", model.equity_volatility())

    test_cases.print("CREDIT SPREAD", model.credit_spread())
    test_cases.print("ASSET TO DEBT", model.asset_to_debt_ratio())
    test_cases.print("RISK NEUTRAL DEFAULT PROB", model.risk_neutral_default_probability())
    test_cases.print("PHYSICAL DEFAULT PROB", model.physical_default_probability())
    test_cases.print("DISTANCE DEFAULT", model.distance_to_default())

    # -------------------------------------------------------------------------
    # Basic Merton balance-sheet identity:
    #
    # A = E + D
    # -------------------------------------------------------------------------

    assert np.allclose(
        model.asset_value(),
        model.equity_value() + model.debt_value(),
        rtol=1.0e-12,
        atol=1.0e-12,
    )

    # Market and basic implementations should agree after inversion.

    assert np.allclose(
        model.equity_value(),
        model_mkt.equity_value(),
        rtol=1.0e-10,
        atol=1.0e-12,
    )

    assert np.allclose(
        model.debt_value(),
        model_mkt.debt_value(),
        rtol=1.0e-10,
        atol=1.0e-12,
    )

    # -------------------------------------------------------------------------
    # Scalar example.
    # -------------------------------------------------------------------------

    asset_value = 140.0
    bond_face = 100.0
    years_to_maturity = 1.0
    risk_free_rate = 0.05
    asset_growth_rate = 0.05
    asset_vol = 0.20

    model = MertonFirm(
        asset_value,
        bond_face,
        years_to_maturity,
        risk_free_rate,
        asset_growth_rate,
        asset_vol,
    )

    test_cases.header("BASIC MERTON MODEL SCALAR", "VALUE")

    test_cases.print("ASSET VALUE", model.asset_value())
    test_cases.print("EQUITY VALUE", model.equity_value())
    test_cases.print("DEBT VALUE", model.debt_value())

    test_cases.print("ASSET VOLATILITY", model.asset_volatility())
    test_cases.print("EQUITY VOLATILITY", model.equity_volatility())

    test_cases.print("CREDIT SPREAD", model.credit_spread())
    test_cases.print("ASSET TO DEBT", model.asset_to_debt_ratio())
    test_cases.print("RISK NEUTRAL DEFAULT PROB", model.risk_neutral_default_probability())
    test_cases.print("PHYSICAL DEFAULT PROB", model.physical_default_probability())
    test_cases.print("DISTANCE DEFAULT", model.distance_to_default())

    assert np.allclose(
        model.asset_value(),
        model.equity_value() + model.debt_value(),
        rtol=1.0e-12,
        atol=1.0e-12,
    )


########################################################################################

test_fin_model_merton_credit()
test_cases.compare_test_cases()
