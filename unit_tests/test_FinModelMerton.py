###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.merton_firm_mkt import MertonFirmMkt
from financepy.models.merton_firm import MertonFirm


def test_merton():
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

    assert [round(x, 4) for x in model.debt_value()] == [
        3.7804, 3.3292, 3.3293, 3.0436, 3.7951]
    assert [round(x*1e4, 4) for x in model.credit_spread()
            ] == [64.6893, 0.2289, .0009, 1.2398, 25.7203]
    assert [round(x, 4) for x in model.leverage()] == [
        1.6052, 1.7174, 2.0875, 1.8720, 1.5808]
    assert [round(x*1e2, 4) for x in model.prob_default()
            ] == [6.3791, .0768, 0.0005, 0.2622, 3.4408]

    asset_value = model._A
    asset_vol = model._vA

    model = MertonFirm(asset_value,
                       bond_face,
                       years_to_maturity,
                       risk_free_rate,
                       asset_growth_rate,
                       asset_vol)

    assert [round(x, 4) for x in model.debt_value()] == [
        3.7804, 3.3292, 3.3293, 3.0436, 3.7951]
    assert [round(x*1e4, 4) for x in model.credit_spread()
            ] == [64.6893, 0.2289, 0.0009, 1.2398, 25.7203]
    assert [round(x, 4) for x in model.leverage()] == [
        1.6052, 1.7174, 2.0875, 1.8720, 1.5808]
    assert [round(x*1e2, 4) for x in model.prob_default()
            ] == [6.3791, 0.0768, 0.0005, 0.2622, 3.4408]
    assert [round(x, 4) for x in model.dist_default()] == [
        1.5237, 3.1679, 4.4298, 2.7916, 1.8196]

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

    assert round(model.debt_value(), 4) == 94.8894
    assert round(model.credit_spread()*10000, 4) == 24.5829
    assert model.leverage() == 1.4
    assert round(model.prob_default(), 4) == 0.0334
    assert round(model.dist_default(), 4) == 1.8324
