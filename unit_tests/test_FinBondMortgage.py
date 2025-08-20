########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from financepy.products.bonds.bond_mortgage import BondMortgageTypes
from financepy.products.bonds.bond_mortgage import BondMortgage
from financepy.utils.date import Date


principal = 130000
start_dt = Date(23, 2, 2018)
end_dt = start_dt.add_tenor("10Y")
mortgage = BondMortgage(start_dt, end_dt, principal)

rate = 0.035


def test_bond_mortgage_repayment():
    mortgage.generate_flows(rate, BondMortgageTypes.REPAYMENT)
    i = 10
    assert mortgage.schedule.adjusted_dts[i] == Date(24, 12, 2018)
    assert round(mortgage.interest_flows[i], 4) == 355.0955
    assert round(mortgage.principal_flows[i], 4) == 930.4208
    assert round(mortgage.principal_remaining[i], 4) == 120816.6155
    assert round(mortgage.total_flows[i], 4) == 1285.5163


def test_bond_mortgage_interest():
    mortgage.generate_flows(rate, BondMortgageTypes.INTEREST_ONLY)
    i = 10
    assert mortgage.schedule.adjusted_dts[i] == Date(24, 12, 2018)
    assert round(mortgage.interest_flows[i], 4) == 379.1667
    assert round(mortgage.principal_flows[i], 4) == 0.0000
    assert round(mortgage.principal_remaining[i], 4) == 130000.0000
    assert round(mortgage.total_flows[i], 4) == 379.1667
