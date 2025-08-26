# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.products.bonds.bond_mortgage import BondMortgageTypes
from financepy.products.bonds.bond_mortgage import BondMortgage
from financepy.utils.date import Date

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_bond_mortgage():

    principal = 130000
    start_dt = Date(23, 2, 2018)
    end_dt = start_dt.add_tenor("10Y")
    mortgage = BondMortgage(start_dt, end_dt, principal)

    rate = 0.035
    mortgage.generate_flows(rate, BondMortgageTypes.REPAYMENT)

    num_flows = len(mortgage.schedule.adjusted_dts)

    test_cases.header("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING", "TOTAL")

    for i in range(0, num_flows):
        test_cases.print(
            mortgage.schedule.adjusted_dts[i],
            mortgage.interest_flows[i],
            mortgage.principal_flows[i],
            mortgage.principal_remaining[i],
            mortgage.total_flows[i],
        )

    mortgage.generate_flows(rate, BondMortgageTypes.INTEREST_ONLY)

    test_cases.header("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING", "TOTAL")

    for i in range(0, num_flows):
        test_cases.print(
            mortgage.schedule.adjusted_dts[i],
            mortgage.interest_flows[i],
            mortgage.principal_flows[i],
            mortgage.principal_remaining[i],
            mortgage.total_flows[i],
        )


########################################################################################

test_bond_mortgage()
test_cases.compare_test_cases()
