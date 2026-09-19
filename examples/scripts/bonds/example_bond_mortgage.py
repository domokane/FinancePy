# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import add_fp_to_path

from financepy.products.bonds.bond_mortgage import BondMortgageTypes
from financepy.products.bonds.bond_mortgage import BondMortgage
from financepy.utils.date import Date



########################################################################################


def test_bond_mortgage():

    principal = 130000
    start_dt = Date(23, 2, 2018)
    end_dt = start_dt.add_tenor("10Y")
    mortgage = BondMortgage(start_dt, end_dt, principal)

    rate = 0.035
    mortgage.generate_flows(rate, BondMortgageTypes.REPAYMENT)

    num_flows = len(mortgage.schedule.adjusted_dts)

    print("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING", "TOTAL")

    for i in range(0, num_flows):
        print(
            mortgage.schedule.adjusted_dts[i],
            mortgage.interest_flows[i],
            mortgage.principal_flows[i],
            mortgage.principal_remaining[i],
            mortgage.total_flows[i],
        )

    mortgage.generate_flows(rate, BondMortgageTypes.INTEREST_ONLY)

    print("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING", "TOTAL")

    for i in range(0, num_flows):
        print(
            mortgage.schedule.adjusted_dts[i],
            mortgage.interest_flows[i],
            mortgage.principal_flows[i],
            mortgage.principal_remaining[i],
            mortgage.total_flows[i],
        )


########################################################################################

test_bond_mortgage()
