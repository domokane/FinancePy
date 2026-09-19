# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.products.bonds.bond_mortgage import BondMortgageTypes
from financepy.products.bonds.bond_mortgage import BondMortgage
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - BondMortgage
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. BOND MORTGAGE
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND MORTGAGE")
print("=" * 78)

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

