# ============================================================================
# FINANCEPY EXAMPLES - BondMortgage
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the cash-flow profiles of two common mortgage
# structures:
#
#   1. Repayment mortgage
#   2. Interest-only mortgage
#
# A repayment mortgage pays both interest and principal over its life,
# progressively reducing the outstanding balance.
#
# An interest-only mortgage pays interest during the life of the mortgage,
# with the principal remaining outstanding until maturity.
# ============================================================================

from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

from financepy.products.bonds.bond_mortgage import BondMortgage
from financepy.products.bonds.bond_mortgage import BondMortgageTypes

set_plot_style()

# ============================================================================
# RESULT FORMATTING
# ============================================================================
#
# Display the generated mortgage cash flows in a fixed-width table so that
# the interest, principal and outstanding-balance profiles can be compared
# easily.
# ============================================================================


def print_flows(mortgage):
    """Print the mortgage cash-flow schedule."""

    print(f"{'PAYMENT DATE':<18}" f"{'INTEREST':>15}" f"{'PRINCIPAL':>15}" f"{'OUTSTANDING':>15}" f"{'TOTAL':>15}")

    print("-" * 78)

    num_flows = len(mortgage.schedule.adjusted_dts)

    for i in range(num_flows):

        payment_dt = mortgage.schedule.adjusted_dts[i]
        interest = mortgage.interest_flows[i]
        principal_flow = mortgage.principal_flows[i]
        outstanding = mortgage.principal_remaining[i]
        total = mortgage.total_flows[i]

        print(
            f"{str(payment_dt):<18}"
            f"{interest:>15.2f}"
            f"{principal_flow:>15.2f}"
            f"{outstanding:>15.2f}"
            f"{total:>15.2f}"
        )


# ============================================================================
# 1. MORTGAGE CONTRACT
# ============================================================================
#
# Define a 10-year mortgage with an initial principal of 130,000 and an
# annual interest rate of 3.5%.
#
# The same contract terms are used for both repayment structures so that
# their cash-flow profiles can be compared directly.
# ============================================================================

principal = 130_000
start_dt = Date(23, 2, 2018)
end_dt = start_dt.add_tenor("10Y")
rate = 0.035

mortgage = BondMortgage(
    start_dt,
    end_dt,
    principal,
)


# ============================================================================
# 2. REPAYMENT MORTGAGE
# ============================================================================
#
# Generate the cash flows for a repayment mortgage.
#
# Each regular payment contains:
#
#   - interest on the outstanding principal
#   - repayment of part of the principal
#
# The outstanding balance therefore declines over the life of the mortgage
# until the principal has been fully repaid.
# ============================================================================

print("\n" + "=" * 78)
print("1. REPAYMENT MORTGAGE")
print("=" * 78)

mortgage.generate_flows(
    rate,
    BondMortgageTypes.REPAYMENT,
)

print_flows(mortgage)


# ============================================================================
# 3. INTEREST-ONLY MORTGAGE
# ============================================================================
#
# Generate the cash flows for the same mortgage assuming an interest-only
# structure.
#
# Periodic payments consist of interest while the principal remains
# outstanding. The principal is then repaid at maturity.
#
# Comparing this schedule with the repayment mortgage illustrates how the
# amortisation structure changes both the periodic cash flows and the
# outstanding principal balance.
# ============================================================================

print("\n" + "=" * 78)
print("2. INTEREST-ONLY MORTGAGE")
print("=" * 78)

mortgage.generate_flows(
    rate,
    BondMortgageTypes.INTEREST_ONLY,
)

print_flows(mortgage)
