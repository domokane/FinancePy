# ============================================================================
# FINANCEPY EXAMPLES - BondFRN
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the pricing and risk analytics of floating-rate
# notes (FRNs).
#
# An FRN pays a floating reference rate plus a quoted margin. Its market
# price can differ from par when the required discount margin differs from
# the contractual quoted margin.
#
# The examples calculate:
#
#   1. Discount margin
#   2. Dirty price
#   3. Accrued interest
#   4. Principal value
#   5. Interest-rate duration
#   6. Macaulay and modified duration
#   7. Convexity
#   8. Credit-spread duration
#
# Two market examples are considered using different coupon frequencies
# and day-count conventions.
# ============================================================================

from financepy.products.bonds.bond_frn import BondFRN
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# RESULT FORMATTING
# ============================================================================


def print_results(
    dm,
    dirty_price,
    last_coupon_dt,
    accrued_days,
    accrued_amount,
    principal,
    dollar_duration,
    modified_duration,
    macaulay_duration,
    convexity,
    dollar_credit_duration,
    modified_credit_duration,
):
    """Print the principal FRN pricing and risk measures."""

    print("-" * 62)
    print(f"{'MEASURE':<32}{'VALUE':>30}")
    print("-" * 62)

    print(f"{'Discount Margin (bp)':<32}" f"{dm * 10000:>30.6f}")

    print(f"{'Dirty Price':<32}" f"{dirty_price:>30.6f}")

    print(f"{'Last Coupon Date':<32}" f"{str(last_coupon_dt):>30}")

    print(f"{'Accrued Days':<32}" f"{accrued_days:>30}")

    print(f"{'Accrued Amount':<32}" f"{accrued_amount:>30.6f}")

    print(f"{'Principal':<32}" f"{principal:>30.6f}")

    print(f"{'Dollar Rate Duration':<32}" f"{dollar_duration:>30.6f}")

    print(f"{'Modified Rate Duration':<32}" f"{modified_duration:>30.6f}")

    print(f"{'Macaulay Duration':<32}" f"{macaulay_duration:>30.6f}")

    print(f"{'Convexity':<32}" f"{convexity:>30.6f}")

    print(f"{'Dollar Credit Duration':<32}" f"{dollar_credit_duration:>30.6f}")

    print(f"{'Modified Credit Duration':<32}" f"{modified_credit_duration:>30.6f}")

    print("-" * 62)


# ============================================================================
# 1. CITIGROUP FRN - BLOOMBERG EXAMPLE
# ============================================================================
#
# Value a Citigroup floating-rate note using quoted market inputs.
#
# The contractual coupon is the reference IBOR rate plus the quoted margin.
# The discount margin is solved so that the model price reproduces the
# observed clean market price.
#
# Once the discount margin has been determined, the example calculates
# price, accrued interest and a collection of interest-rate and credit-risk
# sensitivity measures.
# ============================================================================

print("\n" + "=" * 78)
print("1. CITIGROUP FRN - BLOOMBERG EXAMPLE")
print("=" * 78)


# ============================================================================
# 1.1 CONTRACT TERMS
# ============================================================================

issue_dt = Date(10, 11, 2010)
maturity_dt = Date(10, 11, 2021)

quoted_margin = 0.0025

freq_type = FrequencyTypes.QUARTERLY
dc_type = DayCountTypes.THIRTY_E_360

bond = BondFRN(
    issue_dt,
    maturity_dt,
    quoted_margin,
    freq_type,
    dc_type,
)


# ============================================================================
# 1.2 MARKET INPUTS
# ============================================================================
#
# reset_ibor is the reference rate associated with the current coupon.
#
# current_ibor and future_ibors provide the reference-rate assumptions used
# to value the remaining floating cash flows.
# ============================================================================

settle_dt = Date(21, 7, 2017)

clean_price = 96.793

reset_ibor = 0.0143456 - quoted_margin
current_ibor = 0.0120534
future_ibors = 0.0130522


# ============================================================================
# 1.3 DISCOUNT MARGIN
# ============================================================================
#
# Solve for the discount margin that reproduces the observed clean price.
# The result is reported in basis points.
# ============================================================================

dm = bond.discount_margin(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    clean_price,
)


# ============================================================================
# 1.4 PRICE AND ACCRUED INTEREST
# ============================================================================

dirty_price = bond.dirty_price_from_dm(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

last_coupon_dt = bond._pcd
accrued_days = bond.accrued_days
accrued_amount = bond.accrued_int


# ============================================================================
# 1.5 PRINCIPAL AND INTEREST-RATE RISK
# ============================================================================
#
# Dollar duration measures the first-order price sensitivity to changes in
# the reference interest rate in price units.
#
# Modified duration expresses the corresponding sensitivity on a relative
# price basis.
#
# Macaulay duration measures the present-value-weighted timing of the
# instrument's cash flows.
#
# Convexity captures the second-order curvature of the price response.
# ============================================================================

principal = bond.principal(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

dollar_duration = bond.dollar_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

modified_duration = bond.modified_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

macaulay_duration = bond.macaulay_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

convexity = bond.convexity_from_dm(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)


# ============================================================================
# 1.6 CREDIT-SPREAD RISK
# ============================================================================
#
# Credit duration measures the sensitivity of the FRN price to a change in
# discount margin while holding the reference-rate assumptions unchanged.
# ============================================================================

dollar_credit_duration = bond.dollar_credit_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

modified_credit_duration = bond.modified_credit_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)


# ============================================================================
# 1.7 RESULTS
# ============================================================================

print_results(
    dm,
    dirty_price,
    last_coupon_dt,
    accrued_days,
    accrued_amount,
    principal,
    dollar_duration,
    modified_duration,
    macaulay_duration,
    convexity,
    dollar_credit_duration,
    modified_credit_duration,
)


# ============================================================================
# 2. CITIGROUP FRN - SECOND MARKET EXAMPLE
# ============================================================================
#
# Repeat the analysis using a second Citigroup FRN example.
#
# This instrument differs from the first example in its issue and maturity
# dates, quoted margin, coupon frequency, day-count convention and market
# assumptions.
#
# Reference:
#
#   https://ebrary.net/14293/economics/actual_floater
# ============================================================================

print("\n" + "=" * 78)
print("2. CITIGROUP FRN - SECOND MARKET EXAMPLE")
print("=" * 78)


# ============================================================================
# 2.1 CONTRACT TERMS
# ============================================================================

issue_dt = Date(28, 3, 2000)
maturity_dt = Date(3, 2, 2021)

quoted_margin = 0.0020

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.THIRTY_E_360_ISDA

bond = BondFRN(
    issue_dt,
    maturity_dt,
    quoted_margin,
    freq_type,
    dc_type,
)


# ============================================================================
# 2.2 MARKET INPUTS
# ============================================================================

settle_dt = Date(28, 3, 2014)

clean_price = 93.08

reset_ibor = 0.00537 - quoted_margin
current_ibor = 0.027558
future_ibors = 0.03295


# ============================================================================
# 2.3 DISCOUNT MARGIN
# ============================================================================

dm = bond.discount_margin(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    clean_price,
)


# ============================================================================
# 2.4 PRICE AND ACCRUED INTEREST
# ============================================================================

dirty_price = bond.dirty_price_from_dm(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

last_coupon_dt = bond._pcd
accrued_days = bond.accrued_days
accrued_amount = bond.accrued_int


# ============================================================================
# 2.5 PRINCIPAL AND INTEREST-RATE RISK
# ============================================================================

principal = bond.principal(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

dollar_duration = bond.dollar_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

modified_duration = bond.modified_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

macaulay_duration = bond.macaulay_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

convexity = bond.convexity_from_dm(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)


# ============================================================================
# 2.6 CREDIT-SPREAD RISK
# ============================================================================

dollar_credit_duration = bond.dollar_credit_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)

modified_credit_duration = bond.modified_credit_duration(
    settle_dt,
    reset_ibor,
    current_ibor,
    future_ibors,
    dm,
)


# ============================================================================
# 2.7 RESULTS
# ============================================================================

print_results(
    dm,
    dirty_price,
    last_coupon_dt,
    accrued_days,
    accrued_amount,
    principal,
    dollar_duration,
    modified_duration,
    macaulay_duration,
    convexity,
    dollar_credit_duration,
    modified_credit_duration,
)
