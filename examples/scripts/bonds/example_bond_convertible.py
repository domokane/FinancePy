# ============================================================================
# FINANCEPY EXAMPLES - BondConvertible
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of a callable and puttable
# convertible bond.
#
# A convertible bond combines a conventional bond with an embedded option
# allowing the holder to convert the bond into equity. Its value therefore
# depends on both fixed-income and equity-market inputs.
#
# The example examines:
#
#   1. Convertible-bond valuation with zero dividend yield
#   2. Numerical convergence as the number of tree steps is increased
#   3. The effect of introducing a positive dividend yield
#   4. The computational time required at different tree resolutions
#
# The convertible also contains issuer call and investor put schedules,
# which affect the exercise decisions represented in the valuation model.
# ============================================================================

import time

import numpy as np

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

from financepy.products.bonds.bond_convertible import BondConvertible

set_plot_style()

# ============================================================================
# 1. CONVERTIBLE BOND CONTRACT
# ============================================================================
#
# Define the contractual terms of the convertible bond, including its
# coupon, maturity, conversion ratio and embedded call and put schedules.
#
# The conversion ratio specifies the number of shares received when the
# bond is converted into equity.
# ============================================================================

settle_dt = Date(31, 12, 2003)
start_convert_date = Date(31, 12, 2003)
maturity_dt = Date(15, 3, 2022)

conversion_ratio = 38.4615
coupon = 0.0575

freq_type = FrequencyTypes.SEMI_ANNUAL
accrual_basis = DayCountTypes.ACT_365F


# ============================================================================
# 2. CALL SCHEDULE
# ============================================================================
#
# The call schedule gives the issuer the right to redeem the convertible
# at specified dates and prices. This can limit the value of the
# convertible to the investor when the underlying equity performs strongly.
# ============================================================================

call_price = 1100

call_dts = [
    Date(20, 3, 2007),
    Date(15, 3, 2012),
    Date(15, 3, 2017),
]

call_prices = np.array(
    [
        call_price,
        call_price,
        call_price,
    ]
)


# ============================================================================
# 3. PUT SCHEDULE
# ============================================================================
#
# The put schedule gives the investor the right to sell the convertible
# back to the issuer at specified dates and prices.
# ============================================================================

put_price = 90

put_dts = [
    Date(20, 3, 2007),
    Date(15, 3, 2012),
    Date(15, 3, 2017),
]

put_prices = np.array(
    [
        put_price,
        put_price,
        put_price,
    ]
)


bond = BondConvertible(
    maturity_dt,
    coupon,
    freq_type,
    start_convert_date,
    conversion_ratio,
    call_dts,
    call_prices,
    put_dts,
    put_prices,
    accrual_basis,
)


# ============================================================================
# 4. EQUITY DIVIDEND SCHEDULE
# ============================================================================
#
# Convertible-bond valuation depends on dividends paid by the underlying
# equity because dividends affect the economics of converting the bond
# into shares.
# ============================================================================

dividend_dts = [
    Date(20, 3, 2007),
    Date(15, 3, 2008),
    Date(15, 3, 2009),
    Date(15, 3, 2010),
    Date(15, 3, 2011),
    Date(15, 3, 2012),
    Date(15, 3, 2013),
    Date(15, 3, 2014),
    Date(15, 3, 2015),
    Date(15, 3, 2016),
    Date(15, 3, 2017),
    Date(15, 3, 2018),
    Date(15, 3, 2019),
    Date(15, 3, 2020),
    Date(15, 3, 2021),
    Date(15, 3, 2022),
]


# ============================================================================
# 5. MARKET AND MODEL INPUTS
# ============================================================================
#
# The equity component is driven by the stock price and volatility.
#
# The fixed-income component uses a flat continuously compounded discount
# curve together with a credit spread and recovery-rate assumption.
# ============================================================================

stock_price = 28.5
stock_volatility = 0.370

rate = 0.04

discount_curve = FlatDiscountCurve(
    settle_dt,
    rate,
    FrequencyTypes.CONTINUOUS,
)

credit_spread = 0.00
recovery_rate = 0.40


# ============================================================================
# 6. RESULT FORMATTING
# ============================================================================
#
# BondConvertible.value() returns several model outputs. Displaying the
# individual values in columns makes it easier to compare convergence as
# the number of tree steps is increased.
#
# CB PRICE : full convertible-bond value
# BOND     : straight-bond component
# DELTA    : sensitivity to the underlying stock price
# GAMMA    : change in delta with respect to the stock price
# THETA    : sensitivity to the passage of time
# TIME     : elapsed valuation time in seconds
# ============================================================================


def print_header():
    """Print the valuation-results table header."""

    print(
        f"{'STEPS':>8}"
        f"{'CB PRICE':>14}"
        f"{'BOND':>14}"
        f"{'DELTA':>14}"
        f"{'GAMMA':>14}"
        f"{'THETA':>14}"
        f"{'TIME':>12}"
    )

    print("-" * 90)


def print_result(num_steps_per_year, result, elapsed):
    """Print one convertible-bond valuation result."""

    print(
        f"{num_steps_per_year:8d}"
        f"{result['cbprice']:14.6f}"
        f"{result['bond']:14.6f}"
        f"{result['delta']:14.6f}"
        f"{result['gamma']:14.6f}"
        f"{result['theta']:14.6f}"
        f"{elapsed:12.6f}"
    )


def run_valuations(dividend_yields, steps):
    """Value the convertible for each requested tree resolution."""

    print_header()

    for num_steps_per_year in steps:

        start = time.perf_counter()

        result = bond.value(
            settle_dt,
            stock_price,
            stock_volatility,
            dividend_dts,
            dividend_yields,
            discount_curve,
            credit_spread,
            recovery_rate,
            num_steps_per_year,
        )

        elapsed = time.perf_counter() - start

        print_result(
            num_steps_per_year,
            result,
            elapsed,
        )


# ============================================================================
# 7. VALUATION WITH ZERO DIVIDEND YIELD
# ============================================================================
#
# First assume that the underlying equity pays no dividends.
#
# The valuation is repeated using progressively finer tree discretisations.
# Comparing the results provides a simple convergence check, while the
# elapsed time shows the computational cost of increasing the tree size.
# ============================================================================

dividend_yields = [0.00] * len(dividend_dts)

print("\n" + "=" * 90)
print("1. CONVERTIBLE BOND - ZERO DIVIDEND YIELD")
print("=" * 90)

run_valuations(
    dividend_yields,
    steps=[5, 10, 20],
)


# ============================================================================
# 8. VALUATION WITH 2% DIVIDEND YIELD
# ============================================================================
#
# Repeat the calculation assuming a 2% dividend yield at each dividend
# date, while leaving the other market and contractual inputs unchanged.
#
# A positive dividend yield changes the economics of holding the underlying
# stock and therefore affects the value of the embedded conversion option.
#
# A wider range of tree resolutions is used in this case to examine the
# convergence of the valuation more closely.
# ============================================================================

dividend_yields = [0.02] * len(dividend_dts)

print("\n" + "=" * 90)
print("2. CONVERTIBLE BOND - 2% DIVIDEND YIELD")
print("=" * 90)

run_valuations(
    dividend_yields,
    steps=[5, 20, 80],
)
