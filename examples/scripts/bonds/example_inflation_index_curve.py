# ============================================================================
# FINANCEPY EXAMPLES - InflationIndexCurve
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the construction and use of an inflation index
# curve from historical inflation-index observations.
#
# Inflation indices such as CPI are normally published with a reporting lag.
# The InflationIndexCurve stores:
#
#   - index observation dates
#   - corresponding index values
#   - publication lag
#
# The curve can then be used to calculate:
#
#   - the inflation index value for a reference date
#   - the corresponding index ratio
#
# The example also plots the interpolated index value and index ratio through
# time to illustrate how the curve evolves between the supplied observations.
# ============================================================================

import datetime as dt

import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.products.inflation.inflation_index_curve import (
    InflationIndexCurve,
)
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# GLOBAL OUTPUT FORMAT
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100


# ============================================================================
# 1. INFLATION INDEX CURVE
# ============================================================================
#
# Construct an inflation index curve from three historical index
# observations.
#
# A three-month lag is used. This reflects the delay between the reference
# month of an inflation index and the date on which that index becomes
# available for use in inflation-linked instruments.
# ============================================================================

print("\n" + LINE)
print("1. INFLATION INDEX CURVE")
print(LINE)


# ============================================================================
# 1.1 INDEX OBSERVATIONS
# ============================================================================

index_dates = [
    Date(15, 1, 2008),
    Date(1, 4, 2008),
    Date(1, 5, 2008),
]

index_values = [
    209.49645,
    214.823,
    216.632,
]

lag = 3


print(f"{'OBSERVATION DATE':<25}" f"{'INDEX VALUE':>20}")

print(SUBLINE)

for index_dt, index_value in zip(
    index_dates,
    index_values,
):
    print(f"{str(index_dt):<25}" f"{index_value:20.6f}")

print(SUBLINE)

print(f"{'Inflation Lag (months)':<40}: " f"{lag}")


# ============================================================================
# 1.2 BUILD THE INFLATION INDEX CURVE
# ============================================================================

curve = InflationIndexCurve(
    index_dates,
    index_values,
    lag,
)


# ============================================================================
# 2. INDEX VALUE
# ============================================================================
#
# Calculate the inflation index applicable to a specified reference date.
#
# The curve applies its lag and interpolation rules internally when
# determining the appropriate index value.
# ============================================================================

print("\n" + LINE)
print("2. INDEX VALUE")
print(LINE)

ref_date = Date(22, 7, 2008)

index_value = curve.index_value(
    ref_date,
)

print(f"{'Reference Date':<40}: " f"{ref_date}")

print(f"{'Index Value':<40}: " f"{index_value:12.6f}")


# ============================================================================
# 3. INDEX RATIO
# ============================================================================
#
# The index ratio measures the inflation adjustment relative to the base
# index level used by the curve.
#
# An index ratio:
#
#       > 1.0  indicates an increase in the index
#       = 1.0  indicates no change
#       < 1.0  indicates a decrease in the index
#
# Inflation-linked cash flows can use this ratio to scale principal and
# coupon payments.
# ============================================================================

print("\n" + LINE)
print("3. INDEX RATIO")
print(LINE)

index_ratio = curve.index_ratio(
    ref_date,
)

print(f"{'Reference Date':<40}: " f"{ref_date}")

print(f"{'Index Ratio':<40}: " f"{index_ratio:12.8f}")

print(f"{'Index Change (%)':<40}: " f"{(index_ratio - 1.0) * 100.0:12.6f}%")


# ============================================================================
# 4. INDEX VALUE AND RATIO THROUGH TIME
# ============================================================================
#
# Evaluate the curve over a sequence of reference dates.
#
# This illustrates how the inflation index used by an inflation-linked
# instrument evolves through time once the publication lag and interpolation
# convention have been applied.
# ============================================================================

print("\n" + LINE)
print("4. INDEX VALUE AND RATIO THROUGH TIME")
print(LINE)

start_date = Date(15, 4, 2008)
end_date = Date(15, 8, 2008)

plot_dates = []
plot_index_values = []
plot_index_ratios = []

current_date = start_date

while current_date <= end_date:

    value = curve.index_value(
        current_date,
    )

    ratio = curve.index_ratio(
        current_date,
    )

    plot_dates.append(current_date)

    plot_index_values.append(value)

    plot_index_ratios.append(ratio)

    current_date = current_date.add_days(1)


# ============================================================================
# 4.1 SAMPLE VALUES
# ============================================================================

print(f"{'REFERENCE DATE':<25}" f"{'INDEX VALUE':>20}" f"{'INDEX RATIO':>20}")

print(SUBLINE)

for i in range(
    0,
    len(plot_dates),
    15,
):

    print(f"{str(plot_dates[i]):<25}" f"{plot_index_values[i]:20.6f}" f"{plot_index_ratios[i]:20.8f}")


# ============================================================================
# 5. PLOT INDEX VALUE
# ============================================================================

python_dates = [
    dt.datetime(
        date.y,
        date.m,
        date.d,
    )
    for date in plot_dates
]

plt.figure()

plt.plot(
    python_dates,
    plot_index_values,
)

plt.xlabel("Reference Date")
plt.ylabel("Inflation Index")
plt.title("Inflation Index Value Through Time")
plt.grid(True)


# ============================================================================
# 6. PLOT INDEX RATIO
# ============================================================================
#
# The index ratio is often the economically important quantity because it is
# used to scale the principal and coupon cash flows of inflation-linked
# securities.
# ============================================================================

plt.figure()

plt.plot(
    python_dates,
    plot_index_ratios,
)

plt.axhline(
    1.0,
    linestyle="--",
    label="Base Index Ratio",
)

plt.xlabel("Reference Date")
plt.ylabel("Index Ratio")
plt.title("Inflation Index Ratio Through Time")
plt.legend()
plt.grid(True)


# ============================================================================
# 7. SUMMARY
# ============================================================================

print("\n" + LINE)
print("7. SUMMARY")
print(LINE)

print(f"{'Reference Date':<40}: " f"{ref_date}")

print(f"{'Inflation Lag':<40}: " f"{lag} months")

print(f"{'Index Value':<40}: " f"{index_value:12.6f}")

print(f"{'Index Ratio':<40}: " f"{index_ratio:12.8f}")

print(f"{'Index Change':<40}: " f"{(index_ratio - 1.0) * 100.0:12.6f}%")

print("\n" + LINE)
print("END OF INFLATION INDEX CURVE DEMONSTRATION")
print(LINE)


# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================

plt.show()
