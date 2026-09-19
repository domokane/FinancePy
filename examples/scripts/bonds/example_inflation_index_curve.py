# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.date import Date
from financepy.products.inflation.inflation_index_curve import InflationIndexCurve

# ============================================================================
# FINANCEPY EXAMPLES - InflationIndexCurve
# ============================================================================



########################################################################################




########################################################################################

#    print(curve)

# ============================================================================
# 1. FIN INFLATION INDEX CURVE
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN INFLATION INDEX CURVE")
print("=" * 78)

index_dates = [Date(15, 1, 2008), Date(1, 4, 2008), Date(1, 5, 2008)]
index_values = [209.49645, 214.823, 216.632]
lag = 3  # months

curve = InflationIndexCurve(index_dates, index_values, lag)

ref_date = Date(22, 7, 2008)

print("LABEL", "VALUE")

value = curve.index_value(ref_date)
value = curve.index_value(ref_date)
value = curve.index_value(ref_date)
value = curve.index_value(ref_date)

print(ref_date, value)

index_ratio = curve.index_ratio(ref_date)
print(ref_date, index_ratio)

