# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - FlatDiscountCurve
# ============================================================================




########################################################################################




########################################################################################

# ============================================================================
# 1. FIN FLAT CURVE
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN FLAT CURVE")
print("=" * 78)

curve_dt = Date(1, 1, 2019)
months = range(1, 60, 3)
dates = curve_dt.add_months(months)
print("COMPOUNDING", "DFS")
compounding = FrequencyTypes.CONTINUOUS

flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
dfs = flat_curve.df(dates)
print(compounding, dfs)

compounding = FrequencyTypes.ANNUAL
flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
dfs = flat_curve.df(dates)
print(compounding, dfs)

compounding = FrequencyTypes.SEMI_ANNUAL
flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
dfs = flat_curve.df(dates)
print(compounding, dfs)

compounding = FrequencyTypes.QUARTERLY
flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
dfs = flat_curve.df(dates)
print(compounding, dfs)

compounding = FrequencyTypes.MONTHLY
flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
dfs = flat_curve.df(dates)
print(compounding, dfs)

# =============================================================================
# 2. VISUALISE DISCOUNT FACTORS BY COMPOUNDING CONVENTION
# =============================================================================
# Discount factors decrease with maturity. The compounding convention changes
# the exact discount factor even when the quoted annual rate is the same.
month_values = list(months)
plt.figure()
for plot_compounding in [
    FrequencyTypes.CONTINUOUS,
    FrequencyTypes.ANNUAL,
    FrequencyTypes.SEMI_ANNUAL,
    FrequencyTypes.QUARTERLY,
    FrequencyTypes.MONTHLY,
]:
    plot_curve = FlatDiscountCurve(curve_dt, 0.05, plot_compounding)
    plot_dfs = plot_curve.df(dates)
    plt.plot(month_values, plot_dfs, label=str(plot_compounding))
plt.xlabel("Months from curve date")
plt.ylabel("Discount factor")
plt.title("Flat 5% curve: effect of compounding convention")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
