# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np
import matplotlib.pyplot as plt


from financepy.utils.date import Date
from financepy.market.curves.poly_discount_curve import PolyDiscountCurve

# ============================================================================
# FINANCEPY EXAMPLES - PolyDiscountCurve
# ============================================================================



# TODO
# Inherit from DiscountCurve and add df method
# Put in a convention for the rate
# Use Frequency object

PLOT_GRAPHS = False

########################################################################################




########################################################################################

# ============================================================================
# 1. FIN DISCOUNT CURVE POLYNOMIAL
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN DISCOUNT CURVE POLYNOMIAL")
print("=" * 78)

times = np.linspace(0.00, 10.0, 21)
curve_dt = Date(2, 2, 2019)
dates = curve_dt.add_years(times)
coeffs = [0.0004, -0.0001, 0.00000010]
curve1 = PolyDiscountCurve(curve_dt, coeffs)
zeros = curve1.zero_rate(dates)
fwds = curve1.fwd_rate_inst(dates)

if PLOT_GRAPHS:

    plt.figure(figsize=(6, 4))
    plt.plot(times, zeros, label="Zeros")
    plt.plot(times, fwds, label="Forwards")
    plt.xlabel("Time (years)")
    plt.ylabel("Zero Rate")
    plt.legend(loc="best")

