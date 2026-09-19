# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time
import numpy as np


from financepy.market.curves.zero_rates_discount_curve import ZeroRatesDiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - ZeroRatesDiscountCurve
# ============================================================================




########################################################################################




########################################################################################

# ============================================================================
# 1. FIN DISCOUNT CURVE ZEROS
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN DISCOUNT CURVE ZEROS")
print("=" * 78)

dates = [
    Date(14, 9, 2016),
    Date(14, 12, 2016),
    Date(14, 6, 2017),
    Date(14, 6, 2019),
    Date(14, 6, 2021),
    Date(15, 6, 2026),
    Date(16, 6, 2031),
    Date(16, 6, 2036),
    Date(14, 6, 2046),
]

zero_rates = [
    0.006616,
    0.007049,
    0.007795,
    0.009599,
    0.011203,
    0.015068,
    0.017583,
    0.018998,
    0.020080,
]

start_dt = Date(14, 6, 2016)

times = np.linspace(0.0, 30, 100)

print(f"{'Interpolation':>28s} {'Time':>8s} {'Zero %':>11s} {'Fwd %':>11s} {'Build s':>10s}")
print("-" * 74)

for interp_type in InterpTypes:

    start = time.time()

    freq_type = FrequencyTypes.ANNUAL
    time_dc_type = DayCountTypes.ACT_ACT_ISDA

    curve = ZeroRatesDiscountCurve(
        start_dt,
        dates,
        zero_rates,
        freq_type,
        time_dc_type,
        interp_type,
    )

    zeros_cc = curve.zero_rate_cc_t(times) * 100.0
    fwd_cc = curve.fwd_rate_inst_t(times) * 100.0

    end = time.time()
    period = end - start

    for t, z, f in zip(times, zeros_cc, fwd_cc):
        print(f"{str(interp_type):>28s} {t:8.3f} {z:11.6f} {f:11.6f} {period:10.6f}")

# =============================================================================
# 2. VISUALISE ZERO AND INSTANTANEOUS FORWARD RATES
# =============================================================================
# The interpolation rule affects the shape between quoted curve nodes. Plotting
# the zero and instantaneous forward rates makes smoothness and local changes
# much easier to compare than a long table of numbers.
plt.figure()
for plot_interp_type in InterpTypes:
    plot_curve = ZeroRatesDiscountCurve(start_dt, dates, zero_rates,
                                        FrequencyTypes.ANNUAL,
                                        DayCountTypes.ACT_ACT_ISDA,
                                        plot_interp_type)
    plot_zeros = plot_curve.zero_rate_cc_t(times) * 100.0
    plt.plot(times, plot_zeros, label=str(plot_interp_type))
plt.xlabel("Time from curve date (years)")
plt.ylabel("Continuously compounded zero rate (%)")
plt.title("Zero-rate curve under alternative interpolation methods")
plt.grid(True)
plt.legend(fontsize=8)
plt.tight_layout()
plt.show()
