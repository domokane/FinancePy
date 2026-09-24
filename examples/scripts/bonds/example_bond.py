# ============================================================================
# FINANCEPY EXAMPLES - Bond
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the principal analytics available for a
# fixed-rate bond, including:
#
#   - coupon schedules and cash flows
#   - accrued interest
#   - current yield
#   - yield to maturity
#   - clean and dirty prices
#   - duration and convexity
#   - clean/dirty price behaviour when yield is below, equal to, and above
#     the coupon rate
#   - pull-to-par and accrued-interest effects as maturity approaches
# ============================================================================

import datetime as dt

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.products.bonds.bond import Bond


# ============================================================================
# GLOBAL OUTPUT FORMAT
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100


# ============================================================================
# 1. CREATE A FIXED-RATE BOND
# ============================================================================

print("\n" + LINE)
print("1. CREATE A FIXED-RATE BOND")
print(LINE)

issue_dt = Date(15, 5, 2010)
maturity_dt = Date(15, 5, 2030)

coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

face = 100.0

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)

settle_dt = Date(15, 9, 2026)

print(f"{'Issue Date':<40}: {issue_dt}")
print(f"{'Maturity Date':<40}: {maturity_dt}")
print(f"{'Settlement Date':<40}: {settle_dt}")
print(f"{'Coupon Rate':<40}: {coupon * 100.0:.6f}%")
print(f"{'Frequency':<40}: {freq_type}")
print(f"{'Day Count':<40}: {dc_type}")


# ============================================================================
# 2. COUPON SCHEDULE
# ============================================================================

print("\n" + LINE)
print("2. COUPON SCHEDULE")
print(LINE)

print(f"{'PAYMENT DATE':<20}" f"{'COUPON FLOW':>20}")

print(SUBLINE)

coupon_flow = face * coupon / 2.0

for payment_dt in bond.cpn_dts:

    if payment_dt >= settle_dt:

        print(f"{str(payment_dt):<20}" f"{coupon_flow:20.6f}")


# ============================================================================
# 3. ACCRUED INTEREST
# ============================================================================

print("\n" + LINE)
print("3. ACCRUED INTEREST")
print(LINE)

accrued_interest = bond.accrued_interest(
    settle_dt,
    face,
)

print(f"{'Accrued Interest':<40}: " f"{accrued_interest:12.6f}")


# ============================================================================
# 4. CURRENT YIELD
# ============================================================================
#
# Current yield measures the annual coupon income relative to the bond's
# current clean market price.
#
# It does not include capital gains or losses from the bond pulling toward
# par as maturity approaches.
# ============================================================================

print("\n" + LINE)
print("4. CURRENT YIELD")
print(LINE)

market_clean_price = 96.50

current_yield = bond.current_yield(
    settle_dt,
    market_clean_price,
)

print(f"{'Market Clean Price':<40}: " f"{market_clean_price:12.6f}")

print(f"{'Current Yield':<40}: " f"{current_yield * 100.0:11.6f}%")


# ============================================================================
# 5. YIELD TO MATURITY
# ============================================================================

print("\n" + LINE)
print("5. YIELD TO MATURITY")
print(LINE)

ytm = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
)

print(f"{'Yield to Maturity':<40}: " f"{ytm * 100.0:11.6f}%")


# ============================================================================
# 6. PRICE FROM YIELD
# ============================================================================

print("\n" + LINE)
print("6. PRICE FROM YIELD")
print(LINE)

dirty_price = bond.dirty_price_from_ytm(
    settle_dt,
    ytm,
)

clean_price = bond.clean_price_from_ytm(
    settle_dt,
    ytm,
)

print(f"{'Clean Price':<40}: " f"{clean_price:12.6f}")

print(f"{'Accrued Interest':<40}: " f"{accrued_interest:12.6f}")

print(f"{'Dirty Price':<40}: " f"{dirty_price:12.6f}")

print(f"{'Clean + Accrued':<40}: " f"{clean_price + accrued_interest:12.6f}")


# ============================================================================
# 7. PRICE / YIELD ROUND-TRIP
# ============================================================================

print("\n" + LINE)
print("7. PRICE / YIELD ROUND-TRIP")
print(LINE)

round_trip_ytm = bond.yield_to_maturity(
    settle_dt,
    clean_price,
)

round_trip_price = bond.clean_price_from_ytm(
    settle_dt,
    round_trip_ytm,
)

print(f"{'Original Clean Price':<40}: " f"{clean_price:12.6f}")

print(f"{'Recovered YTM':<40}: " f"{round_trip_ytm * 100.0:11.6f}%")

print(f"{'Recovered Clean Price':<40}: " f"{round_trip_price:12.6f}")


# ============================================================================
# 8. DOLLAR DURATION
# ============================================================================

print("\n" + LINE)
print("8. DOLLAR DURATION")
print(LINE)

dollar_duration = bond.dollar_duration(
    settle_dt,
    ytm,
)

print(f"{'Dollar Duration':<40}: " f"{dollar_duration:12.6f}")


# ============================================================================
# 9. MODIFIED DURATION
# ============================================================================

print("\n" + LINE)
print("9. MODIFIED DURATION")
print(LINE)

modified_duration = bond.modified_duration(
    settle_dt,
    ytm,
)

print(f"{'Modified Duration':<40}: " f"{modified_duration:12.6f}")


# ============================================================================
# 10. MACAULAY DURATION
# ============================================================================

print("\n" + LINE)
print("10. MACAULAY DURATION")
print(LINE)

macaulay_duration = bond.macaulay_duration(
    settle_dt,
    ytm,
)

print(f"{'Macaulay Duration':<40}: " f"{macaulay_duration:12.6f}")


# ============================================================================
# 11. CONVEXITY
# ============================================================================

print("\n" + LINE)
print("11. CONVEXITY")
print(LINE)

convexity = bond.convexity_from_ytm(
    settle_dt,
    ytm,
)

print(f"{'Convexity':<40}: " f"{convexity:12.6f}")


# ============================================================================
# 12. PRICE ACTION: YIELD VERSUS COUPON
# ============================================================================
#
# For a fixed-rate bond:
#
#     YTM < Coupon  -> premium bond
#     YTM = Coupon  -> par bond
#     YTM > Coupon  -> discount bond
#
# Clean price excludes accrued interest.
#
# Dirty price includes accrued interest:
#
#     Dirty Price = Clean Price + Accrued Interest
# ============================================================================

print("\n" + LINE)
print("12. PRICE ACTION: YIELD VERSUS COUPON")
print(LINE)

yield_scenarios = [
    ("YTM < Coupon", 0.03),
    ("YTM = Coupon", coupon),
    ("YTM > Coupon", 0.07),
]

print(
    f"{'SCENARIO':<20}"
    f"{'COUPON (%)':>14}"
    f"{'YTM (%)':>14}"
    f"{'CLEAN PRICE':>16}"
    f"{'ACCRUED':>14}"
    f"{'DIRTY PRICE':>16}"
    f"{'STATUS':>14}"
)

print(SUBLINE)

for scenario_name, scenario_ytm in yield_scenarios:

    scenario_clean = bond.clean_price_from_ytm(
        settle_dt,
        scenario_ytm,
    )

    scenario_dirty = bond.dirty_price_from_ytm(
        settle_dt,
        scenario_ytm,
    )

    scenario_accrued = bond.accrued_interest(
        settle_dt,
        face,
    )

    if scenario_clean > 100.0 + 1.0e-8:
        status = "PREMIUM"

    elif scenario_clean < 100.0 - 1.0e-8:
        status = "DISCOUNT"

    else:
        status = "PAR"

    print(
        f"{scenario_name:<20}"
        f"{coupon * 100.0:14.4f}"
        f"{scenario_ytm * 100.0:14.4f}"
        f"{scenario_clean:16.6f}"
        f"{scenario_accrued:14.6f}"
        f"{scenario_dirty:16.6f}"
        f"{status:>14}"
    )


# ============================================================================
# 13. CLEAN AND DIRTY PRICE THROUGH TIME
# ============================================================================
#
# Hold YTM constant and move the settlement date toward maturity.
#
# Two effects can then be seen clearly:
#
# PULL TO PAR
#
#     YTM < Coupon:
#         The bond trades at a premium and approaches par from above.
#
#     YTM = Coupon:
#         The clean price remains around par.
#
#     YTM > Coupon:
#         The bond trades at a discount and approaches par from below.
#
# ACCRUED INTEREST
#
# Dirty price includes accrued interest. Accrued interest builds between
# coupon dates and resets when a coupon is paid.
#
# FinancePy Date objects are retained for all bond calculations. Separate
# Python datetime objects are constructed for Matplotlib because Matplotlib
# does not understand FinancePy Date objects directly.
# ============================================================================

print("\n" + LINE)
print("13. CLEAN AND DIRTY PRICE THROUGH TIME")
print(LINE)


# ============================================================================
# 13.1 GENERATE DAILY SETTLEMENT DATES
# ============================================================================

time_settle_dts = []

time_dt = settle_dt

while time_dt < maturity_dt:

    time_settle_dts.append(time_dt)

    time_dt = time_dt.add_days(1)


# ============================================================================
# 13.2 CREATE PYTHON DATES FOR MATPLOTLIB
# ============================================================================
#
# time_settle_dts:
#     FinancePy Date objects used for pricing.
#
# plot_settle_dts:
#     Standard Python datetime objects used only for plotting.
# ============================================================================

plot_settle_dts = [
    dt.datetime(
        financepy_dt.y,
        financepy_dt.m,
        financepy_dt.d,
    )
    for financepy_dt in time_settle_dts
]


# ============================================================================
# 13.3 CALCULATE CLEAN, DIRTY AND ACCRUED VALUES
# ============================================================================

time_price_results = {}

for scenario_name, scenario_ytm in yield_scenarios:

    clean_prices = []
    dirty_prices = []
    accrued_values = []

    for time_settle_dt in time_settle_dts:

        clean_price_t = bond.clean_price_from_ytm(
            time_settle_dt,
            scenario_ytm,
        )

        dirty_price_t = bond.dirty_price_from_ytm(
            time_settle_dt,
            scenario_ytm,
        )

        accrued_t = bond.accrued_interest(
            time_settle_dt,
            face,
        )

        clean_prices.append(clean_price_t)

        dirty_prices.append(dirty_price_t)

        accrued_values.append(accrued_t)

    time_price_results[scenario_name] = {
        "ytm": scenario_ytm,
        "clean": clean_prices,
        "dirty": dirty_prices,
        "accrued": accrued_values,
    }


# ============================================================================
# 13.4 SUMMARY OF START AND END PRICES
# ============================================================================

print(
    f"\n{'SCENARIO':<20}"
    f"{'YTM (%)':>12}"
    f"{'START CLEAN':>16}"
    f"{'START DIRTY':>16}"
    f"{'END CLEAN':>16}"
    f"{'END DIRTY':>16}"
)

print(SUBLINE)

for scenario_name, scenario_ytm in yield_scenarios:

    results = time_price_results[scenario_name]

    print(
        f"{scenario_name:<20}"
        f"{scenario_ytm * 100.0:12.4f}"
        f"{results['clean'][0]:16.6f}"
        f"{results['dirty'][0]:16.6f}"
        f"{results['clean'][-1]:16.6f}"
        f"{results['dirty'][-1]:16.6f}"
    )


# ============================================================================
# 14. PREMIUM BOND THROUGH TIME: YTM < COUPON
# ============================================================================
#
# Coupon > YTM.
#
# The bond trades above par because its coupon payments are more attractive
# than the market-required yield.
#
# As maturity approaches, the premium disappears and the clean price pulls
# toward par.
#
# The dirty price additionally contains accrued interest and therefore rises
# between coupon dates before dropping when the coupon is paid.
# ============================================================================

print("\n" + LINE)
print("14. PREMIUM BOND THROUGH TIME: YTM < COUPON")
print(LINE)

scenario_name = "YTM < Coupon"

results = time_price_results[scenario_name]

plt.figure()

plt.plot(
    plot_settle_dts,
    results["clean"],
    label="Clean Price",
)

plt.plot(
    plot_settle_dts,
    results["dirty"],
    label="Dirty Price",
)

plt.axhline(
    100.0,
    linestyle="--",
    label="Par",
)

plt.xlabel("Settlement Date")

plt.ylabel("Price per 100 Face")

plt.title(f"Premium Bond: " f"YTM = {results['ytm'] * 100.0:.2f}% " f"< Coupon = {coupon * 100.0:.2f}%")

plt.grid(True)
plt.legend()
plt.tight_layout()


# ============================================================================
# 15. PAR BOND THROUGH TIME: YTM = COUPON
# ============================================================================
#
# When YTM equals the coupon rate, the clean price is at or very close to par.
#
# The dirty price is not generally equal to par between coupon dates because
# it includes accrued interest.
#
# This distinction is particularly clear in the plot below.
# ============================================================================

print("\n" + LINE)
print("15. PAR BOND THROUGH TIME: YTM = COUPON")
print(LINE)

scenario_name = "YTM = Coupon"

results = time_price_results[scenario_name]

plt.figure()

plt.plot(
    plot_settle_dts,
    results["clean"],
    label="Clean Price",
)

plt.plot(
    plot_settle_dts,
    results["dirty"],
    label="Dirty Price",
)

plt.axhline(
    100.0,
    linestyle="--",
    label="Par",
)

plt.xlabel("Settlement Date")

plt.ylabel("Price per 100 Face")

plt.title(f"Par Bond: " f"YTM = Coupon = {coupon * 100.0:.2f}%")

plt.grid(True)
plt.legend()
plt.tight_layout()


# ============================================================================
# 16. DISCOUNT BOND THROUGH TIME: YTM > COUPON
# ============================================================================
#
# YTM > Coupon.
#
# The bond trades below par because its coupon is less attractive than the
# market-required yield.
#
# As maturity approaches, the discount disappears and the clean price pulls
# toward par from below.
# ============================================================================

print("\n" + LINE)
print("16. DISCOUNT BOND THROUGH TIME: YTM > COUPON")
print(LINE)

scenario_name = "YTM > Coupon"

results = time_price_results[scenario_name]

plt.figure()

plt.plot(
    plot_settle_dts,
    results["clean"],
    label="Clean Price",
)

plt.plot(
    plot_settle_dts,
    results["dirty"],
    label="Dirty Price",
)

plt.axhline(
    100.0,
    linestyle="--",
    label="Par",
)

plt.xlabel("Settlement Date")

plt.ylabel("Price per 100 Face")

plt.title(f"Discount Bond: " f"YTM = {results['ytm'] * 100.0:.2f}% " f"> Coupon = {coupon * 100.0:.2f}%")

plt.grid(True)
plt.legend()
plt.tight_layout()


# ============================================================================
# 17. ACCRUED INTEREST THROUGH TIME
# ============================================================================
#
# Accrued interest is plotted separately to show exactly where the difference
# between clean and dirty prices comes from.
#
# Accrued interest increases between coupon dates and resets around each
# coupon payment date.
#
# It does not depend on the assumed YTM, so only one series is required.
# ============================================================================

print("\n" + LINE)
print("17. ACCRUED INTEREST THROUGH TIME")
print(LINE)

results = time_price_results["YTM = Coupon"]

plt.figure()

plt.plot(
    plot_settle_dts,
    results["accrued"],
    label="Accrued Interest",
)

plt.xlabel("Settlement Date")

plt.ylabel("Accrued Interest per 100 Face")

plt.title("Bond Accrued Interest Through Time")

plt.grid(True)
plt.legend()
plt.tight_layout()


# ============================================================================
# 18. PRICE/YIELD CONVEXITY
# ============================================================================
#
# Bond price falls as yield rises, but the relationship is nonlinear.
#
# The curvature of the price/yield relationship illustrates bond convexity.
#
# The coupon rate and par price are also shown. Their intersection identifies
# the familiar par-bond case where YTM equals the coupon rate.
# ============================================================================

print("\n" + LINE)
print("18. PRICE/YIELD CONVEXITY")
print(LINE)

plot_yields = np.linspace(
    0.01,
    0.09,
    161,
)

plot_prices = [
    bond.clean_price_from_ytm(
        settle_dt,
        plot_ytm,
    )
    for plot_ytm in plot_yields
]

plt.figure()

plt.plot(
    plot_yields * 100.0,
    plot_prices,
    label="Clean Price",
)

plt.scatter(
    [ytm * 100.0],
    [market_clean_price],
    label="Market Point",
)

plt.axvline(
    coupon * 100.0,
    linestyle="--",
    label="Coupon Rate",
)

plt.axhline(
    100.0,
    linestyle="--",
    label="Par",
)

plt.xlabel("Yield to Maturity (%)")

plt.ylabel("Clean Price per 100 Face")

plt.title("Bond Price versus Yield")

plt.grid(True)
plt.legend()
plt.tight_layout()


# ============================================================================
# 19. SUMMARY
# ============================================================================

print("\n" + LINE)
print("19. SUMMARY")
print(LINE)

print(f"{'Market Clean Price':<40}: " f"{market_clean_price:12.6f}")

print(f"{'Accrued Interest':<40}: " f"{accrued_interest:12.6f}")

print(f"{'Current Yield':<40}: " f"{current_yield * 100.0:11.6f}%")

print(f"{'Yield to Maturity':<40}: " f"{ytm * 100.0:11.6f}%")

print(f"{'Dirty Price from YTM':<40}: " f"{dirty_price:12.6f}")

print(f"{'Clean Price from YTM':<40}: " f"{clean_price:12.6f}")

print(f"{'Dollar Duration':<40}: " f"{dollar_duration:12.6f}")

print(f"{'Modified Duration':<40}: " f"{modified_duration:12.6f}")

print(f"{'Macaulay Duration':<40}: " f"{macaulay_duration:12.6f}")

print(f"{'Convexity':<40}: " f"{convexity:12.6f}")

print("\n" + LINE)
print("END OF BOND DEMONSTRATION")
print(LINE)


# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================
#
# All figures are created first and displayed together at the end. This is
# convenient when the example is executed by the FinancePy run-all-examples
# script.
# ============================================================================

plt.show()
