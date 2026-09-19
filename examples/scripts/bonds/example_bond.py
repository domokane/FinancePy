"""
FinancePy Bond Class - Student Demonstration
=============================================

A single, linear Python script demonstrating the main functionality
of FinancePy's Bond class.

There are deliberately NO user-defined functions in this script.
Run it from top to bottom and follow each section in order.
"""

import numpy as np

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes

from financepy.products.bonds.bond import Bond, YTMCalcType, CouponType

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.zero_rates_discount_curve import ZeroRatesDiscountCurve
from financepy.utils.global_types import InterpTypes
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - Bond
# ============================================================================

# =============================================================================
# 0. DISPLAY HELPERS
# =============================================================================
# We are not defining separate functions, so simple formatting strings are
# reused throughout the script.

LINE = "=" * 78
SUBLINE = "-" * 78

print("\n" + LINE)
print("                 FINANCEPY BOND CLASS DEMONSTRATION")
print(LINE)
print("This script demonstrates the main calculations available on a Bond.")
print("Prices are quoted per 100 of face value unless otherwise stated.")
print(LINE)


# =============================================================================
# 1. CREATE A BOND
# =============================================================================

print("\n" + LINE)
print("1. CREATING A STANDARD FIXED-COUPON BOND")
print(LINE)

issue_dt = Date(15, 5, 2020)
maturity_dt = Date(15, 5, 2030)

coupon = 0.05  # 5% annual coupon
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
ex_div_days = 0

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
    ex_div_days,
)

settle_dt = Date(15, 9, 2026)
face = 100.0
market_clean_price = 102.50

print(f"{'Issue date':30s}: {issue_dt}")
print(f"{'Settlement date':30s}: {settle_dt}")
print(f"{'Maturity date':30s}: {maturity_dt}")
print(f"{'Coupon rate':30s}: {coupon * 100:10.4f}%")
print(f"{'Coupon frequency':30s}: {freq_type}")
print(f"{'Day-count convention':30s}: {dc_type}")
print(f"{'Face value':30s}: {face:10.2f}")
print(f"{'Observed clean price':30s}: {market_clean_price:10.4f}")


# =============================================================================
# 2. DISPLAY THE BOND OBJECT
# =============================================================================

print("\n" + LINE)
print("2. BOND OBJECT")
print(LINE)

print(bond)


# =============================================================================
# 3. PAYMENT SCHEDULE
# =============================================================================

print("\n" + LINE)
print("3. FUTURE BOND PAYMENTS")
print(LINE)
print("FinancePy can display the remaining coupon/principal cash flows.\n")

bond.print_payments(settle_dt, face)


# =============================================================================
# 4. ACCRUED INTEREST
# =============================================================================

print("\n" + LINE)
print("4. ACCRUED INTEREST")
print(LINE)

accrued_interest = bond.accrued_interest(settle_dt, face)
accrued_days = bond.accrued_days

print(f"{'Accrued days':35s}: {accrued_days:12.0f}")
print(f"{'Accrued interest':35s}: {accrued_interest:12.6f}")

print("\nAccrued interest is the coupon interest earned since the previous " "coupon date.")


# =============================================================================
# 5. CURRENT YIELD
# =============================================================================

print("\n" + LINE)
print("5. CURRENT YIELD")
print(LINE)

current_yield = bond.current_yield(settle_dt, market_clean_price)

print(f"{'Clean price':35s}: {market_clean_price:12.6f}")
print(f"{'Current yield':35s}: {current_yield * 100:12.6f}%")

print("\nCurrent yield considers the annual coupon relative to the bond's " "current market price.")


# =============================================================================
# 6. YIELD TO MATURITY
# =============================================================================

print("\n" + LINE)
print("6. YIELD TO MATURITY")
print(LINE)

ytm = bond.yield_to_maturity(settle_dt, market_clean_price)

print(f"{'Clean price':35s}: {market_clean_price:12.6f}")
print(f"{'Yield to maturity':35s}: {ytm * 100:12.6f}%")

print(
    "\nYTM is the discount rate that makes the present value of the "
    "bond's cash flows consistent with its market price."
)


# =============================================================================
# 7. DIFFERENT YTM CALCULATION CONVENTIONS
# =============================================================================

print("\n" + LINE)
print("7. YTM CALCULATION CONVENTIONS")
print(LINE)

ytm_uk = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.UK_DMO,
)

ytm_street = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.US_STREET,
)

ytm_treasury = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.US_TREASURY,
)

print(f"{'Convention':25s} {'Yield':>18s}")
print(SUBLINE)
print(f"{'UK DMO':25s} {ytm_uk * 100:17.6f}%")
print(f"{'US Street':25s} {ytm_street * 100:17.6f}%")
print(f"{'US Treasury':25s} {ytm_treasury * 100:17.6f}%")


# =============================================================================
# 8. PRICE FROM YIELD
# =============================================================================

print("\n" + LINE)
print("8. CLEAN AND DIRTY PRICE FROM YTM")
print(LINE)

dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
clean_price = bond.clean_price_from_ytm(settle_dt, ytm)

print(f"{'YTM':35s}: {ytm * 100:12.6f}%")
print(f"{'Dirty price':35s}: {dirty_price:12.6f}")
print(f"{'Clean price':35s}: {clean_price:12.6f}")
print(f"{'Accrued interest':35s}: {accrued_interest:12.6f}")

print("\nRelationship:")
print(f"  Dirty price - accrued interest = {dirty_price - accrued_interest:.6f}")
print(f"  Clean price                    = {clean_price:.6f}")


# =============================================================================
# 9. PRICE/YIELD ROUND-TRIP CHECK
# =============================================================================

print("\n" + LINE)
print("9. PRICE / YIELD ROUND-TRIP CHECK")
print(LINE)

recovered_ytm = bond.yield_to_maturity(settle_dt, clean_price)
recovered_price = bond.clean_price_from_ytm(settle_dt, recovered_ytm)

print(f"{'Original clean price':35s}: {market_clean_price:12.8f}")
print(f"{'Calculated YTM':35s}: {ytm * 100:12.8f}%")
print(f"{'Price recovered from YTM':35s}: {recovered_price:12.8f}")
print(f"{'Pricing difference':35s}: {recovered_price-market_clean_price:12.8f}")


# =============================================================================
# 10. DOLLAR DURATION
# =============================================================================

print("\n" + LINE)
print("10. DOLLAR DURATION")
print(LINE)

dollar_duration = bond.dollar_duration(settle_dt, ytm)

print(f"{'Dollar duration':35s}: {dollar_duration:12.6f}")

print("\nDollar duration measures the approximate price sensitivity of the " "bond to a change in yield.")


# =============================================================================
# 11. MODIFIED DURATION
# =============================================================================

print("\n" + LINE)
print("11. MODIFIED DURATION")
print(LINE)

modified_duration = bond.modified_duration(settle_dt, ytm)

print(f"{'Modified duration':35s}: {modified_duration:12.6f}")

print("\nModified duration measures the approximate percentage change in " "price for a small change in yield.")


# =============================================================================
# 12. MACAULAY DURATION
# =============================================================================

print("\n" + LINE)
print("12. MACAULAY DURATION")
print(LINE)

macaulay_duration = bond.macaulay_duration(settle_dt, ytm)

print(f"{'Macaulay duration':35s}: {macaulay_duration:12.6f} years")

print("\nMacaulay duration is the present-value-weighted average time at " "which the bond's cash flows are received.")


# =============================================================================
# 13. CONVEXITY
# =============================================================================

print("\n" + LINE)
print("13. CONVEXITY")
print(LINE)

convexity = bond.convexity_from_ytm(settle_dt, ytm)

print(f"{'Convexity':35s}: {convexity:12.6f}")

print("\nConvexity captures the curvature in the relationship between " "bond price and yield.")


# =============================================================================
# 14. VERIFY DURATION WITH A SMALL YIELD BUMP
# =============================================================================

print("\n" + LINE)
print("14. PRICE SENSITIVITY TO A 1 BASIS-POINT YIELD MOVE")
print(LINE)

bump = 0.0001  # 1 basis point

price_at_ytm = bond.dirty_price_from_ytm(settle_dt, ytm)
price_yield_up = bond.dirty_price_from_ytm(settle_dt, ytm + bump)
price_yield_down = bond.dirty_price_from_ytm(settle_dt, ytm - bump)

print(f"{'Yield - 1 bp':30s}: {(ytm-bump)*100:11.6f}%" f"   Price = {price_yield_down:12.6f}")
print(f"{'Original yield':30s}: {ytm*100:11.6f}%" f"   Price = {price_at_ytm:12.6f}")
print(f"{'Yield + 1 bp':30s}: {(ytm+bump)*100:11.6f}%" f"   Price = {price_yield_up:12.6f}")

print("\nNotice the inverse price/yield relationship: when yield rises, " "bond price falls.")


# =============================================================================
# 15. PRICE FROM A DISCOUNT CURVE
# =============================================================================

print("\n" + LINE)
print("15. PRICING FROM A DISCOUNT CURVE")
print(LINE)

# Construct a simple flat discount curve at the bond's YTM.
flat_curve = FlatDiscountCurve(
    settle_dt,
    ytm,
    FrequencyTypes.SEMI_ANNUAL,
)

dirty_curve_price = bond.dirty_price_from_discount_curve(
    settle_dt,
    flat_curve,
)

clean_curve_price = bond.clean_price_from_discount_curve(
    settle_dt,
    flat_curve,
)

print(f"{'Dirty price from curve':35s}: {dirty_curve_price:12.6f}")
print(f"{'Clean price from curve':35s}: {clean_curve_price:12.6f}")


# =============================================================================
# 16. OPTION-ADJUSTED SPREAD (OAS)
# =============================================================================

print("\n" + LINE)
print("16. OPTION-ADJUSTED SPREAD (OAS)")
print(LINE)

# Use a deliberately different benchmark rate so the spread is visible.
benchmark_rate = 0.035

benchmark_curve = FlatDiscountCurve(
    settle_dt,
    benchmark_rate,
    FrequencyTypes.SEMI_ANNUAL,
)

oas = bond.option_adjusted_spread(
    settle_dt,
    market_clean_price,
    benchmark_curve,
)

print(f"{'Benchmark curve rate':35s}: {benchmark_rate * 100:12.6f}%")
print(f"{'Bond clean price':35s}: {market_clean_price:12.6f}")
print(f"{'OAS':35s}: {oas * 10000:12.4f} bp")


# =============================================================================
# 17. ASSET SWAP SPREAD
# =============================================================================

print("\n" + LINE)
print("17. ASSET SWAP SPREAD")
print(LINE)

asw = bond.asset_swap_spread(
    settle_dt,
    market_clean_price,
    benchmark_curve,
)

print(f"{'Asset swap spread':35s}: {asw * 10000:12.4f} bp")
print(f"{'Option-adjusted spread':35s}: {oas * 10000:12.4f} bp")


# =============================================================================
# 18. KEY-RATE DURATIONS
# =============================================================================

print("\n" + LINE)
print("18. KEY-RATE DURATIONS")
print(LINE)

key_rate_tenors, key_rate_durations = bond.key_rate_durations(
    settle_dt,
    ytm,
)

print(f"{'Key-rate tenor':>20s} {'Duration':>20s}")
print(SUBLINE)

for tenor, krd in zip(key_rate_tenors, key_rate_durations):
    print(f"{tenor:20.2f} {krd:20.6f}")

print("\nKey-rate duration decomposes interest-rate sensitivity across " "different points on the yield curve.")


# =============================================================================
# 19. KEY-RATE DURATIONS WITH USER-SUPPLIED MARKET RATES
# =============================================================================

print("\n" + LINE)
print("19. KEY-RATE DURATIONS USING SPECIFIED MARKET RATES")
print(LINE)

market_tenors = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0])

market_rates = np.array(
    [
        0.0400,
        0.0390,
        0.0380,
        0.0375,
        0.0370,
        0.0368,
        0.0365,
    ]
)

kr_tenors, kr_durations = bond.key_rate_durations(
    settle_dt,
    ytm,
    key_rate_tenors=market_tenors,
    rates=market_rates,
)

print(f"{'Tenor':>12s} {'Market rate':>18s} {'Key-rate duration':>22s}")
print(SUBLINE)

for tenor, rate, krd in zip(kr_tenors, market_rates, kr_durations):
    print(f"{tenor:12.2f} {rate*100:17.4f}% {krd:22.6f}")


# =============================================================================
# 20. HOLDING-PERIOD RETURN / RATE OF RETURN
# =============================================================================

print("\n" + LINE)
print("20. BOND RATE OF RETURN")
print(LINE)

buy_dt = Date(15, 9, 2024)
sell_dt = Date(15, 9, 2026)

buy_ytm = 0.045
sell_ytm = 0.040

buy_price = bond.dirty_price_from_ytm(
    buy_dt,
    buy_ytm,
    YTMCalcType.US_STREET,
)

sell_price = bond.dirty_price_from_ytm(
    sell_dt,
    sell_ytm,
    YTMCalcType.US_STREET,
)

simple_return, irr, pnl = bond.calc_ror(
    buy_dt,
    sell_dt,
    buy_ytm,
    sell_ytm,
)

print(f"{'Buy date':35s}: {buy_dt}")
print(f"{'Buy YTM':35s}: {buy_ytm * 100:12.6f}%")
print(f"{'Buy dirty price':35s}: {buy_price:12.6f}")
print()
print(f"{'Sell date':35s}: {sell_dt}")
print(f"{'Sell YTM':35s}: {sell_ytm * 100:12.6f}%")
print(f"{'Sell dirty price':35s}: {sell_price:12.6f}")
print()
print(f"{'Simple return':35s}: {simple_return * 100:12.6f}%")
print(f"{'IRR':35s}: {irr * 100:12.6f}%")
print(f"{'P&L':35s}: {pnl:12.6f}")


# =============================================================================
# 21. EX-DIVIDEND BOND
# =============================================================================

print("\n" + LINE)
print("21. EX-DIVIDEND PERIOD")
print(LINE)

ex_div_bond = Bond(
    Date(7, 9, 2020),
    Date(7, 9, 2030),
    0.05,
    FrequencyTypes.SEMI_ANNUAL,
    DayCountTypes.ACT_ACT_ICMA,
    7,  # 7 ex-dividend days
)

ex_div_ytm = 0.05
ex_div_face = 100.0

print("The following shows how accrued interest and price behave as " "settlement approaches a coupon date.\n")

print(f"{'Settlement':>15s}" f"{'Dirty Price':>18s}" f"{'Accrued':>18s}" f"{'Clean Price':>18s}")
print(SUBLINE)

ex_settle_dt = Date(25, 8, 2026)

for _ in range(13):
    ex_settle_dt = ex_settle_dt.add_days(1)

    ex_accrued = ex_div_bond.accrued_interest(
        ex_settle_dt,
        ex_div_face,
    )

    ex_dirty = ex_div_bond.dirty_price_from_ytm(
        ex_settle_dt,
        ex_div_ytm,
    )

    ex_clean = ex_dirty - ex_accrued

    print(f"{str(ex_settle_dt):>15s}" f"{ex_dirty:18.6f}" f"{ex_accrued:18.6f}" f"{ex_clean:18.6f}")


# =============================================================================
# 22. CUSTOM / MANUAL CASH-FLOW SCHEDULE
# =============================================================================

print("\n" + LINE)
print("22. RESETTING THE BOND'S CASH FLOWS MANUALLY")
print(LINE)

custom_bond = Bond(
    issue_dt=Date(1, 1, 2025),
    maturity_dt=Date(1, 1, 2028),
    coupon=0.05,
    freq_type=FrequencyTypes.ANNUAL,
    accrual_dc_type=DayCountTypes.ACT_ACT_ISDA,
    cal_type=CalendarTypes.UNITED_STATES,
)

custom_settle_dt = Date(1, 4, 2025)
custom_ytm = 0.05

print("Automatically generated schedule:")
custom_bond.print_payments(custom_settle_dt, 100.0)

# Coupon dates.
coupon_dates = [
    Date(1, 1, 2025),
    Date(1, 1, 2026),
    Date(1, 1, 2027),
    Date(1, 1, 2028),
]

# Actual payment dates.
payment_dates = [
    Date(1, 1, 2025),
    Date(1, 1, 2026),
    Date(1, 1, 2027),
    Date(1, 1, 2028),
]

# Cash flows are expressed per unit of face value.
flow_amounts = np.array(
    [
        0.00,
        0.05,
        0.05,
        1.05,  # final coupon + principal
    ]
)

custom_bond.reset_flows(
    coupon_dates,
    payment_dates,
    flow_amounts,
)

print("\nAfter reset_flows():")
custom_bond.print_payments(custom_settle_dt, 100.0)

custom_accrued = custom_bond.accrued_interest(
    custom_settle_dt,
    100.0,
)

custom_dirty = custom_bond.dirty_price_from_ytm(
    custom_settle_dt,
    custom_ytm,
)

custom_clean = custom_dirty - custom_accrued

print()
print(f"{'Dirty price':35s}: {custom_dirty:12.6f}")
print(f"{'Accrued interest':35s}: {custom_accrued:12.6f}")
print(f"{'Clean price':35s}: {custom_clean:12.6f}")


# =============================================================================
# 23. FIXED VS ACCRUED COUPON TYPES
# =============================================================================

print("\n" + LINE)
print("23. COUPON TYPES: FIXED VS ACCRUED")
print(LINE)

short_issue_dt = Date(31, 3, 2022)
short_maturity_dt = Date(31, 7, 2023)
short_settle_dt = Date(1, 5, 2022)

short_coupon = 0.0275
short_freq = FrequencyTypes.SEMI_ANNUAL
short_dc = DayCountTypes.ACT_360
short_face = 1_000_000.0

print("\nA) CouponType.FIXED")
print(SUBLINE)

fixed_coupon_bond = Bond(
    short_issue_dt,
    short_maturity_dt,
    short_coupon,
    short_freq,
    short_dc,
    cpn_type=CouponType.FIXED,
)

fixed_coupon_bond.print_payments(
    short_settle_dt,
    short_face,
)

print("\nB) CouponType.ACCRUED")
print(SUBLINE)

accrued_coupon_bond = Bond(
    short_issue_dt,
    short_maturity_dt,
    short_coupon,
    short_freq,
    short_dc,
    cpn_type=CouponType.ACCRUED,
)

accrued_coupon_bond.print_payments(
    short_settle_dt,
    short_face,
)

print("\nThe ACCRUED coupon type allows irregular first/last coupons to " "reflect the actual accrual period.")


# =============================================================================
# 24. NON-FLAT ZERO-RATE CURVE
# =============================================================================

print("\n" + LINE)
print("24. PRICING USING A NON-FLAT ZERO-RATE CURVE")
print(LINE)

curve_value_dt = settle_dt

spot_dates = [
    Date(15, 9, 2027),
    Date(15, 9, 2028),
    Date(15, 9, 2030),
]

spot_rates = [
    0.0300,
    0.0350,
    0.0400,
]

zero_curve = ZeroRatesDiscountCurve(
    curve_value_dt,
    spot_dates,
    spot_rates,
    FrequencyTypes.SEMI_ANNUAL,
    DayCountTypes.ACT_360,
    InterpTypes.LINEAR_ZERO_RATES,
)

zero_dirty = bond.dirty_price_from_discount_curve(
    settle_dt,
    zero_curve,
)

zero_clean = bond.clean_price_from_discount_curve(
    settle_dt,
    zero_curve,
)

print(f"{'Curve point':>20s} {'Zero rate':>20s}")
print(SUBLINE)

for date, rate in zip(spot_dates, spot_rates):
    print(f"{str(date):>20s} {rate*100:19.4f}%")

print()
print(f"{'Dirty price from zero curve':35s}: {zero_dirty:12.6f}")
print(f"{'Clean price from zero curve':35s}: {zero_clean:12.6f}")


# =============================================================================
# 25. SUMMARY
# =============================================================================

print("\n" + LINE)
print("25. SUMMARY OF THE MAIN BOND RESULTS")
print(LINE)

print(f"{'Market clean price':40s}: {market_clean_price:12.6f}")
print(f"{'Accrued interest':40s}: {accrued_interest:12.6f}")
print(f"{'Current yield':40s}: {current_yield * 100:11.6f}%")
print(f"{'Yield to maturity':40s}: {ytm * 100:11.6f}%")
print(f"{'Dirty price from YTM':40s}: {dirty_price:12.6f}")
print(f"{'Clean price from YTM':40s}: {clean_price:12.6f}")
print(f"{'Dollar duration':40s}: {dollar_duration:12.6f}")
print(f"{'Modified duration':40s}: {modified_duration:12.6f}")
print(f"{'Macaulay duration':40s}: {macaulay_duration:12.6f}")
print(f"{'Convexity':40s}: {convexity:12.6f}")
print(f"{'Asset swap spread':40s}: {asw * 10000:11.4f} bp")
print(f"{'Option-adjusted spread':40s}: {oas * 10000:11.4f} bp")

print("\n" + LINE)
print("                     END OF BOND DEMONSTRATION")
print(LINE)

# =============================================================================
# 26. VISUALISE PRICE/YIELD CONVEXITY
# =============================================================================
# A bond's price falls as its yield rises, but the relationship is curved rather
# than linear. This plot makes both the inverse relationship and convexity easy
# to see. The marked point is the market price/yield used above.
plot_yields = np.linspace(max(0.0001, ytm - 0.03), ytm + 0.03, 61)
plot_prices = [bond.clean_price_from_ytm(settle_dt, y) for y in plot_yields]

plt.figure()
plt.plot(plot_yields * 100.0, plot_prices, label="Clean price")
plt.scatter([ytm * 100.0], [market_clean_price], label="Market point")
plt.xlabel("Yield to maturity (%)")
plt.ylabel("Clean price per 100 face")
plt.title("Bond price versus yield")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Key-rate durations show where on the yield curve the bond is most sensitive.
plt.figure()
plt.bar(key_rate_tenors, key_rate_durations)
plt.xlabel("Key-rate tenor (years)")
plt.ylabel("Key-rate duration")
plt.title("Bond key-rate duration profile")
plt.grid(True, axis="y")
plt.tight_layout()
plt.show()
