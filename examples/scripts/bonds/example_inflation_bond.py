# ============================================================================
# FINANCEPY EXAMPLES - InflationBond
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation and risk analytics of
# inflation-linked bonds.
#
# Two examples are considered:
#
#   1. Bloomberg US TIPS example
#   2. Quant Finance US TIPS example
#
# The examples illustrate the distinction between:
#
#   - real yield
#   - real clean and dirty prices
#   - inflation-adjusted accrued interest
#   - inflation-adjusted principal
#   - CPI indexation
#   - inflation zero curves
#   - duration and convexity
#
# Inflation-linked bonds differ from conventional fixed-rate bonds because
# their principal and coupon cash flows are linked to an inflation index.
#
# The ratio
#
#       Reference CPI / Base CPI
#
# determines the inflation adjustment applied to the bond.
# ============================================================================

import datetime as dt

import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.products.bonds.bond import YTMCalcType
from financepy.utils.format_graphs import set_plot_style

from financepy.products.bonds.bond_inflation import BondInflation
from financepy.market.curves.inflation_index_curve import (
    InflationIndexCurve,
)

from financepy.market.curves.zero_rates_discount_curve import (
    ZeroRatesDiscountCurve,
)

from financepy.market.curves.flat_discount_curve import (
    FlatDiscountCurve,
)

set_plot_style()


# ============================================================================
# GLOBAL OUTPUT FORMAT
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100


# ============================================================================
# 1. BLOOMBERG US TIPS EXAMPLE
# ============================================================================
#
# Source example:
#
# Bloomberg US TIPS example used in the original FinancePy test suite.
#
# The example first calculates conventional bond analytics using real yield.
# It then applies CPI indexation to calculate inflation-adjusted accrued
# interest and principal.
# ============================================================================

print("\n" + LINE)
print("1. BLOOMBERG US TIPS EXAMPLE")
print(LINE)


# ============================================================================
# 1.1 BOND DEFINITION
# ============================================================================

settle_dt = Date(21, 7, 2017)
issue_dt = Date(15, 7, 2010)
maturity_dt = Date(15, 7, 2020)

coupon = 0.0125

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

face = 100.0

base_cpi_value = 218.08532

ex_div_days = 0


bond = BondInflation(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
    ex_div_days,
    base_cpi_value,
)


print(f"{'Issue Date':<40}: " f"{issue_dt}")

print(f"{'Maturity Date':<40}: " f"{maturity_dt}")

print(f"{'Settlement Date':<40}: " f"{settle_dt}")

print(f"{'Coupon Rate':<40}: " f"{coupon * 100.0:12.6f}%")

print(f"{'Base CPI':<40}: " f"{base_cpi_value:12.6f}")


# ============================================================================
# 1.2 CURRENT YIELD
# ============================================================================
#
# Current yield compares annual coupon income with the clean market price.
#
# For an inflation-linked bond this is still based on the quoted real bond
# price and does not capture the full inflation adjustment.
# ============================================================================

print("\n" + LINE)
print("1.2 CURRENT YIELD")
print(LINE)

market_clean_price = 104.03502

current_yield = bond.current_yield(
    settle_dt,
    market_clean_price,
)

print(f"{'Clean Price':<40}: " f"{market_clean_price:12.6f}")

print(f"{'Current Yield':<40}: " f"{current_yield * 100.0:12.6f}%")


# ============================================================================
# 1.3 REAL YIELD TO MATURITY
# ============================================================================
#
# Inflation-linked bonds are normally quoted in terms of real yield.
#
# FinancePy supports several bond-market yield conventions. The same clean
# price is therefore converted into real YTM using:
#
#   - UK DMO
#   - US Street
#   - US Treasury
#
# Small differences can arise because the conventions treat coupon periods
# and accrued interest differently.
# ============================================================================

print("\n" + LINE)
print("1.3 REAL YIELD TO MATURITY")
print(LINE)

real_ytm_uk = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.UK_DMO,
)

real_ytm_street = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.US_STREET,
)

real_ytm_treasury = bond.yield_to_maturity(
    settle_dt,
    market_clean_price,
    YTMCalcType.US_TREASURY,
)


print(f"{'YIELD CONVENTION':<30}" f"{'REAL YTM (%)':>20}")

print(SUBLINE)

print(f"{'UK DMO':<30}" f"{real_ytm_uk * 100.0:20.6f}")

print(f"{'US STREET':<30}" f"{real_ytm_street * 100.0:20.6f}")

print(f"{'US TREASURY':<30}" f"{real_ytm_treasury * 100.0:20.6f}")


# Use the US Treasury convention for the remaining Bloomberg calculations.

real_ytm = real_ytm_treasury


# ============================================================================
# 1.4 REAL CLEAN AND DIRTY PRICE
# ============================================================================
#
# The real dirty price includes real accrued coupon interest.
#
# The real clean price removes accrued interest and corresponds to the
# quoted bond price.
# ============================================================================

print("\n" + LINE)
print("1.4 REAL CLEAN AND DIRTY PRICE")
print(LINE)

dirty_price = bond.dirty_price_from_ytm(
    settle_dt,
    real_ytm,
)

clean_price = bond.clean_price_from_ytm(
    settle_dt,
    real_ytm,
)

accrued_days = bond.accrued_days
real_accrued_interest = bond.accrued_int


print(f"{'Real Clean Price':<40}: " f"{clean_price:12.6f}")

print(f"{'Real Accrued Interest':<40}: " f"{real_accrued_interest:12.6f}")

print(f"{'Real Dirty Price':<40}: " f"{dirty_price:12.6f}")

print(f"{'Accrued Days':<40}: " f"{accrued_days}")


# ============================================================================
# 1.5 CPI INDEX RATIO
# ============================================================================
#
# The inflation index ratio compares the current reference CPI with the
# base CPI established when the bond was issued.
#
#       Index Ratio = Reference CPI / Base CPI
#
# An index ratio greater than one indicates cumulative inflation since the
# bond's base CPI date.
# ============================================================================

print("\n" + LINE)
print("1.5 CPI INDEX RATIO")
print(LINE)

ref_cpi_value = 244.65884

index_ratio = ref_cpi_value / base_cpi_value

cumulative_inflation = index_ratio - 1.0

print(f"{'Base CPI':<40}: " f"{base_cpi_value:12.6f}")

print(f"{'Reference CPI':<40}: " f"{ref_cpi_value:12.6f}")

print(f"{'Index Ratio':<40}: " f"{index_ratio:12.6f}")

print(f"{'Cumulative Inflation':<40}: " f"{cumulative_inflation * 100.0:12.6f}%")


# ============================================================================
# 1.6 INFLATION-ADJUSTED ACCRUED INTEREST
# ============================================================================
#
# Real accrued interest is scaled by the CPI index ratio to obtain the
# inflation-adjusted accrued interest.
# ============================================================================

print("\n" + LINE)
print("1.6 INFLATION-ADJUSTED ACCRUED INTEREST")
print(LINE)

inflation_accrued = bond.inflation_accrued_interest(
    settle_dt,
    face,
    ref_cpi_value,
)

print(f"{'Real Accrued Interest':<40}: " f"{real_accrued_interest:12.6f}")

print(f"{'Inflation Accrued Interest':<40}: " f"{inflation_accrued:12.6f}")


# ============================================================================
# 1.7 FLAT PRICE
# ============================================================================
#
# The inflation-linked flat price uses the CPI value associated with the
# previous coupon date.
# ============================================================================

print("\n" + LINE)
print("1.7 FLAT PRICE")
print(LINE)

last_cpn_cpi_value = 244.61839

flat_price = bond.flat_price_from_yield_to_maturity(
    settle_dt,
    real_ytm,
    last_cpn_cpi_value,
    YTMCalcType.US_TREASURY,
)

print(f"{'Last Coupon CPI':<40}: " f"{last_cpn_cpi_value:12.6f}")

print(f"{'Flat Price':<40}: " f"{flat_price:12.6f}")


# ============================================================================
# 1.8 INFLATION-ADJUSTED PRINCIPAL
# ============================================================================
#
# Inflation protection increases the effective principal value as the
# reference CPI rises relative to the base CPI.
# ============================================================================

print("\n" + LINE)
print("1.8 INFLATION-ADJUSTED PRINCIPAL")
print(LINE)

inflation_principal = bond.inflation_principal(
    settle_dt,
    face,
    real_ytm,
    ref_cpi_value,
    YTMCalcType.US_TREASURY,
)

print(f"{'Face Amount':<40}: " f"{face:12.6f}")

print(f"{'Inflation Principal':<40}: " f"{inflation_principal:12.6f}")


# ============================================================================
# 1.9 REAL-YIELD RISK MEASURES
# ============================================================================
#
# These measures describe sensitivity to changes in real yield.
# ============================================================================

print("\n" + LINE)
print("1.9 REAL-YIELD RISK MEASURES")
print(LINE)

dollar_duration = bond.dollar_duration(
    settle_dt,
    real_ytm,
)

modified_duration = bond.modified_duration(
    settle_dt,
    real_ytm,
)

macaulay_duration = bond.macaulay_duration(
    settle_dt,
    real_ytm,
)

convexity = bond.convexity_from_ytm(
    settle_dt,
    real_ytm,
)


print(f"{'MEASURE':<35}" f"{'VALUE':>20}")

print(SUBLINE)

print(f"{'Dollar Duration':<35}" f"{dollar_duration:20.6f}")

print(f"{'Modified Duration':<35}" f"{modified_duration:20.6f}")

print(f"{'Macaulay Duration':<35}" f"{macaulay_duration:20.6f}")

print(f"{'Convexity':<35}" f"{convexity:20.6f}")


# ============================================================================
# 2. QUANT FINANCE US TIPS EXAMPLE
# ============================================================================
#
# This example introduces explicit inflation-market data:
#
#   - historical CPI fixings
#   - an inflation index curve
#   - zero-coupon inflation swap rates
#   - an inflation zero curve
#   - a nominal discount curve
#
# These are the market-data components required for more complete
# inflation-linked valuation.
# ============================================================================

print("\n" + LINE)
print("2. QUANT FINANCE US TIPS EXAMPLE")
print(LINE)


# ============================================================================
# 2.1 BOND DEFINITION
# ============================================================================

settle_dt_2 = Date(23, 8, 2019)
issue_dt_2 = Date(25, 9, 2013)
maturity_dt_2 = Date(22, 3, 2068)

coupon_2 = 0.00125

freq_type_2 = FrequencyTypes.SEMI_ANNUAL
dc_type_2 = DayCountTypes.ACT_ACT_ICMA

base_cpi_value_2 = 249.70
ref_cpi_value_2 = 244.65884

ex_div_days_2 = 0


bond_2 = BondInflation(
    issue_dt_2,
    maturity_dt_2,
    coupon_2,
    freq_type_2,
    dc_type_2,
    ex_div_days_2,
    base_cpi_value_2,
)


print(f"{'Issue Date':<40}: " f"{issue_dt_2}")

print(f"{'Maturity Date':<40}: " f"{maturity_dt_2}")

print(f"{'Settlement Date':<40}: " f"{settle_dt_2}")

print(f"{'Coupon Rate':<40}: " f"{coupon_2 * 100.0:12.6f}%")

print(f"{'Base CPI':<40}: " f"{base_cpi_value_2:12.6f}")


# ============================================================================
# 2.2 NOMINAL DISCOUNT CURVE
# ============================================================================
#
# A flat nominal discount curve is used in the original example.
# ============================================================================

nominal_rate = 0.01033692

discount_curve = FlatDiscountCurve(
    settle_dt_2,
    nominal_rate,
    FrequencyTypes.ANNUAL,
    DayCountTypes.ACT_ACT_ISDA,
)

print("\n" + LINE)
print("2.2 NOMINAL DISCOUNT CURVE")
print(LINE)

print(f"{'Nominal Rate':<40}: " f"{nominal_rate * 100.0:12.6f}%")


# ============================================================================
# 2.3 CPI FIXINGS
# ============================================================================
#
# Inflation indices are published with a lag. The original example assumes
# a three-month lag.
# ============================================================================

print("\n" + LINE)
print("2.3 CPI FIXINGS")
print(LINE)

lag = 3

months = range(
    0,
    12,
    1,
)

fixing_dates = Date(
    31,
    8,
    2018,
).add_months(months)

fixing_rates = [
    284.2,
    284.1,
    284.5,
    284.6,
    285.6,
    283.0,
    285.0,
    285.1,
    288.2,
    289.2,
    289.6,
    289.5,
]


inflation_index = InflationIndexCurve(
    fixing_dates,
    fixing_rates,
    lag,
)


print(f"{'FIXING DATE':<20}" f"{'CPI':>15}")

print(SUBLINE)

for fixing_dt, fixing_rate in zip(
    fixing_dates,
    fixing_rates,
):

    print(f"{str(fixing_dt):<20}" f"{fixing_rate:15.4f}")


# ============================================================================
# 2.4 ZERO-COUPON INFLATION SWAP DATA
# ============================================================================

print("\n" + LINE)
print("2.4 ZERO-COUPON INFLATION SWAP DATA")
print(LINE)

zciis_data = [
    (Date(31, 7, 2020), 3.1500000000137085),
    (Date(31, 7, 2021), 3.547500000013759),
    (Date(31, 7, 2022), 3.675000000013573),
    (Date(31, 7, 2023), 3.7250000000134342),
    (Date(31, 7, 2024), 3.750000000013265),
    (Date(31, 7, 2025), 3.7430000000129526),
    (Date(31, 7, 2026), 3.741200000012679),
    (Date(31, 7, 2027), 3.7337000000123632),
    (Date(31, 7, 2028), 3.725000000011902),
    (Date(31, 7, 2029), 3.720000000011603),
    (Date(31, 7, 2030), 3.712517289063011),
    (Date(31, 7, 2031), 3.7013000000108764),
    (Date(31, 7, 2032), 3.686986039205209),
    (Date(31, 7, 2033), 3.671102614032895),
    (Date(31, 7, 2034), 3.655000000009778),
    (Date(31, 7, 2035), 3.6394715951305834),
    (Date(31, 7, 2036), 3.624362044800966),
    (Date(31, 7, 2037), 3.6093619727979087),
    (Date(31, 7, 2038), 3.59421438364369),
    (Date(31, 7, 2039), 3.5787000000081948),
    (Date(31, 7, 2040), 3.5626192748395624),
    (Date(31, 7, 2041), 3.545765016376823),
    (Date(31, 7, 2042), 3.527943521613608),
    (Date(31, 7, 2043), 3.508977137925462),
    (Date(31, 7, 2044), 3.48870000000685),
    (Date(31, 7, 2045), 3.467083068721011),
    (Date(31, 7, 2046), 3.4445738220594935),
    (Date(31, 7, 2047), 3.4216470902302065),
    (Date(31, 7, 2048), 3.3986861494999188),
    (Date(31, 7, 2049), 3.376000000005752),
    (Date(31, 7, 2050), 3.3538412080641233),
    (Date(31, 7, 2051), 3.3324275806807746),
    (Date(31, 7, 2052), 3.311938788306623),
    (Date(31, 7, 2053), 3.2925208131865835),
    (Date(31, 7, 2054), 3.274293040759302),
    (Date(31, 7, 2055), 3.2573541974782794),
    (Date(31, 7, 2056), 3.241787355503245),
    (Date(31, 7, 2057), 3.227664186159851),
    (Date(31, 7, 2058), 3.2150486140060774),
    (Date(31, 7, 2059), 3.204000000004159),
    (Date(31, 7, 2060), 3.1945334946674064),
    (Date(31, 7, 2061), 3.1865047145143377),
    (Date(31, 7, 2062), 3.179753073456304),
    (Date(31, 7, 2063), 3.1741427790361154),
    (Date(31, 7, 2064), 3.1695593261025223),
    (Date(31, 7, 2065), 3.1659065919088736),
    (Date(31, 7, 2066), 3.163104428386987),
    (Date(31, 7, 2067), 3.1610866681252903),
    (Date(31, 7, 2068), 3.1597994770515836),
    (Date(31, 7, 2069), 3.159200000003204),
    (Date(31, 7, 2070), 3.159242349440139),
    (Date(31, 7, 2071), 3.1598400898057433),
    (Date(31, 7, 2072), 3.16090721831932),
    (Date(31, 7, 2073), 3.162369676612098),
    (Date(31, 7, 2074), 3.1641636543027207),
]


zc_dates = []
zc_rates = []

for zc_dt, zc_rate in zciis_data:

    zc_dates.append(zc_dt)

    zc_rates.append(zc_rate / 100.0)


inflation_zero_curve = ZeroRatesDiscountCurve(
    settle_dt_2,
    zc_dates,
    zc_rates,
    FrequencyTypes.ANNUAL,
)


print(f"{'MATURITY':<20}" f"{'ZERO INFLATION (%)':>22}")

print(SUBLINE)

for zc_dt, zc_rate in zciis_data:

    print(f"{str(zc_dt):<20}" f"{zc_rate:22.6f}")


# ============================================================================
# 2.5 CURRENT YIELD
# ============================================================================

print("\n" + LINE)
print("2.5 CURRENT YIELD")
print(LINE)

market_clean_price_2 = 104.03502

current_yield_2 = bond_2.current_yield(
    settle_dt_2,
    market_clean_price_2,
)

print(f"{'Clean Price':<40}: " f"{market_clean_price_2:12.6f}")

print(f"{'Current Yield':<40}: " f"{current_yield_2 * 100.0:12.6f}%")


# ============================================================================
# 3. VISUALISE CPI FIXINGS
# ============================================================================
#
# Historical CPI fixings show the observed inflation index used by the
# inflation index curve.
# ============================================================================

print("\n" + LINE)
print("3. VISUALISE CPI FIXINGS")
print(LINE)

plot_fixing_dates = [
    dt.datetime(
        fixing_dt.y,
        fixing_dt.m,
        fixing_dt.d,
    )
    for fixing_dt in fixing_dates
]

plt.figure()

plt.plot(
    plot_fixing_dates,
    fixing_rates,
    marker="o",
)

plt.xlabel("Fixing Date")

plt.ylabel("CPI Index")

plt.title("Historical CPI Fixings")

plt.grid(True)


# ============================================================================
# 4. VISUALISE ZERO INFLATION TERM STRUCTURE
# ============================================================================
#
# The ZCIIS rates describe the market's zero-coupon inflation term structure.
#
# The shape of the curve shows how the inflation rate implied by the supplied
# market data varies by maturity.
# ============================================================================

print("\n" + LINE)
print("4. VISUALISE ZERO INFLATION TERM STRUCTURE")
print(LINE)

plot_zc_dates = [
    dt.datetime(
        zc_dt.y,
        zc_dt.m,
        zc_dt.d,
    )
    for zc_dt in zc_dates
]

plot_zc_rates = [rate * 100.0 for rate in zc_rates]

plt.figure()

plt.plot(
    plot_zc_dates,
    plot_zc_rates,
)

plt.xlabel("Maturity")

plt.ylabel("Zero Inflation Rate (%)")

plt.title("Zero-Coupon Inflation Term Structure")

plt.grid(True)


# ============================================================================
# 5. CPI INDEXATION EFFECT
# ============================================================================
#
# Show how a nominal face amount of 100 changes as the CPI index ratio changes.
#
# This isolates the central economic mechanism of an inflation-linked bond:
#
#       Indexed Principal = Face × CPI / Base CPI
# ============================================================================

print("\n" + LINE)
print("5. CPI INDEXATION EFFECT")
print(LINE)

indexed_principal = face * ref_cpi_value / base_cpi_value

print(f"{'Base CPI':<40}: " f"{base_cpi_value:12.6f}")

print(f"{'Reference CPI':<40}: " f"{ref_cpi_value:12.6f}")

print(f"{'Face Amount':<40}: " f"{face:12.6f}")

print(f"{'Simple CPI-Indexed Principal':<40}: " f"{indexed_principal:12.6f}")


# ============================================================================
# 6. SUMMARY
# ============================================================================

print("\n" + LINE)
print("6. SUMMARY")
print(LINE)

print(f"{'Bloomberg Clean Price':<40}: " f"{market_clean_price:12.6f}")

print(f"{'US Treasury Real YTM':<40}: " f"{real_ytm * 100.0:12.6f}%")

print(f"{'Real Dirty Price':<40}: " f"{dirty_price:12.6f}")

print(f"{'Real Accrued Interest':<40}: " f"{real_accrued_interest:12.6f}")

print(f"{'CPI Index Ratio':<40}: " f"{index_ratio:12.6f}")

print(f"{'Inflation Accrued Interest':<40}: " f"{inflation_accrued:12.6f}")

print(f"{'Inflation Principal':<40}: " f"{inflation_principal:12.6f}")

print(f"{'Modified Duration':<40}: " f"{modified_duration:12.6f}")

print(f"{'Convexity':<40}: " f"{convexity:12.6f}")

print("\n" + LINE)
print("END OF INFLATION BOND DEMONSTRATION")
print(LINE)


# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================

plt.show()
