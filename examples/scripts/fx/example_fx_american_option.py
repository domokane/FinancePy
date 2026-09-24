# ============================================================================
# FINANCEPY EXAMPLES - FXVanillaOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import numpy as np
import matplotlib.pyplot as plt

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(13, 2, 2018)
expiry_dt = Date(13, 2, 2019)

# In the FX Black-Scholes model, the FX rate is the domestic-currency
# price of one unit of foreign currency.
#
# EURUSD = 1.20 means:
#     1 EUR = 1.20 USD
#
# Domestic currency = USD
# Foreign currency  = EUR
#
# The foreign interest rate plays a role analogous to a continuous
# dividend yield in the equity Black-Scholes model.

foreign_ccy = "EUR"
domestic_ccy = "USD"
currency_pair = foreign_ccy + domestic_ccy

foreign_rate = 0.030
domestic_rate = 0.025

spot_fx_rate = 1.20
strike_fx_rate = 1.25
volatility = 0.10

notional = 1_000_000
notional_ccy = "USD"

domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)
foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)
model = BlackScholes(volatility)


# ============================================================================
# 1. EUROPEAN CALL AND PUT
# ============================================================================

print("\n" + "=" * 78)
print("1. EUROPEAN FX CALL AND PUT")
print("=" * 78)

european_call = FXVanillaOption(
    expiry_dt,
    strike_fx_rate,
    currency_pair,
    OptionTypes.EUROPEAN_CALL,
    notional,
    notional_ccy,
)

european_put = FXVanillaOption(
    expiry_dt,
    strike_fx_rate,
    currency_pair,
    OptionTypes.EUROPEAN_PUT,
    notional,
    notional_ccy,
)

call_results = european_call.value(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
)

put_results = european_put.value(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
)

call_value = call_results["v"]
put_value = put_results["v"]

print(f"Currency pair        : {currency_pair}")
print(f"Spot FX rate         : {spot_fx_rate:.6f}")
print(f"Strike FX rate       : {strike_fx_rate:.6f}")
print(f"Domestic rate        : {domestic_rate:.6f}")
print(f"Foreign rate         : {foreign_rate:.6f}")
print(f"Volatility           : {volatility:.6f}")
print(f"European call value  : {call_value:.10f}")
print(f"European put value   : {put_value:.10f}")


# ============================================================================
# 2. EUROPEAN VS AMERICAN CALL
# ============================================================================

print("\n" + "=" * 78)
print("2. EUROPEAN VS AMERICAN CALL")
print("=" * 78)

spot_fx_rates = np.arange(0.50, 2.01, 0.10)

call_european_values = []
call_american_values = []
call_early_exercise = []

print(
    f"{'SPOT':>10s}"
    f"{'EUROPEAN':>16s}"
    f"{'AMERICAN':>16s}"
    f"{'A - E':>16s}"
)
print("-" * 58)

for spot in spot_fx_rates:

    european_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_CALL,
        notional,
        notional_ccy,
    )

    american_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.AMERICAN_CALL,
        notional,
        notional_ccy,
    )

    value_european = european_option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    value_american = american_option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    premium = value_american - value_european

    call_european_values.append(value_european)
    call_american_values.append(value_american)
    call_early_exercise.append(premium)

    print(
        f"{spot:10.4f}"
        f"{value_european:16.8f}"
        f"{value_american:16.8f}"
        f"{premium:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    spot_fx_rates,
    call_european_values,
    marker="o",
    label="European Call",
)
plt.plot(
    spot_fx_rates,
    call_american_values,
    marker="o",
    label="American Call",
)
plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.2f}",
)
plt.xlabel("EUR/USD Spot Rate")
plt.ylabel("Option Value")
plt.title("FX Call: European vs American Exercise")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(
    spot_fx_rates,
    call_early_exercise,
    marker="o",
    label="American - European",
)
plt.axhline(0.0, linestyle="--")
plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.2f}",
)
plt.xlabel("EUR/USD Spot Rate")
plt.ylabel("American Value - European Value")
plt.title("FX Call: Early-Exercise Premium")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. EUROPEAN VS AMERICAN PUT
# ============================================================================

print("\n" + "=" * 78)
print("3. EUROPEAN VS AMERICAN PUT")
print("=" * 78)

put_european_values = []
put_american_values = []
put_early_exercise = []

print(
    f"{'SPOT':>10s}"
    f"{'EUROPEAN':>16s}"
    f"{'AMERICAN':>16s}"
    f"{'A - E':>16s}"
)
print("-" * 58)

for spot in spot_fx_rates:

    european_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_PUT,
        notional,
        notional_ccy,
    )

    american_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.AMERICAN_PUT,
        notional,
        notional_ccy,
    )

    value_european = european_option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    value_american = american_option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    premium = value_american - value_european

    put_european_values.append(value_european)
    put_american_values.append(value_american)
    put_early_exercise.append(premium)

    print(
        f"{spot:10.4f}"
        f"{value_european:16.8f}"
        f"{value_american:16.8f}"
        f"{premium:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    spot_fx_rates,
    put_european_values,
    marker="o",
    label="European Put",
)
plt.plot(
    spot_fx_rates,
    put_american_values,
    marker="o",
    label="American Put",
)
plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.2f}",
)
plt.xlabel("EUR/USD Spot Rate")
plt.ylabel("Option Value")
plt.title("FX Put: European vs American Exercise")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(
    spot_fx_rates,
    put_early_exercise,
    marker="o",
    label="American - European",
)
plt.axhline(0.0, linestyle="--")
plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.2f}",
)
plt.xlabel("EUR/USD Spot Rate")
plt.ylabel("American Value - European Value")
plt.title("FX Put: Early-Exercise Premium")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. EUROPEAN OPTION SPOT SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("4. EUROPEAN OPTION SPOT SENSITIVITY")
print("=" * 78)

print(
    f"{'SPOT':>10s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)
print("-" * 42)

for spot, call, put in zip(
    spot_fx_rates,
    call_european_values,
    put_european_values,
):
    print(
        f"{spot:10.4f}"
        f"{call:16.8f}"
        f"{put:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    spot_fx_rates,
    call_european_values,
    marker="o",
    label="European Call",
)
plt.plot(
    spot_fx_rates,
    put_european_values,
    marker="o",
    label="European Put",
)
plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.2f}",
)
plt.xlabel("EUR/USD Spot Rate")
plt.ylabel("Option Value")
plt.title("European FX Option Value vs Spot")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. EUROPEAN PUT-CALL PARITY
# ============================================================================

print("\n" + "=" * 78)
print("5. EUROPEAN FX PUT-CALL PARITY")
print("=" * 78)

# The standard unit-price FX put-call parity relationship is
#
#     C - P = S * exp(-r_f T) - K * exp(-r_d T)
#
# where
#
#     r_d = domestic interest rate
#     r_f = foreign interest rate
#
# This test is appropriate only if the returned "v" represents the
# corresponding unit option value. We therefore compare "v" directly
# with the unit-price parity relationship rather than dividing it by
# the contract notional.

t = (expiry_dt - value_dt) / 365.0

df_domestic = np.exp(-domestic_rate * t)
df_foreign = np.exp(-foreign_rate * t)

parity_errors = []

print(
    f"{'SPOT':>10s}"
    f"{'C - P':>16s}"
    f"{'PARITY RHS':>16s}"
    f"{'ERROR':>16s}"
)
print("-" * 58)

for spot in spot_fx_rates:

    call = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_CALL,
        notional,
        notional_ccy,
    )

    put = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_PUT,
        notional,
        notional_ccy,
    )

    call_v = call.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    put_v = put.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    lhs = call_v - put_v

    rhs = (
        spot * df_foreign
        - strike_fx_rate * df_domestic
    )

    error = lhs - rhs

    parity_errors.append(error)

    print(
        f"{spot:10.4f}"
        f"{lhs:16.10f}"
        f"{rhs:16.10f}"
        f"{error:16.10f}"
    )

max_parity_error = np.max(np.abs(parity_errors))

print()
print(
    f"Maximum absolute parity error: "
    f"{max_parity_error:.10e}"
)


# ============================================================================
# 6. VOLATILITY SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("6. VOLATILITY SENSITIVITY")
print("=" * 78)

volatilities = np.arange(0.05, 0.31, 0.025)

call_vol_values = []
put_vol_values = []

print(
    f"{'VOL':>10s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)
print("-" * 42)

for vol in volatilities:

    vol_model = BlackScholes(vol)

    call_v = european_call.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        vol_model,
    )["v"]

    put_v = european_put.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        vol_model,
    )["v"]

    call_vol_values.append(call_v)
    put_vol_values.append(put_v)

    print(
        f"{vol:10.4f}"
        f"{call_v:16.8f}"
        f"{put_v:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    volatilities,
    call_vol_values,
    marker="o",
    label="European Call",
)
plt.plot(
    volatilities,
    put_vol_values,
    marker="o",
    label="European Put",
)
plt.axvline(
    volatility,
    linestyle="--",
    label=f"Base Vol = {volatility:.2f}",
)
plt.xlabel("Volatility")
plt.ylabel("Option Value")
plt.title("European FX Option Value vs Volatility")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. STRIKE SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("7. STRIKE SENSITIVITY")
print("=" * 78)

strike_fx_rates = np.arange(0.80, 1.61, 0.05)

call_strike_values = []
put_strike_values = []

print(
    f"{'STRIKE':>10s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)
print("-" * 42)

for strike in strike_fx_rates:

    call = FXVanillaOption(
        expiry_dt,
        strike,
        currency_pair,
        OptionTypes.EUROPEAN_CALL,
        notional,
        notional_ccy,
    )

    put = FXVanillaOption(
        expiry_dt,
        strike,
        currency_pair,
        OptionTypes.EUROPEAN_PUT,
        notional,
        notional_ccy,
    )

    call_v = call.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    put_v = put.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    call_strike_values.append(call_v)
    put_strike_values.append(put_v)

    print(
        f"{strike:10.4f}"
        f"{call_v:16.8f}"
        f"{put_v:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    strike_fx_rates,
    call_strike_values,
    marker="o",
    label="European Call",
)
plt.plot(
    strike_fx_rates,
    put_strike_values,
    marker="o",
    label="European Put",
)
plt.axvline(
    spot_fx_rate,
    linestyle="--",
    label=f"Spot = {spot_fx_rate:.2f}",
)
plt.xlabel("Strike FX Rate")
plt.ylabel("Option Value")
plt.title("European FX Option Value vs Strike")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. DOMESTIC INTEREST-RATE SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("8. DOMESTIC INTEREST-RATE SENSITIVITY")
print("=" * 78)

domestic_rates = np.arange(0.00, 0.061, 0.005)

call_domestic_rate_values = []
put_domestic_rate_values = []

print(
    f"{'DOM RATE':>10s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)
print("-" * 42)

for rate in domestic_rates:

    curve = FlatDiscountCurve(value_dt, rate)

    call_v = european_call.value(
        value_dt,
        spot_fx_rate,
        curve,
        foreign_curve,
        model,
    )["v"]

    put_v = european_put.value(
        value_dt,
        spot_fx_rate,
        curve,
        foreign_curve,
        model,
    )["v"]

    call_domestic_rate_values.append(call_v)
    put_domestic_rate_values.append(put_v)

    print(
        f"{rate:10.4f}"
        f"{call_v:16.8f}"
        f"{put_v:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    domestic_rates,
    call_domestic_rate_values,
    marker="o",
    label="European Call",
)
plt.plot(
    domestic_rates,
    put_domestic_rate_values,
    marker="o",
    label="European Put",
)
plt.axvline(
    domestic_rate,
    linestyle="--",
    label=f"Base Rate = {domestic_rate:.3f}",
)
plt.xlabel("Domestic Interest Rate")
plt.ylabel("Option Value")
plt.title("European FX Option Value vs Domestic Interest Rate")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. FOREIGN INTEREST-RATE SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("9. FOREIGN INTEREST-RATE SENSITIVITY")
print("=" * 78)

foreign_rates = np.arange(0.00, 0.061, 0.005)

call_foreign_rate_values = []
put_foreign_rate_values = []

print(
    f"{'FOR RATE':>10s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)
print("-" * 42)

for rate in foreign_rates:

    curve = FlatDiscountCurve(value_dt, rate)

    call_v = european_call.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        curve,
        model,
    )["v"]

    put_v = european_put.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        curve,
        model,
    )["v"]

    call_foreign_rate_values.append(call_v)
    put_foreign_rate_values.append(put_v)

    print(
        f"{rate:10.4f}"
        f"{call_v:16.8f}"
        f"{put_v:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    foreign_rates,
    call_foreign_rate_values,
    marker="o",
    label="European Call",
)
plt.plot(
    foreign_rates,
    put_foreign_rate_values,
    marker="o",
    label="European Put",
)
plt.axvline(
    foreign_rate,
    linestyle="--",
    label=f"Base Rate = {foreign_rate:.3f}",
)
plt.xlabel("Foreign Interest Rate")
plt.ylabel("Option Value")
plt.title("European FX Option Value vs Foreign Interest Rate")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 10. NOTIONAL CONVENTION
# ============================================================================

print("\n" + "=" * 78)
print("10. NOTIONAL CONVENTION")
print("=" * 78)

# This section deliberately examines the behaviour of the returned "v"
# rather than assuming that "v" is the total cash PV of the contract.
#
# If "v" is a unit option value, changing the contract notional should
# leave it unchanged.

notionals = [
    100_000,
    500_000,
    1_000_000,
    2_000_000,
    5_000_000,
]

notional_values = []

print(
    f"{'NOTIONAL':>14s}"
    f"{'v':>20s}"
)
print("-" * 34)

for test_notional in notionals:

    option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_CALL,
        test_notional,
        notional_ccy,
    )

    results = option.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
    )

    option_value = results["v"]
    notional_values.append(option_value)

    print(
        f"{test_notional:14,.0f}"
        f"{option_value:20.10f}"
    )

notional_values = np.asarray(notional_values)

notional_range = (
    np.max(notional_values)
    - np.min(notional_values)
)

print()
print(
    f"Range of returned 'v' values : "
    f"{notional_range:.10e}"
)


# ============================================================================
# 11. VALUE() OUTPUT DIAGNOSTIC
# ============================================================================

print("\n" + "=" * 78)
print("11. FXVanillaOption.value() OUTPUT")
print("=" * 78)

# Showing the complete dictionary is useful because FXVanillaOption may
# return several representations of value in addition to "v".

results = european_call.value(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
)

for key, value in results.items():
    print(f"{key:<30s}: {value}")


# ============================================================================
# 12. AMERICAN EXERCISE DIAGNOSTICS
# ============================================================================

print("\n" + "=" * 78)
print("12. AMERICAN EXERCISE DIAGNOSTICS")
print("=" * 78)

# The theoretical American option value cannot be below the corresponding
# European value because the American holder has all European exercise
# opportunities plus the right to exercise earlier.
#
# The two prices may, however, be produced by different numerical methods.
# We therefore inspect the differences before applying a numerical
# tolerance in the sanity tests.

call_exercise_diff = (
    np.asarray(call_american_values)
    - np.asarray(call_european_values)
)

put_exercise_diff = (
    np.asarray(put_american_values)
    - np.asarray(put_european_values)
)

print(
    f"{'SPOT':>10s}"
    f"{'CALL A-E':>18s}"
    f"{'PUT A-E':>18s}"
)
print("-" * 46)

for spot, call_diff, put_diff in zip(
    spot_fx_rates,
    call_exercise_diff,
    put_exercise_diff,
):
    print(
        f"{spot:10.4f}"
        f"{call_diff:18.10f}"
        f"{put_diff:18.10f}"
    )

min_call_exercise_diff = np.min(call_exercise_diff)
min_put_exercise_diff = np.min(put_exercise_diff)

print()
print(
    f"Minimum call A-E difference : "
    f"{min_call_exercise_diff:.10e}"
)
print(
    f"Minimum put A-E difference  : "
    f"{min_put_exercise_diff:.10e}"
)


# ============================================================================
# 13. NUMERICAL SANITY TESTS
# ============================================================================

print("\n" + "=" * 78)
print("13. NUMERICAL SANITY TESTS")
print("=" * 78)

tol = 1.0e-8

# American options may be valued using a numerical method, so use a
# separate tolerance for comparisons with analytic European values.
american_tol = 5.0e-5

tests = []


# --------------------------------------------------------------------------
# Non-negative option values
# --------------------------------------------------------------------------

tests.append((
    "European call values are non-negative",
    np.min(call_european_values) >= -tol,
))

tests.append((
    "European put values are non-negative",
    np.min(put_european_values) >= -tol,
))


# --------------------------------------------------------------------------
# American options should not be worth less than European options.
# --------------------------------------------------------------------------

tests.append((
    "American call >= European call (within tolerance)",
    min_call_exercise_diff >= -american_tol,
))

tests.append((
    "American put >= European put (within tolerance)",
    min_put_exercise_diff >= -american_tol,
))


# --------------------------------------------------------------------------
# Spot monotonicity
# --------------------------------------------------------------------------

tests.append((
    "European call increases with spot",
    np.all(
        np.diff(call_european_values) >= -tol
    ),
))

tests.append((
    "European put decreases with spot",
    np.all(
        np.diff(put_european_values) <= tol
    ),
))


# --------------------------------------------------------------------------
# Positive vega
# --------------------------------------------------------------------------

tests.append((
    "European call increases with volatility",
    np.all(
        np.diff(call_vol_values) >= -tol
    ),
))

tests.append((
    "European put increases with volatility",
    np.all(
        np.diff(put_vol_values) >= -tol
    ),
))


# --------------------------------------------------------------------------
# Strike monotonicity
# --------------------------------------------------------------------------

tests.append((
    "European call decreases with strike",
    np.all(
        np.diff(call_strike_values) <= tol
    ),
))

tests.append((
    "European put increases with strike",
    np.all(
        np.diff(put_strike_values) >= -tol
    ),
))


# --------------------------------------------------------------------------
# Put-call parity
# --------------------------------------------------------------------------

tests.append((
    "European put-call parity",
    max_parity_error < 1.0e-8,
))


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

print(
    f"{'TEST':<60s}"
    f"{'RESULT':>10s}"
)
print("-" * 70)

for description, passed in tests:

    result = "PASS" if passed else "FAIL"

    print(
        f"{description:<60s}"
        f"{result:>10s}"
    )


# ============================================================================
# 14. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("14. SUMMARY")
print("=" * 78)

num_passed = sum(bool(passed) for _, passed in tests)
num_tests = len(tests)

print(f"European call value           : {call_value:.10f}")
print(f"European put value            : {put_value:.10f}")

print(
    f"Maximum put-call parity error : "
    f"{max_parity_error:.10e}"
)

print(
    f"Minimum call A-E difference   : "
    f"{min_call_exercise_diff:.10e}"
)

print(
    f"Minimum put A-E difference    : "
    f"{min_put_exercise_diff:.10e}"
)

print(
    f"Maximum call exercise premium : "
    f"{np.max(call_exercise_diff):.10f}"
)

print(
    f"Maximum put exercise premium  : "
    f"{np.max(put_exercise_diff):.10f}"
)

print(
    f"Range of 'v' across notionals : "
    f"{notional_range:.10e}"
)

print(
    f"Sanity tests passed           : "
    f"{num_passed}/{num_tests}"
)
