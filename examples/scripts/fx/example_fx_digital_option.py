# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - FXDigitalOption
# ============================================================================
#
# This example demonstrates:
#
# 1. Baseline domestic-currency digital call
# 2. Analytic versus Monte Carlo for all digital types
# 3. Domestic digital call and put versus spot
# 4. Monte Carlo versus analytic across spot
# 5. Monte Carlo path convergence
# 6. Call / put complementarity
# 7. Strike sensitivity
# 8. Volatility sensitivity
# 9. Notional scaling
# 10. Domestic versus foreign payout conventions
# 11. Numerical sanity tests
#
# For currency pair FORDOM:
#
#     S = units of domestic currency per unit of foreign currency
#
# Example:
#
#     EURUSD = 1.20
#
# means:
#
#     1 EUR = 1.20 USD
#
# Domestic-currency digital:
#
#     Call = N * DF_dom * N(d2)
#     Put  = N * DF_dom * N(-d2)
#
# Foreign-currency digital:
#
#     Call = N * S * DF_for * N(d1)
#     Put  = N * S * DF_for * N(-d1)
#
# Hence:
#
#     Domestic Call + Domestic Put = N * DF_dom
#
#     Foreign Call + Foreign Put   = N * S * DF_for
#
# The Monte Carlo implementation simulates only the terminal FX rate because
# a European digital option is path independent.
# ============================================================================


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(13, 2, 2018)
expiry_dt = Date(13, 2, 2019)

ccy1 = "EUR"
ccy2 = "USD"

currency_pair = ccy1 + ccy2

foreign_rate = 0.030
domestic_rate = 0.025

spot_fx_rate = 1.20
strike_fx_rate = 1.25
volatility = 0.10

notional = 1.0

domestic_curve = FlatDiscountCurve(
    value_dt,
    domestic_rate,
)

foreign_curve = FlatDiscountCurve(
    value_dt,
    foreign_rate,
)

model = BlackScholes(volatility)

seed = 4242
num_paths = 100000


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def make_option(
    option_type,
    strike,
    premium_currency,
    option_notional=1.0,
):
    return FXDigitalOption(
        expiry_dt,
        strike,
        currency_pair,
        option_type,
        option_notional,
        premium_currency,
    )


def analytic_value(
    option_type,
    spot,
    strike,
    premium_currency,
    option_notional=1.0,
    option_model=None,
):
    if option_model is None:
        option_model = model

    option = make_option(
        option_type,
        strike,
        premium_currency,
        option_notional,
    )

    return option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        option_model,
    )


def mc_value(
    option_type,
    spot,
    strike,
    premium_currency,
    option_notional=1.0,
    option_model=None,
    paths=num_paths,
    mc_seed=seed,
):
    if option_model is None:
        option_model = model

    option = make_option(
        option_type,
        strike,
        premium_currency,
        option_notional,
    )

    return option.value_mc(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        option_model,
        num_paths=paths,
        seed=mc_seed,
    )


def print_test(name, condition):
    result = "PASS" if condition else "FAIL"
    print(f"{name:<68s}{result:>8s}")
    return bool(condition)


# ============================================================================
# 1. BASELINE DOMESTIC-CURRENCY DIGITAL CALL
# ============================================================================

print("\n" + "=" * 78)
print("1. BASELINE DOMESTIC-CURRENCY DIGITAL CALL")
print("=" * 78)

digital_call = make_option(
    OptionTypes.DIGITAL_CALL,
    strike_fx_rate,
    ccy2,
    notional,
)

analytic_call = digital_call.value(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
)

mc_call = digital_call.value_mc(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
    num_paths=num_paths,
    seed=seed,
)

print(f"Currency pair       : {currency_pair}")
print(f"Spot FX rate        : {spot_fx_rate:.8f}")
print(f"Strike FX rate      : {strike_fx_rate:.8f}")
print(f"Volatility          : {volatility:.8f}")
print(f"Domestic rate       : {domestic_rate:.8f}")
print(f"Foreign rate        : {foreign_rate:.8f}")
print(f"Notional            : {notional:.8f}")
print(f"Premium currency    : {ccy2}")
print(f"Monte Carlo paths   : {num_paths}")
print()
print(f"Analytic value      : {analytic_call:.8f}")
print(f"Monte Carlo value   : {mc_call:.8f}")
print(f"MC - analytic       : {mc_call - analytic_call:.8f}")


# ============================================================================
# 2. ANALYTIC VERSUS MONTE CARLO - ALL FOUR CASES
# ============================================================================

print("\n" + "=" * 78)
print("2. ANALYTIC VERSUS MONTE CARLO - ALL FOUR CASES")
print("=" * 78)

cases = [
    (OptionTypes.DIGITAL_CALL, ccy2),
    (OptionTypes.DIGITAL_PUT, ccy2),
    (OptionTypes.DIGITAL_CALL, ccy1),
    (OptionTypes.DIGITAL_PUT, ccy1),
]

all_type_results = {}

print(
    f"{'TYPE':<28s}"
    f"{'CCY':>8s}"
    f"{'ANALYTIC':>16s}"
    f"{'MC VALUE':>16s}"
    f"{'DIFF':>16s}"
)

print("-" * 84)

for option_type, premium_currency in cases:

    option = make_option(
        option_type,
        strike_fx_rate,
        premium_currency,
        notional,
    )

    value = option.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
    )

    value_mc = option.value_mc(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
        num_paths=num_paths,
        seed=seed,
    )

    diff = value_mc - value

    all_type_results[
        (option_type, premium_currency)
    ] = (value, value_mc)

    print(
        f"{str(option_type):<28s}"
        f"{premium_currency:>8s}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{diff:16.8f}"
    )


# ============================================================================
# 3. DOMESTIC DIGITAL CALL AND PUT VERSUS SPOT
# ============================================================================

print("\n" + "=" * 78)
print("3. DOMESTIC DIGITAL CALL AND PUT VERSUS SPOT")
print("=" * 78)

spot_fx_rates = np.linspace(
    0.75,
    1.75,
    41,
)

domestic_call = make_option(
    OptionTypes.DIGITAL_CALL,
    strike_fx_rate,
    ccy2,
)

domestic_put = make_option(
    OptionTypes.DIGITAL_PUT,
    strike_fx_rate,
    ccy2,
)

domestic_call_values = domestic_call.value(
    value_dt,
    spot_fx_rates,
    domestic_curve,
    foreign_curve,
    model,
)

domestic_put_values = domestic_put.value(
    value_dt,
    spot_fx_rates,
    domestic_curve,
    foreign_curve,
    model,
)

print(
    f"{'SPOT':>14s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
    f"{'SUM':>16s}"
)

print("-" * 62)

for spot, call, put in zip(
    spot_fx_rates,
    domestic_call_values,
    domestic_put_values,
):
    print(
        f"{spot:14.8f}"
        f"{call:16.8f}"
        f"{put:16.8f}"
        f"{call + put:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    spot_fx_rates,
    domestic_call_values,
    label="Digital Call",
)

plt.plot(
    spot_fx_rates,
    domestic_put_values,
    label="Digital Put",
)

plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Spot FX Rate")
plt.ylabel("Option Value")
plt.title("Domestic-Currency FX Digitals versus Spot")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. ANALYTIC VERSUS MONTE CARLO ACROSS SPOT
# ============================================================================

print("\n" + "=" * 78)
print("4. ANALYTIC VERSUS MONTE CARLO ACROSS SPOT")
print("=" * 78)

mc_spots = np.linspace(
    0.90,
    1.60,
    15,
)

spot_analytic_values = []
spot_mc_values = []
spot_mc_errors = []

print(
    f"{'SPOT':>14s}"
    f"{'ANALYTIC':>16s}"
    f"{'MC VALUE':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 62)

for spot in mc_spots:

    value = analytic_value(
        OptionTypes.DIGITAL_CALL,
        spot,
        strike_fx_rate,
        ccy2,
    )

    value_mc = mc_value(
        OptionTypes.DIGITAL_CALL,
        spot,
        strike_fx_rate,
        ccy2,
        paths=num_paths,
    )

    error = value_mc - value

    spot_analytic_values.append(value)
    spot_mc_values.append(value_mc)
    spot_mc_errors.append(error)

    print(
        f"{spot:14.8f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{error:16.8f}"
    )

spot_analytic_values = np.asarray(
    spot_analytic_values
)

spot_mc_values = np.asarray(
    spot_mc_values
)

spot_mc_errors = np.asarray(
    spot_mc_errors
)

plt.figure(figsize=(10, 6))

plt.plot(
    mc_spots,
    spot_analytic_values,
    label="Analytic",
)

plt.plot(
    mc_spots,
    spot_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Spot FX Rate")
plt.ylabel("Digital Call Value")
plt.title("FX Digital Call: Analytic versus Monte Carlo")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    mc_spots,
    spot_mc_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Spot FX Rate")
plt.ylabel("MC Value - Analytic Value")
plt.title("FX Digital Call: Monte Carlo Error versus Spot")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. MONTE CARLO PATH CONVERGENCE
# ============================================================================

print("\n" + "=" * 78)
print("5. MONTE CARLO PATH CONVERGENCE")
print("=" * 78)

num_paths_list = [
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
    200000,
]

path_mc_values = []
path_errors = []

analytic_baseline = analytic_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
)

print(
    f"{'PATHS':>12s}"
    f"{'MC VALUE':>16s}"
    f"{'ANALYTIC':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 60)

for paths in num_paths_list:

    value_mc = mc_value(
        OptionTypes.DIGITAL_CALL,
        spot_fx_rate,
        strike_fx_rate,
        ccy2,
        paths=paths,
    )

    error = value_mc - analytic_baseline

    path_mc_values.append(value_mc)
    path_errors.append(error)

    print(
        f"{paths:12d}"
        f"{value_mc:16.8f}"
        f"{analytic_baseline:16.8f}"
        f"{error:16.8f}"
    )

path_mc_values = np.asarray(
    path_mc_values
)

path_errors = np.asarray(
    path_errors
)

plt.figure(figsize=(10, 6))

plt.semilogx(
    num_paths_list,
    path_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.axhline(
    analytic_baseline,
    linestyle="--",
    label="Analytic",
)

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Digital Call Value")
plt.title("FX Digital Call: Monte Carlo Path Convergence")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.semilogx(
    num_paths_list,
    path_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("MC Value - Analytic Value")
plt.title("FX Digital Call: Monte Carlo Error")
plt.grid(True)
plt.show()


# ============================================================================
# 6. CALL / PUT COMPLEMENTARITY
# ============================================================================

print("\n" + "=" * 78)
print("6. CALL / PUT COMPLEMENTARITY")
print("=" * 78)

t = (expiry_dt - value_dt) / 365.0

domestic_df = domestic_curve.df_t(t)
foreign_df = foreign_curve.df_t(t)

domestic_sum = (
    domestic_call_values
    + domestic_put_values
)

domestic_expected = (
    notional
    * domestic_df
)

domestic_errors = (
    domestic_sum
    - domestic_expected
)

foreign_call = make_option(
    OptionTypes.DIGITAL_CALL,
    strike_fx_rate,
    ccy1,
)

foreign_put = make_option(
    OptionTypes.DIGITAL_PUT,
    strike_fx_rate,
    ccy1,
)

foreign_call_values = foreign_call.value(
    value_dt,
    spot_fx_rates,
    domestic_curve,
    foreign_curve,
    model,
)

foreign_put_values = foreign_put.value(
    value_dt,
    spot_fx_rates,
    domestic_curve,
    foreign_curve,
    model,
)

foreign_sum = (
    foreign_call_values
    + foreign_put_values
)

foreign_expected = (
    notional
    * spot_fx_rates
    * foreign_df
)

foreign_errors = (
    foreign_sum
    - foreign_expected
)

print(
    f"{'SPOT':>14s}"
    f"{'DOM ERROR':>20s}"
    f"{'FOR ERROR':>20s}"
)

print("-" * 54)

for spot, dom_error, for_error in zip(
    spot_fx_rates,
    domestic_errors,
    foreign_errors,
):
    print(
        f"{spot:14.8f}"
        f"{dom_error:20.10e}"
        f"{for_error:20.10e}"
    )


# ============================================================================
# 7. MONTE CARLO CALL / PUT COMPLEMENTARITY
# ============================================================================

print("\n" + "=" * 78)
print("7. MONTE CARLO CALL / PUT COMPLEMENTARITY")
print("=" * 78)

mc_dom_call = mc_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    paths=num_paths,
)

mc_dom_put = mc_value(
    OptionTypes.DIGITAL_PUT,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    paths=num_paths,
)

mc_for_call = mc_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy1,
    paths=num_paths,
)

mc_for_put = mc_value(
    OptionTypes.DIGITAL_PUT,
    spot_fx_rate,
    strike_fx_rate,
    ccy1,
    paths=num_paths,
)

mc_dom_expected = (
    notional
    * domestic_df
)

mc_for_expected = (
    notional
    * spot_fx_rate
    * foreign_df
)

print(
    f"Domestic MC call              : "
    f"{mc_dom_call:.8f}"
)

print(
    f"Domestic MC put               : "
    f"{mc_dom_put:.8f}"
)

print(
    f"Domestic MC call + put        : "
    f"{mc_dom_call + mc_dom_put:.8f}"
)

print(
    f"Domestic expected             : "
    f"{mc_dom_expected:.8f}"
)

print(
    f"Domestic MC parity error      : "
    f"{mc_dom_call + mc_dom_put - mc_dom_expected:.10e}"
)

print()

print(
    f"Foreign MC call               : "
    f"{mc_for_call:.8f}"
)

print(
    f"Foreign MC put                : "
    f"{mc_for_put:.8f}"
)

print(
    f"Foreign MC call + put         : "
    f"{mc_for_call + mc_for_put:.8f}"
)

print(
    f"Foreign expected              : "
    f"{mc_for_expected:.8f}"
)

print(
    f"Foreign MC parity error       : "
    f"{mc_for_call + mc_for_put - mc_for_expected:.10e}"
)


# ============================================================================
# 8. STRIKE SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("8. STRIKE SENSITIVITY")
print("=" * 78)

strike_fx_rates = np.linspace(
    0.80,
    1.70,
    37,
)

strike_call_values = []
strike_put_values = []

print(
    f"{'STRIKE':>14s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)

print("-" * 46)

for strike in strike_fx_rates:

    call_value = analytic_value(
        OptionTypes.DIGITAL_CALL,
        spot_fx_rate,
        strike,
        ccy2,
    )

    put_value = analytic_value(
        OptionTypes.DIGITAL_PUT,
        spot_fx_rate,
        strike,
        ccy2,
    )

    strike_call_values.append(call_value)
    strike_put_values.append(put_value)

    print(
        f"{strike:14.8f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

strike_call_values = np.asarray(
    strike_call_values
)

strike_put_values = np.asarray(
    strike_put_values
)

plt.figure(figsize=(10, 6))

plt.plot(
    strike_fx_rates,
    strike_call_values,
    label="Digital Call",
)

plt.plot(
    strike_fx_rates,
    strike_put_values,
    label="Digital Put",
)

plt.axvline(
    spot_fx_rate,
    linestyle="--",
    label="Spot",
)

plt.xlabel("Strike FX Rate")
plt.ylabel("Option Value")
plt.title("Domestic-Currency FX Digitals versus Strike")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. VOLATILITY SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("9. VOLATILITY SENSITIVITY")
print("=" * 78)

volatilities = np.linspace(
    0.02,
    0.40,
    20,
)

vol_call_values = []
vol_put_values = []

print(
    f"{'VOL':>14s}"
    f"{'CALL':>16s}"
    f"{'PUT':>16s}"
)

print("-" * 46)

for vol in volatilities:

    vol_model = BlackScholes(vol)

    call_value = analytic_value(
        OptionTypes.DIGITAL_CALL,
        spot_fx_rate,
        strike_fx_rate,
        ccy2,
        option_model=vol_model,
    )

    put_value = analytic_value(
        OptionTypes.DIGITAL_PUT,
        spot_fx_rate,
        strike_fx_rate,
        ccy2,
        option_model=vol_model,
    )

    vol_call_values.append(call_value)
    vol_put_values.append(put_value)

    print(
        f"{vol:14.8f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

vol_call_values = np.asarray(
    vol_call_values
)

vol_put_values = np.asarray(
    vol_put_values
)

plt.figure(figsize=(10, 6))

plt.plot(
    volatilities,
    vol_call_values,
    label="Digital Call",
)

plt.plot(
    volatilities,
    vol_put_values,
    label="Digital Put",
)

plt.xlabel("Volatility")
plt.ylabel("Option Value")
plt.title("Domestic-Currency FX Digitals versus Volatility")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 10. NOTIONAL SCALING
# ============================================================================

print("\n" + "=" * 78)
print("10. NOTIONAL SCALING")
print("=" * 78)

notionals = [
    1.0,
    10.0,
    100.0,
    1000.0,
    1_000_000.0,
]

print(
    f"{'NOTIONAL':>16s}"
    f"{'ANALYTIC':>20s}"
    f"{'MC VALUE':>20s}"
)

print("-" * 56)

for test_notional in notionals:

    value = analytic_value(
        OptionTypes.DIGITAL_CALL,
        spot_fx_rate,
        strike_fx_rate,
        ccy2,
        option_notional=test_notional,
    )

    value_mc = mc_value(
        OptionTypes.DIGITAL_CALL,
        spot_fx_rate,
        strike_fx_rate,
        ccy2,
        option_notional=test_notional,
        paths=num_paths,
    )

    print(
        f"{test_notional:16.8f}"
        f"{value:20.8f}"
        f"{value_mc:20.8f}"
    )


# ============================================================================
# 11. DOMESTIC VERSUS FOREIGN PAYOUT CONVENTIONS
# ============================================================================

print("\n" + "=" * 78)
print("11. DOMESTIC VERSUS FOREIGN PAYOUT CONVENTIONS")
print("=" * 78)

print(
    f"{'SPOT':>14s}"
    f"{'DOM CALL':>16s}"
    f"{'FOR CALL':>16s}"
    f"{'DOM PUT':>16s}"
    f"{'FOR PUT':>16s}"
)

print("-" * 78)

for i, spot in enumerate(spot_fx_rates):

    print(
        f"{spot:14.8f}"
        f"{domestic_call_values[i]:16.8f}"
        f"{foreign_call_values[i]:16.8f}"
        f"{domestic_put_values[i]:16.8f}"
        f"{foreign_put_values[i]:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    spot_fx_rates,
    domestic_call_values,
    label="Domestic Digital Call",
)

plt.plot(
    spot_fx_rates,
    foreign_call_values,
    label="Foreign Digital Call",
)

plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Spot FX Rate")
plt.ylabel("Option Value")
plt.title("Domestic versus Foreign Digital Calls")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 12. NUMERICAL SANITY TESTS
# ============================================================================

print("\n" + "=" * 78)
print("12. NUMERICAL SANITY TESTS")
print("=" * 78)

tol = 1.0e-10

tests = []

tests.append((
    "Domestic digital call values are non-negative",
    np.all(domestic_call_values >= -tol),
))

tests.append((
    "Domestic digital put values are non-negative",
    np.all(domestic_put_values >= -tol),
))

tests.append((
    "Foreign digital call values are non-negative",
    np.all(foreign_call_values >= -tol),
))

tests.append((
    "Foreign digital put values are non-negative",
    np.all(foreign_put_values >= -tol),
))

tests.append((
    "Domestic digital call increases with spot",
    np.all(np.diff(domestic_call_values) >= -tol),
))

tests.append((
    "Domestic digital put decreases with spot",
    np.all(np.diff(domestic_put_values) <= tol),
))

tests.append((
    "Domestic digital call decreases with strike",
    np.all(np.diff(strike_call_values) <= tol),
))

tests.append((
    "Domestic digital put increases with strike",
    np.all(np.diff(strike_put_values) >= -tol),
))

tests.append((
    "Domestic analytic call + put identity",
    np.max(np.abs(domestic_errors)) < tol,
))

tests.append((
    "Foreign analytic call + put identity",
    np.max(np.abs(foreign_errors)) < tol,
))

unit_analytic = analytic_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    option_notional=1.0,
)

scaled_analytic = analytic_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    option_notional=1000.0,
)

tests.append((
    "Analytic value scales linearly with notional",
    abs(
        scaled_analytic
        - 1000.0 * unit_analytic
    ) < tol,
))

unit_mc = mc_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    option_notional=1.0,
    paths=num_paths,
)

scaled_mc = mc_value(
    OptionTypes.DIGITAL_CALL,
    spot_fx_rate,
    strike_fx_rate,
    ccy2,
    option_notional=1000.0,
    paths=num_paths,
)

tests.append((
    "Monte Carlo value scales linearly with notional",
    abs(
        scaled_mc
        - 1000.0 * unit_mc
    ) < 1.0e-5,
))

# MC accuracy should not be tested using machine precision.
#
# A digital payoff is discontinuous, so Monte Carlo convergence is relatively
# noisy. This tolerance is deliberately much looser than the analytic
# identities above.

mc_abs_tolerance = 0.01

tests.append((
    "Baseline MC is reasonably close to analytic",
    abs(mc_call - analytic_call) < mc_abs_tolerance,
))

print(
    f"{'TEST':<68s}"
    f"{'RESULT':>8s}"
)

print("-" * 76)

num_passed = 0

for name, condition in tests:
    if print_test(name, condition):
        num_passed += 1


# ============================================================================
# 13. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("13. SUMMARY")
print("=" * 78)

print(
    f"Baseline analytic value      : "
    f"{analytic_call:.8f}"
)

print(
    f"Baseline Monte Carlo value   : "
    f"{mc_call:.8f}"
)

print(
    f"Baseline MC error            : "
    f"{mc_call - analytic_call:.8f}"
)

print(
    f"Maximum spot MC error        : "
    f"{np.max(np.abs(spot_mc_errors)):.8f}"
)

print(
    f"Maximum domestic parity err  : "
    f"{np.max(np.abs(domestic_errors)):.10e}"
)

print(
    f"Maximum foreign parity err   : "
    f"{np.max(np.abs(foreign_errors)):.10e}"
)

print(
    f"Sanity tests passed          : "
    f"{num_passed}/{len(tests)}"
)
