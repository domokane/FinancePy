# ============================================================================
# FINANCEPY EXAMPLES - EquityVanillaOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example is deliberately more verbose than a unit test.  It is intended
# to be read, modified, and run interactively by somebody learning how
# EquityVanillaOption works in FinancePy.
#
# The script demonstrates:
#
#   1. Black-Scholes valuation of European calls and puts
#   2. Monte Carlo valuation and convergence with the number of paths
#   3. Pseudo-random versus Sobol Monte Carlo sampling
#   4. Option value as a function of the stock price
#   5. Intrinsic value and time value
#   6. Call and put Greeks as functions of the stock price
#   7. Strike sensitivity
#   8. Time-to-expiry sensitivity
#   9. Implied-volatility recovery
#  10. Several plots that make the numerical results easier to interpret
#
# Notes
# -----
# * The discount curve represents the risk-free interest rate.
# * The dividend curve represents the continuously compounded dividend yield.
# * BlackScholes(volatility) supplies the constant volatility used for pricing.
# * Monte Carlo estimates contain sampling error, so they need not equal the
#   analytic Black-Scholes value exactly.
# * Sobol sampling is a low-discrepancy alternative to pseudo-random sampling.
#
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption

# ============================================================================
# 0. COMMON MARKET DATA
# ============================================================================
#
# We use one common market setup for most examples so that results from
# different sections can be compared directly.
#
# Valuation date : 1 January 2015
# Expiry date    : 1 July 2015
# Spot price     : 100
# Strike         : 100
# Volatility     : 30%
# Risk-free rate : 5%
# Dividend yield : 1%
#
# The option is therefore initially at-the-money.

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 7, 2015)

stock_price = 100.0
strike_price = 100.0
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01

# FinancePy represents rates using discount curves.  A FlatDiscountCurve is
# convenient for examples because the continuously compounded rate is constant.
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

# The Black-Scholes model object contains the volatility assumption.
model = BlackScholes(volatility)

# Construct one European call and one European put with the same strike/expiry.
call_option = EquityVanillaOption(
    expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
)
put_option = EquityVanillaOption(
    expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT
)


def section(title):
    """Print a consistent heading for each example section."""
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ============================================================================
# 1. BASIC BLACK-SCHOLES VALUATION
# ============================================================================

section("1. BASIC BLACK-SCHOLES VALUATION")

call_value = call_option.value(
    value_dt, stock_price, discount_curve, dividend_curve, model
)
put_value = put_option.value(
    value_dt, stock_price, discount_curve, dividend_curve, model
)

print(f"Stock price          : {stock_price:12.4f}")
print(f"Strike price         : {strike_price:12.4f}")
print(f"Volatility           : {volatility:12.4%}")
print(f"Interest rate        : {interest_rate:12.4%}")
print(f"Dividend yield       : {dividend_yield:12.4%}")
print(f"European call value  : {call_value:12.8f}")
print(f"European put value   : {put_value:12.8f}")


# ============================================================================
# 2. INTRINSIC VALUE AND TIME VALUE
# ============================================================================
#
# An option value can be thought of as:
#
#       option value = intrinsic value + time value
#
# Before expiry, time value is generally positive because there is still a
# possibility that future stock-price movements improve the payoff.

section("2. INTRINSIC VALUE AND TIME VALUE")

call_intrinsic = call_option.intrinsic(
    value_dt, stock_price, discount_curve, dividend_curve
)
put_intrinsic = put_option.intrinsic(
    value_dt, stock_price, discount_curve, dividend_curve
)

print(f"{'OPTION':>12s}{'VALUE':>16s}{'INTRINSIC':>16s}{'TIME VALUE':>16s}")
print("-" * 60)
print(
    f"{'CALL':>12s}{call_value:16.8f}{call_intrinsic:16.8f}"
    f"{call_value - call_intrinsic:16.8f}"
)
print(
    f"{'PUT':>12s}{put_value:16.8f}{put_intrinsic:16.8f}"
    f"{put_value - put_intrinsic:16.8f}"
)


# ============================================================================
# 3. MONTE CARLO PATH CONVERGENCE
# ============================================================================
#
# The analytic Black-Scholes value is our benchmark.  We now increase the
# number of Monte Carlo paths and observe how the estimate behaves.
#
# Because Monte Carlo is statistical, convergence is not necessarily monotonic:
# doubling the number of paths does not guarantee that every single estimate
# will be closer to the analytic value.

section("3. MONTE CARLO PATH CONVERGENCE")

num_paths_list = [1000, 2000, 5000, 10000, 20000, 50000, 100000]

mc_call_values = []
mc_call_errors = []
mc_times = []

print(
    f"{'PATHS':>12s}{'BS VALUE':>16s}{'MC VALUE':>16s}"
    f"{'MC - BS':>16s}{'TIME (S)':>14s}"
)
print("-" * 74)

for num_paths in num_paths_list:
    start = time.time()

    value_mc = call_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
    )

    elapsed = time.time() - start
    error = value_mc - call_value

    mc_call_values.append(value_mc)
    mc_call_errors.append(error)
    mc_times.append(elapsed)

    print(
        f"{num_paths:12d}{call_value:16.8f}{value_mc:16.8f}"
        f"{error:16.8f}{elapsed:14.6f}"
    )

# Plot the Monte Carlo estimate against the analytic benchmark.
plt.figure(figsize=(10, 6))
plt.semilogx(
    num_paths_list,
    mc_call_values,
    marker="o",
    label="Monte Carlo",
)
plt.axhline(call_value, linestyle="--", label="Black-Scholes")
plt.xlabel("Number of Monte Carlo paths")
plt.ylabel("European call value")
plt.title("European Call: Monte Carlo Path Convergence")
plt.grid(True)
plt.legend()
plt.show()

# Plot the pricing error directly.  The zero line corresponds to an exact match
# with the analytic Black-Scholes value.
plt.figure(figsize=(10, 6))
plt.semilogx(
    num_paths_list,
    mc_call_errors,
    marker="o",
    label="MC - Black-Scholes",
)
plt.axhline(0.0, linestyle="--", label="Zero error")
plt.xlabel("Number of Monte Carlo paths")
plt.ylabel("Pricing error")
plt.title("European Call: Monte Carlo Pricing Error")
plt.grid(True)
plt.legend()
plt.show()

# Runtime is also useful to inspect: more paths normally improve statistical
# precision, but they require more computation.
plt.figure(figsize=(10, 6))
plt.semilogx(num_paths_list, mc_times, marker="o")
plt.xlabel("Number of Monte Carlo paths")
plt.ylabel("Runtime (seconds)")
plt.title("European Call: Monte Carlo Runtime")
plt.grid(True)
plt.show()


# ============================================================================
# 4. PSEUDO-RANDOM VERSUS SOBOL MONTE CARLO
# ============================================================================
#
# FinancePy's value_mc method can use ordinary pseudo-random numbers or Sobol
# low-discrepancy numbers.  Here we compare both methods over a range of spots.
#
# This is a useful practical example because the analytic Black-Scholes value
# lets us see the error produced by each simulation method.

section("4. PSEUDO-RANDOM VERSUS SOBOL MONTE CARLO")

comparison_spots = np.arange(80.0, 121.0, 10.0)
num_paths = 100000

call_bs_values = []
call_mc_random = []
call_mc_sobol = []
call_random_errors = []
call_sobol_errors = []

print(
    f"{'SPOT':>10s}{'BS':>14s}{'MC RANDOM':>14s}{'MC SOBOL':>14s}"
    f"{'ERR RANDOM':>14s}{'ERR SOBOL':>14s}"
)
print("-" * 80)

for spot in comparison_spots:
    value_bs = call_option.value(
        value_dt, spot, discount_curve, dividend_curve, model
    )

    value_random = call_option.value_mc(
        value_dt,
        spot,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
        False,
    )

    value_sobol = call_option.value_mc(
        value_dt,
        spot,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
        True,
    )

    err_random = value_random - value_bs
    err_sobol = value_sobol - value_bs

    call_bs_values.append(value_bs)
    call_mc_random.append(value_random)
    call_mc_sobol.append(value_sobol)
    call_random_errors.append(err_random)
    call_sobol_errors.append(err_sobol)

    print(
        f"{spot:10.2f}{value_bs:14.8f}{value_random:14.8f}"
        f"{value_sobol:14.8f}{err_random:14.8f}{err_sobol:14.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(comparison_spots, call_bs_values, marker="o", label="Black-Scholes")
plt.plot(comparison_spots, call_mc_random, marker="o", label="MC pseudo-random")
plt.plot(comparison_spots, call_mc_sobol, marker="o", label="MC Sobol")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("European call value")
plt.title("European Call: Analytic vs Monte Carlo")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(
    comparison_spots,
    call_random_errors,
    marker="o",
    label="Pseudo-random error",
)
plt.plot(
    comparison_spots,
    call_sobol_errors,
    marker="o",
    label="Sobol error",
)
plt.axhline(0.0, linestyle="--", label="Zero error")
plt.xlabel("Stock price")
plt.ylabel("MC value - Black-Scholes value")
plt.title("European Call: Monte Carlo Error by Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. CALL AND PUT VALUES VERSUS STOCK PRICE
# ============================================================================
#
# A call becomes more valuable as the stock price rises.
# A put becomes more valuable as the stock price falls.
#
# The curves are nonlinear because option payoffs are asymmetric.

section("5. CALL AND PUT VALUES VERSUS STOCK PRICE")

spot_grid = np.linspace(60.0, 140.0, 81)

call_values = np.array(
    [
        call_option.value(
            value_dt, spot, discount_curve, dividend_curve, model
        )
        for spot in spot_grid
    ]
)

put_values = np.array(
    [
        put_option.value(
            value_dt, spot, discount_curve, dividend_curve, model
        )
        for spot in spot_grid
    ]
)

plt.figure(figsize=(10, 6))
plt.plot(spot_grid, call_values, label="European call")
plt.plot(spot_grid, put_values, label="European put")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Option value")
plt.title("European Option Value versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. PAYOFF AT EXPIRY
# ============================================================================
#
# At expiry there is no remaining time value:
#
#       Call payoff = max(S_T - K, 0)
#       Put payoff  = max(K - S_T, 0)
#
# This plot is not a FinancePy valuation call; it is a direct illustration of
# the contract payoff and is useful for understanding the option's shape.

section("6. PAYOFF AT EXPIRY")

terminal_spots = np.linspace(60.0, 140.0, 161)
call_payoff = np.maximum(terminal_spots - strike_price, 0.0)
put_payoff = np.maximum(strike_price - terminal_spots, 0.0)

plt.figure(figsize=(10, 6))
plt.plot(terminal_spots, call_payoff, label="Call payoff")
plt.plot(terminal_spots, put_payoff, label="Put payoff")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price at expiry")
plt.ylabel("Payoff")
plt.title("European Call and Put Payoffs at Expiry")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. GREEKS VERSUS STOCK PRICE
# ============================================================================
#
# Greeks measure how the option value responds to changes in market inputs.
#
# Delta : sensitivity to the stock price
# Vega  : sensitivity to volatility
# Theta : sensitivity to the passage of time
# Rho   : sensitivity to the interest rate
# Vanna : cross-sensitivity involving spot and volatility
#
# FinancePy provides these directly on EquityVanillaOption.

section("7. GREEKS VERSUS STOCK PRICE")

greek_spots = np.linspace(70.0, 130.0, 61)

call_delta = []
call_vega = []
call_theta = []
call_rho = []
call_vanna = []

put_delta = []
put_vega = []
put_theta = []
put_rho = []
put_vanna = []

for spot in greek_spots:
    call_delta.append(
        call_option.delta(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    call_vega.append(
        call_option.vega(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    call_theta.append(
        call_option.theta(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    call_rho.append(
        call_option.rho(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    call_vanna.append(
        call_option.vanna(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )

    put_delta.append(
        put_option.delta(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    put_vega.append(
        put_option.vega(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    put_theta.append(
        put_option.theta(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    put_rho.append(
        put_option.rho(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )
    put_vanna.append(
        put_option.vanna(
            value_dt, spot, discount_curve, dividend_curve, model
        )
    )

# Delta
plt.figure(figsize=(10, 6))
plt.plot(greek_spots, call_delta, label="Call delta")
plt.plot(greek_spots, put_delta, label="Put delta")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Delta")
plt.title("European Option Delta versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()

# Vega
plt.figure(figsize=(10, 6))
plt.plot(greek_spots, call_vega, label="Call vega")
plt.plot(greek_spots, put_vega, label="Put vega")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Vega")
plt.title("European Option Vega versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()

# Theta
plt.figure(figsize=(10, 6))
plt.plot(greek_spots, call_theta, label="Call theta")
plt.plot(greek_spots, put_theta, label="Put theta")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Theta")
plt.title("European Option Theta versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()

# Rho
plt.figure(figsize=(10, 6))
plt.plot(greek_spots, call_rho, label="Call rho")
plt.plot(greek_spots, put_rho, label="Put rho")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Rho")
plt.title("European Option Rho versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()

# Vanna
plt.figure(figsize=(10, 6))
plt.plot(greek_spots, call_vanna, label="Call vanna")
plt.plot(greek_spots, put_vanna, label="Put vanna")
plt.axvline(strike_price, linestyle="--", label="Strike")
plt.xlabel("Stock price")
plt.ylabel("Vanna")
plt.title("European Option Vanna versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. OPTION VALUE VERSUS STRIKE
# ============================================================================
#
# Holding spot fixed:
# * a higher strike generally reduces a call's value;
# * a higher strike generally increases a put's value.
#
# We create a new option for every strike because strike is a contract term.

section("8. OPTION VALUE VERSUS STRIKE")

strike_grid = np.linspace(60.0, 140.0, 81)
call_values_by_strike = []
put_values_by_strike = []

for strike in strike_grid:
    call = EquityVanillaOption(
        expiry_dt, strike, OptionTypes.EUROPEAN_CALL
    )
    put = EquityVanillaOption(
        expiry_dt, strike, OptionTypes.EUROPEAN_PUT
    )

    call_values_by_strike.append(
        call.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )
    put_values_by_strike.append(
        put.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

plt.figure(figsize=(10, 6))
plt.plot(strike_grid, call_values_by_strike, label="European call")
plt.plot(strike_grid, put_values_by_strike, label="European put")
plt.axvline(stock_price, linestyle="--", label="Current spot")
plt.xlabel("Strike price")
plt.ylabel("Option value")
plt.title("European Option Value versus Strike")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. OPTION VALUE VERSUS TIME TO EXPIRY
# ============================================================================
#
# Here the market inputs are held fixed while expiry is moved farther into the
# future.  This is an intuitive way to see the effect of remaining maturity.

section("9. OPTION VALUE VERSUS TIME TO EXPIRY")

times_to_expiry = np.array(
    [0.01, 0.03, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00]
)

call_values_by_time = []
put_values_by_time = []

for time_to_expiry in times_to_expiry:
    maturity = value_dt.add_years(float(time_to_expiry))

    call = EquityVanillaOption(
        maturity, strike_price, OptionTypes.EUROPEAN_CALL
    )
    put = EquityVanillaOption(
        maturity, strike_price, OptionTypes.EUROPEAN_PUT
    )

    call_values_by_time.append(
        call.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )
    put_values_by_time.append(
        put.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

plt.figure(figsize=(10, 6))
plt.plot(times_to_expiry, call_values_by_time, marker="o", label="Call")
plt.plot(times_to_expiry, put_values_by_time, marker="o", label="Put")
plt.xlabel("Time to expiry (years)")
plt.ylabel("Option value")
plt.title("European Option Value versus Time to Expiry")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 10. IMPLIED VOLATILITY - SINGLE EXAMPLE
# ============================================================================
#
# Implied volatility reverses the pricing problem:
#
#   forward problem:  volatility -> option value
#   inverse problem:  option value -> implied volatility
#
# We first price an option using a known 30% volatility, then ask FinancePy to
# recover the volatility from that option value.

section("10. IMPLIED VOLATILITY - SINGLE EXAMPLE")

market_value = call_option.value(
    value_dt, stock_price, discount_curve, dividend_curve, model
)

implied_vol = call_option.implied_volatility(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    market_value,
)

# The recovery errors below are extremely small (around 1e-12).
# They arise from floating-point arithmetic and the numerical root solver
# used to invert the Black-Scholes price.  The visible structure in the
# graph is therefore numerical rather than economically significant.
#
# Plotting the error on its natural scale is useful because it demonstrates
# how accurately implied_volatility() recovers the volatility used to
# generate the option price.

print(f"Input volatility     : {volatility:.10f}")
print(f"Option value         : {market_value:.10f}")
print(f"Recovered implied vol: {implied_vol:.10f}")
print(f"Recovery error       : {implied_vol - volatility:.10e}")


# ============================================================================
# 11. IMPLIED VOLATILITY RECOVERY ACROSS INPUT VOLATILITIES
# ============================================================================
#
# This example repeats the previous calculation over a range of volatilities.
# If the pricing and inversion are numerically consistent, recovered implied
# volatility should lie very close to the 45-degree line y = x.

section("11. IMPLIED VOLATILITY RECOVERY ACROSS INPUT VOLATILITIES")

input_vols = np.linspace(0.05, 0.80, 16)
recovered_vols = []
recovery_errors = []

print(f"{'INPUT VOL':>14s}{'OPTION VALUE':>18s}{'IMPLIED VOL':>18s}{'ERROR':>16s}")
print("-" * 66)

for input_vol in input_vols:
    input_model = BlackScholes(input_vol)

    option_value = call_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        input_model,
    )

    recovered_vol = call_option.implied_volatility(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        option_value,
    )

    error = recovered_vol - input_vol

    recovered_vols.append(recovered_vol)
    recovery_errors.append(error)

    print(
        f"{input_vol:14.8f}{option_value:18.8f}"
        f"{recovered_vol:18.8f}{error:16.8e}"
    )

plt.figure(figsize=(10, 6))
plt.plot(input_vols, recovered_vols, marker="o", label="Recovered implied vol")
plt.plot(input_vols, input_vols, linestyle="--", label="Perfect recovery")
plt.xlabel("Input volatility")
plt.ylabel("Recovered implied volatility")
plt.title("Implied Volatility Recovery")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(input_vols, recovery_errors, marker="o")
plt.axhline(0.0, linestyle="--")
plt.xlabel("Input volatility")
plt.ylabel("Recovered vol - input vol")
plt.title("Implied Volatility Recovery Error")
plt.grid(True)
plt.show()


# ============================================================================
# 12. IMPLIED VOLATILITY ACROSS STRIKES
# ============================================================================
#
# Under a constant-volatility Black-Scholes model, if option prices are
# generated with one constant volatility and then inverted consistently, the
# implied volatility should be approximately flat across strike.
#
# This is a useful baseline before studying real market volatility smiles.

section("12. IMPLIED VOLATILITY ACROSS STRIKES")

iv_strikes = np.linspace(70.0, 130.0, 25)
implied_vols_by_strike = []

for strike in iv_strikes:
    option = EquityVanillaOption(
        expiry_dt, strike, OptionTypes.EUROPEAN_CALL
    )

    option_value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    intrinsic = option.intrinsic(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
    )

    # Deep in-the-money or very short-dated options can have almost no time
    # value.  Implied-volatility inversion is numerically less informative in
    # that limit, so we only invert when there is meaningful time value.
    if option_value - intrinsic > 1.0e-10:
        iv = option.implied_volatility(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            option_value,
        )
    else:
        iv = np.nan

    implied_vols_by_strike.append(iv)

plt.figure(figsize=(10, 6))
plt.plot(
    iv_strikes,
    implied_vols_by_strike,
    marker="o",
    label="Recovered implied volatility",
)
plt.axhline(volatility, linestyle="--", label="Input volatility")
plt.axvline(stock_price, linestyle="--", label="Current spot")
plt.xlabel("Strike price")
plt.ylabel("Implied volatility")
plt.title("Black-Scholes Implied Volatility versus Strike")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 13. SMALL GREEKS TABLE AROUND THE STRIKE
# ============================================================================
#
# A compact table is often easier to inspect than a plot when checking exact
# numbers.  These spots show the transition from out-of-the-money through
# at-the-money to in-the-money for the call.

section("13. SMALL GREEKS TABLE AROUND THE STRIKE")

table_spots = [80.0, 90.0, 100.0, 110.0, 120.0]

print(
    f"{'SPOT':>8s}{'VALUE':>14s}{'DELTA':>14s}{'VEGA':>14s}"
    f"{'THETA':>14s}{'RHO':>14s}{'VANNA':>14s}"
)
print("-" * 92)

for spot in table_spots:
    value = call_option.value(
        value_dt, spot, discount_curve, dividend_curve, model
    )
    delta = call_option.delta(
        value_dt, spot, discount_curve, dividend_curve, model
    )
    vega = call_option.vega(
        value_dt, spot, discount_curve, dividend_curve, model
    )
    theta = call_option.theta(
        value_dt, spot, discount_curve, dividend_curve, model
    )
    rho = call_option.rho(
        value_dt, spot, discount_curve, dividend_curve, model
    )
    vanna = call_option.vanna(
        value_dt, spot, discount_curve, dividend_curve, model
    )

    print(
        f"{spot:8.2f}{value:14.8f}{delta:14.8f}{vega:14.8f}"
        f"{theta:14.8f}{rho:14.8f}{vanna:14.8f}"
    )


# ============================================================================
# 14. SUMMARY
# ============================================================================

section("14. SUMMARY")

print(f"European call value        : {call_value:.10f}")
print(f"European put value         : {put_value:.10f}")
print(f"Call intrinsic value       : {call_intrinsic:.10f}")
print(f"Put intrinsic value        : {put_intrinsic:.10f}")
print(f"Last MC call estimate      : {mc_call_values[-1]:.10f}")
print(f"Last MC call error         : {mc_call_errors[-1]:.10f}")
print(f"Recovered implied vol      : {implied_vol:.10f}")
print()
print("The figures above illustrate valuation, Monte Carlo convergence,")
print("payoffs, Greeks, strike/maturity sensitivity, and implied volatility.")
