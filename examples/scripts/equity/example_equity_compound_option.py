# ============================================================================
# FINANCEPY EXAMPLES - EquityCompoundOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

# Allow this example to run directly from its category folder.

import numpy as np
import matplotlib.pyplot as plt

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style


set_plot_style()


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(1, 1, 2015)

# The compound option expires first.
expiry_dt1 = Date(1, 1, 2017)

# The underlying option expires later.
expiry_dt2 = Date(1, 1, 2018)

# Strike of the compound option.
k1 = 5.0

# Strike of the underlying option.
k2 = 95.0

stock_price = 85.0
volatility = 0.15
interest_rate = 0.035
dividend_yield = 0.01

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)


# ============================================================================
# 1. EUROPEAN COMPOUND OPTION VALUES
# ============================================================================
# What this section demonstrates:
# Values the four possible European compound-option combinations:
#
#     Call on Call
#     Call on Put
#     Put on Call
#     Put on Put
#
# The first option type describes the compound option. The second describes
# the underlying option.
#
# At expiry_dt1 the holder decides whether the underlying option is worth
# acquiring or delivering according to the compound-option payoff.

print("\n" + "=" * 78)
print("1. EUROPEAN COMPOUND OPTION VALUES")
print("=" * 78)

european_types = [
    (
        "CALL ON CALL",
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_CALL,
    ),
    (
        "CALL ON PUT",
        OptionTypes.EUROPEAN_CALL,
        OptionTypes.EUROPEAN_PUT,
    ),
    (
        "PUT ON CALL",
        OptionTypes.EUROPEAN_PUT,
        OptionTypes.EUROPEAN_CALL,
    ),
    (
        "PUT ON PUT",
        OptionTypes.EUROPEAN_PUT,
        OptionTypes.EUROPEAN_PUT,
    ),
]

print(
    f"{'TYPE':<20}"
    f"{'K1':>10}"
    f"{'K2':>10}"
    f"{'S':>10}"
    f"{'VALUE':>16}"
)

print("-" * 66)

for label, opt_type1, opt_type2 in european_types:

    cmpd_option = EquityCompoundOption(
        expiry_dt1,
        opt_type1,
        k1,
        expiry_dt2,
        opt_type2,
        k2,
    )

    value = cmpd_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    print(
        f"{label:<20}"
        f"{k1:10.2f}"
        f"{k2:10.2f}"
        f"{stock_price:10.2f}"
        f"{value:16.8f}"
    )


# ============================================================================
# 2. EUROPEAN TREE CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Compares the numerical binomial-tree value with the analytical European
# compound-option value.
#
# A call-on-call is used as the representative convergence test. As the
# number of tree steps increases, the numerical value should approach the
# analytical value.
#
# Plotting the pricing error makes convergence much easier to see than
# plotting the two absolute prices on the same graph.

print("\n" + "=" * 78)
print("2. EUROPEAN TREE CONVERGENCE")
print("=" * 78)

cmpd_option = EquityCompoundOption(
    expiry_dt1,
    OptionTypes.EUROPEAN_CALL,
    k1,
    expiry_dt2,
    OptionTypes.EUROPEAN_CALL,
    k2,
)

analytic_value = cmpd_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

num_steps_list = [
    100,
    200,
    300,
    400,
    500,
    600,
    800,
    1000,
    1200,
    1500,
    2000,
    2500,
    3000,
]

tree_values = []
tree_errors = []

print(
    f"{'STEPS':>10}"
    f"{'ANALYTIC':>16}"
    f"{'TREE':>16}"
    f"{'ERROR':>16}"
)

print("-" * 58)

for num_steps in num_steps_list:

    tree_result = cmpd_option.value_tree(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        int(num_steps),
    )

    tree_value = tree_result[0]
    tree_error = tree_value - analytic_value

    tree_values.append(tree_value)
    tree_errors.append(tree_error)

    print(
        f"{num_steps:10d}"
        f"{analytic_value:16.8f}"
        f"{tree_value:16.8f}"
        f"{tree_error:16.8f}"
    )

tree_values = np.asarray(tree_values)
tree_errors = np.asarray(tree_errors)

plt.figure()

plt.plot(
    num_steps_list,
    tree_errors,
    marker="o",
    label="Tree Error",
)

plt.axhline(
    0.0,
    linestyle="--",
    label="Zero Error",
)

plt.xlabel("Number of Tree Steps")
plt.ylabel("Pricing Error (Tree - Analytic)")
plt.title("European Compound Option Tree Convergence Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. AMERICAN COMPOUND OPTIONS
# ============================================================================
# What this section demonstrates:
# Values American compound options using the tree.
#
# Unlike the European case, early exercise means that the analytical
# European compound-option formula is no longer the appropriate benchmark.
# We therefore examine numerical convergence as the tree is refined.

print("\n" + "=" * 78)
print("3. AMERICAN COMPOUND OPTIONS")
print("=" * 78)

american_types = [
    (
        "CALL ON CALL",
        OptionTypes.AMERICAN_CALL,
        OptionTypes.AMERICAN_CALL,
    ),
    (
        "CALL ON PUT",
        OptionTypes.AMERICAN_CALL,
        OptionTypes.AMERICAN_PUT,
    ),
    (
        "PUT ON CALL",
        OptionTypes.AMERICAN_PUT,
        OptionTypes.AMERICAN_CALL,
    ),
    (
        "PUT ON PUT",
        OptionTypes.AMERICAN_PUT,
        OptionTypes.AMERICAN_PUT,
    ),
]

american_steps = [
    100,
    200,
    500,
    1000,
]

print(
    f"{'TYPE':<20}"
    f"{'STEPS':>10}"
    f"{'VALUE':>16}"
)

print("-" * 46)

for label, opt_type1, opt_type2 in american_types:

    cmpd_option = EquityCompoundOption(
        expiry_dt1,
        opt_type1,
        k1,
        expiry_dt2,
        opt_type2,
        k2,
    )

    for num_steps in american_steps:

        tree_result = cmpd_option.value_tree(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_steps,
        )

        print(
            f"{label:<20}"
            f"{num_steps:10d}"
            f"{tree_result[0]:16.8f}"
        )


# ============================================================================
# 4. COMPOUND OPTION VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how each of the four European compound-option structures responds to
# the current stock price.
#
# This is particularly useful for understanding structures such as a
# call-on-put or put-on-call, whose economic behaviour is less obvious than
# that of a vanilla option.

print("\n" + "=" * 78)
print("4. COMPOUND OPTION VALUE VERSUS STOCK PRICE")
print("=" * 78)

stock_prices = np.linspace(
    60.0,
    120.0,
    31,
)

compound_values = {}

for label, opt_type1, opt_type2 in european_types:

    cmpd_option = EquityCompoundOption(
        expiry_dt1,
        opt_type1,
        k1,
        expiry_dt2,
        opt_type2,
        k2,
    )

    values = []

    for s in stock_prices:

        value = cmpd_option.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            model,
        )

        values.append(value)

    compound_values[label] = np.asarray(values)

print(
    f"{'STOCK':>10}"
    f"{'CALL/CALL':>14}"
    f"{'CALL/PUT':>14}"
    f"{'PUT/CALL':>14}"
    f"{'PUT/PUT':>14}"
)

print("-" * 66)

for i in range(
    0,
    len(stock_prices),
    5,
):

    print(
        f"{stock_prices[i]:10.2f}"
        f"{compound_values['CALL ON CALL'][i]:14.6f}"
        f"{compound_values['CALL ON PUT'][i]:14.6f}"
        f"{compound_values['PUT ON CALL'][i]:14.6f}"
        f"{compound_values['PUT ON PUT'][i]:14.6f}"
    )

plt.figure()

for label, _, _ in european_types:

    plt.plot(
        stock_prices,
        compound_values[label],
        label=label.title(),
    )

plt.xlabel("Stock Price")
plt.ylabel("Compound Option Value")
plt.title("European Compound Option Value versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. COMPOUND OPTION GREEKS VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Examines the delta, vega and theta of a representative European
# call-on-call compound option.
#
# A compound option has exposure both to the stock and to the value of the
# underlying option, so its Greeks can behave differently from those of a
# simple vanilla option.

print("\n" + "=" * 78)
print("5. COMPOUND OPTION GREEKS VERSUS STOCK PRICE")
print("=" * 78)

cmpd_option = EquityCompoundOption(
    expiry_dt1,
    OptionTypes.EUROPEAN_CALL,
    k1,
    expiry_dt2,
    OptionTypes.EUROPEAN_CALL,
    k2,
)

stock_prices = np.linspace(
    60.0,
    120.0,
    31,
)

values = []
deltas = []
vegas = []
thetas = []

print(
    f"{'STOCK':>10}"
    f"{'VALUE':>14}"
    f"{'DELTA':>14}"
    f"{'VEGA':>14}"
    f"{'THETA':>14}"
)

print("-" * 66)

for s in stock_prices:

    value = cmpd_option.value(
        value_dt,
        s,
        discount_curve,
        dividend_curve,
        model,
    )

    delta = cmpd_option.delta(
        value_dt,
        s,
        discount_curve,
        dividend_curve,
        model,
    )

    vega = cmpd_option.vega(
        value_dt,
        s,
        discount_curve,
        dividend_curve,
        model,
    )

    theta = cmpd_option.theta(
        value_dt,
        s,
        discount_curve,
        dividend_curve,
        model,
    )

    values.append(value)
    deltas.append(delta)
    vegas.append(vega)
    thetas.append(theta)

    print(
        f"{s:10.2f}"
        f"{value:14.6f}"
        f"{delta:14.6f}"
        f"{vega:14.6f}"
        f"{theta:14.6f}"
    )

values = np.asarray(values)
deltas = np.asarray(deltas)
vegas = np.asarray(vegas)
thetas = np.asarray(thetas)

plt.figure()

plt.plot(
    stock_prices,
    deltas,
    label="Delta",
)

plt.xlabel("Stock Price")
plt.ylabel("Delta")
plt.title("Call-on-Call Delta versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


plt.figure()

plt.plot(
    stock_prices,
    vegas,
    label="Vega",
)

plt.xlabel("Stock Price")
plt.ylabel("Vega")
plt.title("Call-on-Call Vega versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


plt.figure()

plt.plot(
    stock_prices,
    thetas,
    label="Theta",
)

plt.xlabel("Stock Price")
plt.ylabel("Theta")
plt.title("Call-on-Call Theta versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. COMPOUND OPTION VALUE VERSUS VOLATILITY
# ============================================================================
# What this section demonstrates:
# Shows how the value of the four European compound-option structures changes
# with volatility.
#
# A compound option contains optionality on an underlying option, making
# volatility an especially important input.

print("\n" + "=" * 78)
print("6. COMPOUND OPTION VALUE VERSUS VOLATILITY")
print("=" * 78)

stock_price = 85.0

volatilities = np.linspace(
    0.05,
    0.40,
    15,
)

volatility_values = {
    label: []
    for label, _, _ in european_types
}

print(
    f"{'VOL':>10}"
    f"{'CALL/CALL':>14}"
    f"{'CALL/PUT':>14}"
    f"{'PUT/CALL':>14}"
    f"{'PUT/PUT':>14}"
)

print("-" * 66)

for volatility in volatilities:

    test_model = BlackScholes(
        volatility,
    )

    row_values = {}

    for label, opt_type1, opt_type2 in european_types:

        cmpd_option = EquityCompoundOption(
            expiry_dt1,
            opt_type1,
            k1,
            expiry_dt2,
            opt_type2,
            k2,
        )

        value = cmpd_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            test_model,
        )

        volatility_values[label].append(value)
        row_values[label] = value

    print(
        f"{volatility:10.4f}"
        f"{row_values['CALL ON CALL']:14.6f}"
        f"{row_values['CALL ON PUT']:14.6f}"
        f"{row_values['PUT ON CALL']:14.6f}"
        f"{row_values['PUT ON PUT']:14.6f}"
    )

plt.figure()

for label, _, _ in european_types:

    plt.plot(
        volatilities,
        volatility_values[label],
        label=label.title(),
    )

plt.xlabel("Volatility")
plt.ylabel("Compound Option Value")
plt.title("European Compound Option Value versus Volatility")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. VALUE VERSUS COMPOUND OPTION STRIKE
# ============================================================================
# What this section demonstrates:
# Varies K1, the strike of the compound option.
#
# For a call-on-call, K1 is the amount paid at the first expiry date to
# acquire the underlying call option. Increasing this cost should reduce the
# value of the compound call.

print("\n" + "=" * 78)
print("7. VALUE VERSUS COMPOUND OPTION STRIKE")
print("=" * 78)

volatility = 0.15

model = BlackScholes(
    volatility,
)

k1_values = np.linspace(
    1.0,
    15.0,
    29,
)

call_on_call_values = []
put_on_call_values = []

print(
    f"{'K1':>10}"
    f"{'CALL ON CALL':>18}"
    f"{'PUT ON CALL':>18}"
)

print("-" * 46)

for test_k1 in k1_values:

    call_on_call = EquityCompoundOption(
        expiry_dt1,
        OptionTypes.EUROPEAN_CALL,
        test_k1,
        expiry_dt2,
        OptionTypes.EUROPEAN_CALL,
        k2,
    )

    put_on_call = EquityCompoundOption(
        expiry_dt1,
        OptionTypes.EUROPEAN_PUT,
        test_k1,
        expiry_dt2,
        OptionTypes.EUROPEAN_CALL,
        k2,
    )

    call_value = call_on_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    put_value = put_on_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    call_on_call_values.append(call_value)
    put_on_call_values.append(put_value)

    print(
        f"{test_k1:10.4f}"
        f"{call_value:18.8f}"
        f"{put_value:18.8f}"
    )

plt.figure()

plt.plot(
    k1_values,
    call_on_call_values,
    label="Call on Call",
)

plt.plot(
    k1_values,
    put_on_call_values,
    label="Put on Call",
)

plt.xlabel("Compound Option Strike K1")
plt.ylabel("Compound Option Value")
plt.title("Compound Option Value versus Compound Strike")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. VALUE VERSUS UNDERLYING OPTION STRIKE
# ============================================================================
# What this section demonstrates:
# Varies K2, the strike of the underlying option.
#
# This experiment separates the strike of the compound option from the strike
# of the underlying vanilla option and shows how changing the underlying
# option economics feeds through into the compound option value.

print("\n" + "=" * 78)
print("8. VALUE VERSUS UNDERLYING OPTION STRIKE")
print("=" * 78)

k2_values = np.linspace(
    70.0,
    120.0,
    26,
)

call_on_call_values = []
call_on_put_values = []

print(
    f"{'K2':>10}"
    f"{'CALL ON CALL':>18}"
    f"{'CALL ON PUT':>18}"
)

print("-" * 46)

for test_k2 in k2_values:

    call_on_call = EquityCompoundOption(
        expiry_dt1,
        OptionTypes.EUROPEAN_CALL,
        k1,
        expiry_dt2,
        OptionTypes.EUROPEAN_CALL,
        test_k2,
    )

    call_on_put = EquityCompoundOption(
        expiry_dt1,
        OptionTypes.EUROPEAN_CALL,
        k1,
        expiry_dt2,
        OptionTypes.EUROPEAN_PUT,
        test_k2,
    )

    call_call_value = call_on_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    call_put_value = call_on_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    call_on_call_values.append(
        call_call_value,
    )

    call_on_put_values.append(
        call_put_value,
    )

    print(
        f"{test_k2:10.4f}"
        f"{call_call_value:18.8f}"
        f"{call_put_value:18.8f}"
    )

plt.figure()

plt.plot(
    k2_values,
    call_on_call_values,
    label="Call on Call",
)

plt.plot(
    k2_values,
    call_on_put_values,
    label="Call on Put",
)

plt.xlabel("Underlying Option Strike K2")
plt.ylabel("Compound Option Value")
plt.title("Compound Option Value versus Underlying Option Strike")
plt.grid(True)
plt.legend()
plt.show()
