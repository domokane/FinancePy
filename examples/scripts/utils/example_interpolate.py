# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import math
import time
import numpy as np
import matplotlib.pyplot as plt


from financepy.market.curves.interpolator import Interpolator, InterpTypes

# ============================================================================
# FINANCEPY EXAMPLES - Interpolate
# ============================================================================



PLOT_GRAPHS = False


x_values = np.array([0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 10.0])
a = -0.1
b = 0.002

y_values = []
for x in x_values:
    y = math.exp(a * x + b * x * x)
    y_values.append(y)

y_values = np.array(y_values)

x_interpolate_values = np.linspace(0.0, 10.0, 20)


########################################################################################




########################################################################################




########################################################################################




########################################################################################



########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. FIN INTERPOLATE
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN INTERPOLATE")
print("=" * 78)

print("METHOD", "X", "Y_INTERPOLATED")

for interp_type in InterpTypes:

    y_interp_values = []
    start = time.time()

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    for x in x_interpolate_values:
        y_int = interpolator.interpolate(x)
        print(interp_type, x, y_int)
        y_interp_values.append(y_int)

    end = time.time()

    if PLOT_GRAPHS:
        plt.figure(figsize=(12, 10))
        plt.plot(x_values, y_values, color="r", marker="o")
        plt.plot(
            x_interpolate_values,
            y_interp_values,
            color="b",
            label=str(interp_type),
        )
        plt.legend()

xp = np.array([0.2, 0.4, 0.45, 0.6, 0.82, 0.93, 0.99])
yp = np.array([0.4, 0.9, 0.32, 0.2, 0.22, 0.10, 0.28])
n = 10000

print("LABEL", "TIME")
interpolator = Interpolator(interp_type)
interpolator.fit(xp, yp)

start = time.time()
for _ in range(0, n):
    interpolator.interpolate(0.8)
end = time.time()
print("10000 Interpolations", end - start)

# ============================================================================
# 2. FIN INTERPOLATE  RUNS
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. FIN INTERPOLATE  RUNS")
print("=" * 78)

for interp_type in InterpTypes:

    y_interp_values = []

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    for x in x_interpolate_values:
        y_int = interpolator.interpolate(x)
        y_interp_values.append(y_int)

# ============================================================================
# 3. FIN INTERPOLATE  RECOVERS  INPUTS
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("3. FIN INTERPOLATE  RECOVERS  INPUTS")
print("=" * 78)

for interp_type in InterpTypes:

    y_interp_values = []

    interpolator = Interpolator(interp_type)
    interpolator.fit(x_values, y_values)

    for x in x_values:
        y_int = interpolator.interpolate(x)
        y_interp_values.append(y_int)

    y_interp_values = np.array(y_interp_values)
    assert (
        np.linalg.norm(y_values - y_interp_values)
        / np.linalg.norm(y_values)
        <= 1e-6
    )

# ============================================================================
# 4. FLAT FWD RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("4. FLAT FWD RATES")
print("=" * 78)

interp_type = InterpTypes.FLAT_FWD_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 0
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)

assert round(x, 4) == 0.0
assert round(y_int, 4) == 1.0

index = 5
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 2.6316
assert round(y_int, 4) == 0.7797

index = 10
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 5.2632
assert round(y_int, 4) == 0.6260

# ============================================================================
# 5. LINEAR ZERO RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("5. LINEAR ZERO RATES")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ZERO_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 8
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 4.2105
assert round(y_int, 4) == 0.6800

index = 13
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 6.8421
assert round(y_int, 4) == 0.5540

index = 18
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 9.4737
assert round(y_int, 4) == 0.4640

# ============================================================================
# 6. FINCUBIC ZERO RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("6. FINCUBIC ZERO RATES")
print("=" * 78)

interp_type = InterpTypes.FINCUBIC_ZERO_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 1
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 0.5263
assert round(y_int, 4) == 0.9493

index = 6
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 3.1579
assert round(y_int, 4) == 0.7439

index = 11
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 5.7895
assert round(y_int, 4) == 0.6007

# ============================================================================
# 7. NATCUBIC LOG DISCOUNT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("7. NATCUBIC LOG DISCOUNT")
print("=" * 78)

interp_type = InterpTypes.NATCUBIC_LOG_DISCOUNT

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 4
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 2.1053
assert round(y_int, 4) == 0.8174

index = 9
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 4.7368
assert round(y_int, 4) == 0.6512

index = 14
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 7.3684
assert round(y_int, 4) == 0.5355

# ============================================================================
# 8. NATCUBIC ZERO RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("8. NATCUBIC ZERO RATES")
print("=" * 78)

interp_type = InterpTypes.NATCUBIC_ZERO_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 2
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 1.0526
assert round(y_int, 4) == 0.9021

index = 7
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 3.6842
assert round(y_int, 4) == 0.7109

index = 12
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 6.3158
assert round(y_int, 4) == 0.5759

# ============================================================================
# 9. PCHIP ZERO RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("9. PCHIP ZERO RATES")
print("=" * 78)

interp_type = InterpTypes.PCHIP_ZERO_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 0
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 0.0
assert round(y_int, 4) == 1.0

index = 5
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 2.6316
assert round(y_int, 4) == 0.7793

index = 10
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 5.2632
assert round(y_int, 4) == 0.6244

# ============================================================================
# 10. PCHIP LOG DISCOUNT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("10. PCHIP LOG DISCOUNT")
print("=" * 78)

interp_type = InterpTypes.PCHIP_LOG_DISCOUNT

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 3
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 1.5789
assert round(y_int, 4) == 0.8582

index = 8
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 4.2105
assert round(y_int, 4) == 0.6796

index = 13
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(x, 4) == 6.8421
assert round(y_int, 4) == 0.5551

# ============================================================================
# 11. LINEAR ONFWD RATES EMPTY FIT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("11. LINEAR ONFWD RATES EMPTY FIT")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ONFWD_RATES

interpolator = Interpolator(interp_type)

interpolator.fit([], [])
assert round(interpolator.interpolate(1.0), 4) == 1.0

# ============================================================================
# 12. LINEAR ONFWD RATES SINGLE VALUE AT ORIGIN
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("12. LINEAR ONFWD RATES SINGLE VALUE AT ORIGIN")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ONFWD_RATES

interpolator = Interpolator(interp_type)
interpolator.fit([0.0], [1.0])
assert round(interpolator.interpolate(1.0), 4) == 1.0

# ============================================================================
# 13. LINEAR ONFWD RATES SINGLE VALUE NOT AT ORIGIN
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("13. LINEAR ONFWD RATES SINGLE VALUE NOT AT ORIGIN")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ONFWD_RATES

interpolator = Interpolator(interp_type)
interpolator.fit([0.1], [0.9])
assert round(interpolator.interpolate(0.0), 4) == 1.0
assert round(interpolator.interpolate(0.05), 4) == 0.9487
assert round(interpolator.interpolate(0.1), 4) == 0.9
assert round(interpolator.interpolate(1.0), 4) == 0.3487

# ============================================================================
# 14. LINEAR ONFWD RATES TWO VALUES INCLUDING ORIGIN
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("14. LINEAR ONFWD RATES TWO VALUES INCLUDING ORIGIN")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ONFWD_RATES

interpolator = Interpolator(interp_type)
interpolator.fit([0.0, 0.1], [1.0, 0.9])
assert round(interpolator.interpolate(0.0), 4) == 1.0
assert round(interpolator.interpolate(0.05), 4) == 0.9487
assert round(interpolator.interpolate(0.1), 4) == 0.9
assert round(interpolator.interpolate(1.0), 4) == 0.3487

# ============================================================================
# 15. LINEAR ONFWD RATES
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("15. LINEAR ONFWD RATES")
print("=" * 78)

interp_type = InterpTypes.LINEAR_ONFWD_RATES

interpolator = Interpolator(interp_type)
interpolator.fit(x_values, y_values)

index = 3
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(y_int, 4) == 0.8583

index = 8
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(y_int, 4) == 0.6802

index = 13
x = x_interpolate_values[index]
y_int = interpolator.interpolate(x)
assert round(y_int, 4) == 0.5537

