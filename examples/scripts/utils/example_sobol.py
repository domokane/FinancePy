
# Allow this example to run directly from its category folder.
import time


from financepy.models.sobol import get_uniform_sobol, get_gaussian_sobol

# ============================================================================
# FINANCEPY EXAMPLES - Sobol
# ============================================================================


########################################################################################






########################################################################################

# ============================================================================
# 1. FIN SOBOL
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN SOBOL")
print("=" * 78)

num_points = 1000
dimensions = 3

points = get_uniform_sobol(num_points, dimensions)

for d in range(dimensions):
    av = 0.0
    var = 0.0

    for point in points[:, d]:
        av += point
        var += point**2

    av /= num_points
    var /= num_points

    av_error = abs(av - (1 / 2))
    var_error = abs(var - (1 / 3))
    assert av_error < 0.002
    assert var_error < 0.002

num_repeats = 100
num_dimensions = 10

print("LABEL", "TIME")
start = time.time()
for _ in range(num_repeats):
    get_uniform_sobol(1000, num_dimensions)
end = time.time()
print("Average time taken", (end - start) / num_repeats)

start = time.time()
for _ in range(num_repeats):
    get_gaussian_sobol(1000, num_dimensions)
end = time.time()
print("Average time taken", (end - start) / num_repeats)

# ============================================================================
# 2. FIN SOBOL CACHE
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("2. FIN SOBOL CACHE")
print("=" * 78)



