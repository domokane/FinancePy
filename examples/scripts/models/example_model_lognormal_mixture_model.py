
# Allow this example to run directly from its category folder.
import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.models.black import Black
from financepy.models.lognormal_mixture_model import LognormalMixtureModel
from financepy.models.lognormal_mixture_surface import LognormalMixtureSurface

# ============================================================================
# FINANCEPY EXAMPLES - LognormalMixtureModel
# ============================================================================




















#########################################################################################




#########################################################################################

# ============================================================================
# 1. FORWARD CONSTRAINT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. FORWARD CONSTRAINT")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

p = 0.4
displacement = -0.05

F1, F2 = model._component_forwards(p, displacement)

mixture_forward = p * F1 + (1.0 - p) * F2

assert np.isclose(mixture_forward, model.F, atol=1e-12)

# ============================================================================
# 2. MIXTURE REDUCES TO BLACK
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("2. MIXTURE REDUCES TO BLACK")
print("=" * 78)

"""
If both mixture components have the same forward and volatility,
the mixture price must equal the FinancePy Black price.
"""

F = 100.0
T = 1.0
r = 0.03
sigma = 0.20

strikes = np.array([80.0, 90.0, 100.0, 110.0, 120.0])

model = LognormalMixtureModel(F=F, T=T, r=r)

mixture_prices = model.price(K=strikes, p=0.35, displacement=0.0, sigma1=sigma, sigma2=sigma)

black = Black(sigma)

discount_factor = np.exp(-r * T)

black_prices = np.array([black.value(F, K, T, discount_factor, OptionTypes.EUROPEAN_CALL) for K in strikes])

assert np.allclose(mixture_prices, black_prices, atol=1e-10)

# ============================================================================
# 3. MIXTURE PRICE DECREASES WITH STRIKE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("3. MIXTURE PRICE DECREASES WITH STRIKE")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

strikes = np.linspace(70.0, 130.0, 50)

prices = model.price(K=strikes, p=0.4, displacement=-0.05, sigma1=0.15, sigma2=0.30)

assert np.all(np.diff(prices) <= 1e-12)

# ============================================================================
# 4. MIXTURE PRICE CONVEX IN STRIKE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("4. MIXTURE PRICE CONVEX IN STRIKE")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0, r=0.03)

strikes = np.linspace(70.0, 130.0, 100)

prices = model.price(K=strikes, p=0.4, displacement=-0.05, sigma1=0.15, sigma2=0.30)

second_diff = prices[:-2] - 2.0 * prices[1:-1] + prices[2:]

assert np.all(second_diff >= -1e-10)

# ============================================================================
# 5. DENSITY POSITIVE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("5. DENSITY POSITIVE")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0)

model.p = 0.4
model.displacement = -0.05
model.sigma1 = 0.15
model.sigma2 = 0.30

model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

s = np.linspace(1e-4, 400.0, 10000)

density = model.density(s)

assert np.all(density >= 0.0)

# ============================================================================
# 6. DENSITY INTEGRATES TO ONE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("6. DENSITY INTEGRATES TO ONE")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0)

model.p = 0.4
model.displacement = -0.05
model.sigma1 = 0.15
model.sigma2 = 0.30

model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

s = np.linspace(1e-4, 500.0, 100000)

density = model.density(s)

integral = np.trapezoid(density, s)

assert np.isclose(integral, 1.0, atol=1e-4)

# ============================================================================
# 7. DENSITY MEAN EQUALS FORWARD
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("7. DENSITY MEAN EQUALS FORWARD")
print("=" * 78)

model = LognormalMixtureModel(F=100.0, T=1.0)

model.p = 0.4
model.displacement = -0.05
model.sigma1 = 0.15
model.sigma2 = 0.30

model.F1, model.F2 = model._component_forwards(model.p, model.displacement)

s = np.linspace(1e-4, 500.0, 100000)

density = model.density(s)

expected_s = np.trapezoid(s * density, s)

assert np.isclose(expected_s, model.F, atol=1e-3)

# ============================================================================
# 8. CALIBRATION RECOVERS SYNTHETIC SMILE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("8. CALIBRATION RECOVERS SYNTHETIC SMILE")
print("=" * 78)

"""
Generate a synthetic smile from known mixture parameters
and check that calibration reproduces the smile.
"""

F = 100.0
T = 1.0
r = 0.03

true_model = LognormalMixtureModel(F=F, T=T, r=r)

true_model.p = 0.35
true_model.displacement = -0.06
true_model.sigma1 = 0.14
true_model.sigma2 = 0.28

true_model.F1, true_model.F2 = true_model._component_forwards(true_model.p, true_model.displacement)

strikes = np.linspace(75.0, 130.0, 20)

market_vols = true_model.implied_vol(strikes)

calibrated_model = LognormalMixtureModel(F=F, T=T, r=r)

calibrated_model.calibrate(strikes, market_vols)

fitted_vols = calibrated_model.implied_vol(strikes)

assert np.max(np.abs(fitted_vols - market_vols)) < 1e-4

# ============================================================================
# 9. SURFACE CALIBRATION
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("9. SURFACE CALIBRATION")
print("=" * 78)

maturities = np.array([0.25, 0.50, 1.00, 2.00])

forwards = np.array([100.5, 101.0, 102.0, 104.0])

rates = np.array([0.03, 0.03, 0.032, 0.035])

strikes = np.array([80.0, 90.0, 100.0, 110.0, 120.0])

# Actual volatility surface
vol_surface = np.array(
    [
        [0.3000, 0.2600, 0.2300, 0.2200, 0.2300],
        [0.2800, 0.2400, 0.2100, 0.2000, 0.2100],
        [0.2500, 0.2200, 0.2000, 0.1900, 0.2000],
        [0.2300, 0.2100, 0.1900, 0.1900, 0.2000],
    ]
)

# --------------------------------------------------------
# Calibrate to the volatility surface
# --------------------------------------------------------

surface = LognormalMixtureSurface(
    maturities=maturities, forwards=forwards, strikes=strikes, vol_surface=vol_surface, rates=rates
)

surface.calibrate()

# --------------------------------------------------------
# Reconstruct fitted surface
# --------------------------------------------------------

fitted_surface = surface.fitted_vol_surface()

errors = fitted_surface - vol_surface

max_error = np.max(np.abs(errors))

rmse = np.sqrt(np.mean(errors**2))

# --------------------------------------------------------
# Regression checks
# --------------------------------------------------------

assert np.all(np.isfinite(fitted_surface))

assert np.all(fitted_surface > 0.0)

# Reasonable calibration quality
assert max_error < 0.01
assert rmse < 0.005

# Forward condition at each maturity
for i, model in enumerate(surface.models):

    assert np.isclose(model.forward_check, forwards[i], atol=1e-10)

# --------------------------------------------------------
# Diagnostics
# --------------------------------------------------------

if 1 == 0:
    print("Market surface:")
    print(vol_surface)

    print("Fitted surface:")
    print(fitted_surface)

    print("Errors:")
    print(errors)

    print("Maximum calibration error:", max_error)
    print("Surface RMSE:", rmse)

    for i, model in enumerate(surface.models):

        print(maturities[i], model.p, model.F1, model.F2, model.sigma1, model.sigma2)

