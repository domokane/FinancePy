
# Allow this example to run directly from its category folder.
import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.models.merton_jump_diffusion import MertonJumpDiffusion
from financepy.models.black_scholes_analytic import value

# ============================================================================
# FINANCEPY EXAMPLES - MertonJumpDiffusion
# ============================================================================




TOL = 1.0e-10


###############################################################################
# Pricing tests
###############################################################################




###############################################################################




###############################################################################




###############################################################################




###############################################################################




###############################################################################
# Implied volatility tests
###############################################################################




###############################################################################




###############################################################################




###############################################################################
# Calibration tests
###############################################################################




###############################################################################




###############################################################################
# Behavioural tests
###############################################################################




###############################################################################
# Parameter validation
###############################################################################

# ============================================================================
# 1. MERTON REDUCES TO BLACK SCHOLES ZERO INTENSITY
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. MERTON REDUCES TO BLACK SCHOLES ZERO INTENSITY")
print("=" * 78)

s = 100.0
t = 1.25
k = 105.0
r = 0.04
q = 0.015
sigma = 0.22

model = MertonJumpDiffusion(
    sigma=sigma,
    jump_intensity=0.0,
    jump_mean=-0.10,
    jump_volatility=0.25,
)

for option_type in [
    OptionTypes.EUROPEAN_CALL,
    OptionTypes.EUROPEAN_PUT,
]:

    v_merton = model.value(
        s,
        t,
        k,
        r,
        q,
        option_type,
    )

    v_bs = value(
        s,
        t,
        k,
        r,
        q,
        sigma,
        option_type.value,
    )

    assert abs(v_merton - v_bs) < TOL

# ============================================================================
# 2. MERTON REDUCES TO BLACK SCHOLES ZERO JUMP SIZE
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. MERTON REDUCES TO BLACK SCHOLES ZERO JUMP SIZE")
print("=" * 78)

s = 100.0
t = 0.75
k = 95.0
r = 0.035
q = 0.01
sigma = 0.20

model = MertonJumpDiffusion(
    sigma=sigma,
    jump_intensity=2.0,
    jump_mean=0.0,
    jump_volatility=0.0,
)

for option_type in [
    OptionTypes.EUROPEAN_CALL,
    OptionTypes.EUROPEAN_PUT,
]:

    v_merton = model.value(
        s,
        t,
        k,
        r,
        q,
        option_type,
    )

    v_bs = value(
        s,
        t,
        k,
        r,
        q,
        sigma,
        option_type.value,
    )

    assert abs(v_merton - v_bs) < 1.0e-9

# ============================================================================
# 3. MERTON PUT CALL PARITY
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("3. MERTON PUT CALL PARITY")
print("=" * 78)

s = 100.0
t = 1.50
k = 110.0
r = 0.045
q = 0.02

model = MertonJumpDiffusion(
    sigma=0.18,
    jump_intensity=0.80,
    jump_mean=-0.12,
    jump_volatility=0.25,
)

call = model.value(
    s,
    t,
    k,
    r,
    q,
    OptionTypes.EUROPEAN_CALL,
)

put = model.value(
    s,
    t,
    k,
    r,
    q,
    OptionTypes.EUROPEAN_PUT,
)

parity = s * np.exp(-q * t) - k * np.exp(-r * t)

assert abs(call - put - parity) < 1.0e-9

# ============================================================================
# 4. MERTON PRICE POSITIVE
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("4. MERTON PRICE POSITIVE")
print("=" * 78)

model = MertonJumpDiffusion(
    sigma=0.20,
    jump_intensity=1.0,
    jump_mean=-0.15,
    jump_volatility=0.30,
)

call = model.value(
    100.0,
    1.0,
    120.0,
    0.04,
    0.01,
    OptionTypes.EUROPEAN_CALL,
)

put = model.value(
    100.0,
    1.0,
    80.0,
    0.04,
    0.01,
    OptionTypes.EUROPEAN_PUT,
)

assert call > 0.0
assert put > 0.0

# ============================================================================
# 5. MERTON CALL MONOTONIC IN STRIKE
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("5. MERTON CALL MONOTONIC IN STRIKE")
print("=" * 78)

model = MertonJumpDiffusion(
    sigma=0.18,
    jump_intensity=0.75,
    jump_mean=-0.12,
    jump_volatility=0.25,
)

strikes = np.array(
    [
        70.0,
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
        130.0,
    ]
)

values = np.array(
    [
        model.value(
            100.0,
            1.0,
            k,
            0.04,
            0.02,
            OptionTypes.EUROPEAN_CALL,
        )
        for k in strikes
    ]
)

assert np.all(np.diff(values) < 0.0)

# ============================================================================
# 6. MERTON IMPLIED VOL BLACK SCHOLES LIMIT
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("6. MERTON IMPLIED VOL BLACK SCHOLES LIMIT")
print("=" * 78)

sigma = 0.235

model = MertonJumpDiffusion(
    sigma=sigma,
    jump_intensity=0.0,
    jump_mean=-0.10,
    jump_volatility=0.20,
)

vol = model.implied_volatility(
    100.0,
    1.0,
    105.0,
    0.04,
    0.015,
    OptionTypes.EUROPEAN_CALL,
)

assert abs(vol - sigma) < 1.0e-8

# ============================================================================
# 7. MERTON SMILE DIMENSIONS
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("7. MERTON SMILE DIMENSIONS")
print("=" * 78)

model = MertonJumpDiffusion(
    sigma=0.18,
    jump_intensity=0.8,
    jump_mean=-0.10,
    jump_volatility=0.25,
)

strikes = np.array(
    [
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
    ]
)

vols = model.volatility_smile(
    100.0,
    1.0,
    strikes,
    0.04,
    0.02,
    OptionTypes.EUROPEAN_CALL,
)

assert vols.shape == strikes.shape
assert np.all(np.isfinite(vols))
assert np.all(vols > 0.0)

# ============================================================================
# 8. MERTON SURFACE DIMENSIONS
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("8. MERTON SURFACE DIMENSIONS")
print("=" * 78)

model = MertonJumpDiffusion(
    sigma=0.18,
    jump_intensity=0.8,
    jump_mean=-0.10,
    jump_volatility=0.25,
)

strikes = np.array(
    [
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
    ]
)

expiries = np.array(
    [
        0.25,
        0.50,
        1.00,
        2.00,
    ]
)

vols = model.volatility_surface(
    100.0,
    strikes,
    expiries,
    0.04,
    0.02,
    OptionTypes.EUROPEAN_CALL,
)

assert vols.shape == (
    len(expiries),
    len(strikes),
)

assert np.all(np.isfinite(vols))
assert np.all(vols > 0.0)

# ============================================================================
# 9. MERTON CALIBRATION RECOVERS SYNTHETIC SURFACE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("9. MERTON CALIBRATION RECOVERS SYNTHETIC SURFACE")
print("=" * 78)

s = 100.0
r = 0.04
q = 0.015

strikes = np.array(
    [
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
    ]
)

expiries = np.array(
    [
        0.25,
        0.50,
        1.00,
        2.00,
    ]
)

true_model = MertonJumpDiffusion(
    sigma=0.18,
    jump_intensity=0.70,
    jump_mean=-0.12,
    jump_volatility=0.22,
)

market_vols = true_model.volatility_surface(
    s,
    strikes,
    expiries,
    r,
    q,
    OptionTypes.EUROPEAN_CALL,
)

result = MertonJumpDiffusion.calibrate(
    stock_price=s,
    strikes=strikes,
    expiries=expiries,
    market_vols=market_vols,
    risk_free_rates=r,
    dividend_yields=q,
    option_type=OptionTypes.EUROPEAN_CALL,
    initial_guess=np.array(
        [
            0.20,
            0.50,
            -0.10,
            0.20,
        ]
    ),
    calibration_type="VOL",
)

assert result.success

assert abs(result.sigma - true_model.sigma) < 1.0e-4

assert abs(result.jump_intensity - true_model.jump_intensity) < 1.0e-3

assert abs(result.jump_mean - true_model.jump_mean) < 1.0e-4

assert abs(result.jump_volatility - true_model.jump_volatility) < 1.0e-4

assert result.rmse < 1.0e-7

# ============================================================================
# 10. MERTON VEGA CALIBRATION RECOVERS SURFACE
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("10. MERTON VEGA CALIBRATION RECOVERS SURFACE")
print("=" * 78)

s = 100.0
r = 0.04
q = 0.015

strikes = np.array(
    [
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
    ]
)

expiries = np.array(
    [
        0.25,
        0.50,
        1.00,
        2.00,
    ]
)

true_model = MertonJumpDiffusion(
    sigma=0.19,
    jump_intensity=0.60,
    jump_mean=-0.10,
    jump_volatility=0.20,
)

market_vols = true_model.volatility_surface(
    s,
    strikes,
    expiries,
    r,
    q,
    OptionTypes.EUROPEAN_CALL,
)

result = MertonJumpDiffusion.calibrate(
    stock_price=s,
    strikes=strikes,
    expiries=expiries,
    market_vols=market_vols,
    risk_free_rates=r,
    dividend_yields=q,
    option_type=OptionTypes.EUROPEAN_CALL,
    initial_guess=np.array(
        [
            0.20,
            0.50,
            -0.08,
            0.18,
        ]
    ),
    calibration_type="VEGA",
)

assert result.success
assert result.rmse < 1.0e-6

# ============================================================================
# 11. NEGATIVE JUMPS GENERATE EQUITY SKEW
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("11. NEGATIVE JUMPS GENERATE EQUITY SKEW")
print("=" * 78)

model = MertonJumpDiffusion(
    sigma=0.17,
    jump_intensity=1.0,
    jump_mean=-0.15,
    jump_volatility=0.20,
)

vol_low_strike = model.implied_volatility(
    100.0,
    0.50,
    80.0,
    0.03,
    0.01,
    OptionTypes.EUROPEAN_CALL,
)

vol_high_strike = model.implied_volatility(
    100.0,
    0.50,
    120.0,
    0.03,
    0.01,
    OptionTypes.EUROPEAN_CALL,
)

assert vol_low_strike > vol_high_strike

# ============================================================================
# 12. MERTON INVALID PARAMETERS
# ============================================================================
# What this section demonstrates:
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("12. MERTON INVALID PARAMETERS")
print("=" * 78)

try:
    MertonJumpDiffusion(
        sigma=-0.1,
        jump_intensity=1.0,
        jump_mean=-0.1,
        jump_volatility=0.2,
    )

    assert False

except ValueError:
    pass

try:
    MertonJumpDiffusion(
        sigma=0.2,
        jump_intensity=-1.0,
        jump_mean=-0.1,
        jump_volatility=0.2,
    )

    assert False

except ValueError:
    pass

try:
    MertonJumpDiffusion(
        sigma=0.2,
        jump_intensity=1.0,
        jump_mean=-0.1,
        jump_volatility=-0.2,
    )

    assert False

except ValueError:
    pass

