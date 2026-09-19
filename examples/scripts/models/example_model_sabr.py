# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.models.sabr import vol_function_sabr
from financepy.models.sabr import SABR
from financepy.utils.global_types import OptionTypes
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - Model Sabr
# ============================================================================



alpha = 0.28
beta = 1.0
rho = -0.09
nu = 0.21

f = 0.043
k = 0.050
t = 2.0

########################################################################################




########################################################################################




########################################################################################

# ============================================================================
# 1. SABR
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. SABR")
print("=" * 78)

print("ALPHA", "BETA", "RHO", "VOL")

for alpha in [0.1, 0.2, 0.3]:
    for beta in [0.5, 1.0, 2.0]:
        for rho in [-0.8, 0.0, 0.8]:
            params = np.array([alpha, beta, rho, nu])
            vol = vol_function_sabr(params, f, k, t)
            print(alpha, beta, rho, vol)

# ============================================================================
# 2. SABR  CALIBRATION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("2. SABR  CALIBRATION")
print("=" * 78)

beta = 0.5
rho = -0.09
nu = 0.1

strike_vol = 0.1

f = 0.043
k = 0.050
r = 0.03
t_exp = 2.0

call_option_type = OptionTypes.EUROPEAN_CALL
put_option_type = OptionTypes.EUROPEAN_PUT

df = np.exp(-r * t_exp)

print("TEST", "CALIBRATION ERROR")

# Make SABR equivalent to lognormal (Black) model
# (i.e. alpha = 0, beta = 1, rho = 0, nu = 0, shift = 0)
model_sabr_01 = SABR(0.0, 1.0, 0.0, 0.0)
model_sabr_01.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

implied_lognormal_vol = model_sabr_01.black_vol(f, k, t_exp)
implied_atm_lognormal_vol = model_sabr_01.black_vol(k, k, t_exp)
implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol

assert implied_lognormal_smile == 0.0, "In lognormal model, smile should be flat"
calibration_error = round(strike_vol - implied_lognormal_vol, 12)
print("LOGNORMAL CASE", calibration_error)

# Volatility: pure SABR dynamics
model_sabr_02 = SABR(alpha, beta, rho, nu)
model_sabr_02.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

implied_lognormal_vol = model_sabr_02.black_vol(f, k, t_exp)
implied_atm_lognormal_vol = model_sabr_02.black_vol(k, k, t_exp)
implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol
calibration_error = round(strike_vol - implied_lognormal_vol, 12)
print("SABR CASE", calibration_error)

# Valuation: pure SABR dynamics
value_call = model_sabr_02.value(f, k, t_exp, df, call_option_type)
value_put = model_sabr_02.value(f, k, t_exp, df, put_option_type)
assert round(value_call - value_put, 12) == round(
    df * (f - k), 12
), "The method called 'value()' doesn't comply with Call-Put parity"

# =============================================================================
# 3. VISUALISE THE SABR VOLATILITY SMILE
# =============================================================================
# The original nested loop compares parameter combinations numerically. Here we
# hold alpha, beta and nu fixed and vary rho. Multiple lines are clearer than a
# surface plot and show how correlation changes smile/skew across strikes.
plot_strikes = np.linspace(0.025, 0.070, 80)
plt.figure()
for plot_rho in [-0.8, 0.0, 0.8]:
    plot_params = np.array([0.2, 0.5, plot_rho, nu])
    plot_vols = [vol_function_sabr(plot_params, f, x, t) for x in plot_strikes]
    plt.plot(plot_strikes * 100.0, np.array(plot_vols) * 100.0, label=f"rho = {plot_rho:+.1f}")
plt.axvline(f * 100.0, linestyle="--", label="Forward")
plt.xlabel("Strike (%)")
plt.ylabel("Implied volatility (%)")
plt.title("SABR implied-volatility smile for different rho values")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
