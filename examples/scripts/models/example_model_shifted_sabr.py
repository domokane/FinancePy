# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.models.sabr_shifted import SABRShifted
from financepy.utils.global_types import OptionTypes

# ============================================================================
# FINANCEPY EXAMPLES - SABRShifted
# ============================================================================



########################################################################################



    # TODO: adding Call-Put parity test for all sensitivities


########################################################################################

# ============================================================================
# 1. SHIFTED SABR
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. SHIFTED SABR")
print("=" * 78)

print("TEST", "CALIBRATION ERROR")

alpha = 0.0
beta = 0.5
rho = -0.09
nu = 0.1
shift = 0.02

strike_vol = 0.1

f = 0.043
k = 0.050
r = 0.03
t_exp = 2.0

call_option_type = OptionTypes.EUROPEAN_CALL
put_option_type = OptionTypes.EUROPEAN_PUT

df = np.exp(-r * t_exp)

# SABR equivalent to lognormal (Black) model (i.e. beta = 1, rho = 0, nu = 0, shift = 0)
model_sabr_01 = SABRShifted(0.0, 1.0, 0.0, 0.0, 0.0)
model_sabr_01.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

implied_lognormal_vol = model_sabr_01.black_vol(f, k, t_exp)
implied_atm_lognormal_vol = model_sabr_01.black_vol(k, k, t_exp)
implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol

assert implied_lognormal_smile == 0.0, "In lognormal model, smile should be flat"
calibration_error = round(strike_vol - implied_lognormal_vol, 12)
print("LOGNORMAL CASE", calibration_error)

# Volatility: pure SABR dynamics
model_sabr_02 = SABRShifted(alpha, beta, rho, nu, shift)
model_sabr_02.set_alpha_from_black_vol(strike_vol, f, k, t_exp)

implied_lognormal_vol = model_sabr_02.black_vol(f, k, t_exp)
implied_atm_lognormal_vol = model_sabr_02.black_vol(k, k, t_exp)
implied_lognormal_smile = implied_lognormal_vol - implied_atm_lognormal_vol
calibration_error = round(strike_vol - implied_lognormal_vol, 12)
print("SABR CASE", calibration_error)

# Valuation: pure SABR dynamics
value_call = model_sabr_02.value(f, k, t_exp, df, call_option_type)
value_put = model_sabr_02.value(f, k, t_exp, df, put_option_type)
assert round((value_call - value_put), 12) == round(
    df * (f - k), 12
), "The method called 'value()' doesn't comply with Call-Put parity"

