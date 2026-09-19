# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.models.black import Black
from financepy.utils.global_types import OptionTypes
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - Black
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. BLACK
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# Measures first-order sensitivity of value to the underlying market variable.
# Measures how delta itself changes as the underlying market variable changes.
# Measures sensitivity of value to the passage of time.
# Measures sensitivity of value to volatility.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. BLACK")
print("=" * 78)

forward = 0.034
strike = 0.050
risk_free_ir = 0.00
t_exp = 2.0
volatility = 0.20

print("ITEM", "CALL", "PUT")

call_option_type = OptionTypes.EUROPEAN_CALL
put_option_type = OptionTypes.EUROPEAN_PUT

df = np.exp(-risk_free_ir * t_exp)
model = Black(volatility)

dp = 12  # Precision

try:

    value_call = model.value(forward, strike, t_exp, df, call_option_type)
    value_put = model.value(forward, strike, t_exp, df, put_option_type)

    assert round((value_call - value_put), dp) == round(
        df * (forward - strike), dp
    ), "The method called 'value()' doesn't comply with Call-Put parity"

    print("VALUE", value_call, value_put)

    delta_call = model.delta(forward, strike, t_exp, df, call_option_type)
    delta_put = model.delta(forward, strike, t_exp, df, put_option_type)

    assert (
        round((1 / df) * (delta_call - delta_put), dp) == 1.0
    ), "The method called 'delta()' doesn't comply with Call-put parity"

    print("DELTA", delta_call, delta_put)

    gamma_call = model.gamma(forward, strike, t_exp, df, call_option_type)
    gamma_put = model.gamma(forward, strike, t_exp, df, put_option_type)

    assert (
        round(gamma_call - gamma_put, dp) == 0.0
    ), "The method called 'gamma()' doesn't comply with Call-Put parity"

    print("GAMMA", gamma_call, gamma_put)

    theta_call = model.theta(forward, strike, t_exp, df, call_option_type)
    theta_put = model.theta(forward, strike, t_exp, df, put_option_type)

    assert round((theta_call - theta_put), dp) == round(
        (risk_free_ir * t_exp) * (forward - strike) * df, dp
    ), "The method called 'theta()' doesn't comply with Call-Put parity"

    print("THETA", theta_call, theta_put)

    vega_call = model.vega(forward, strike, t_exp, df, call_option_type)
    vega_put = model.vega(forward, strike, t_exp, df, put_option_type)

    assert (
        round(vega_call - vega_put, dp) == 0.0
    ), "The method called 'vega()' doesn't comply with Call-Put parity"

    print("VEGA", vega_call, vega_put)

except AssertionError as err:
    raise err

# =============================================================================
# 2. VISUALISE OPTION VALUE VERSUS FORWARD
# =============================================================================
# Moving the forward through the strike shows the familiar call/put behaviour:
# calls become more valuable as the forward rises, while puts become less valuable.
forwards = np.linspace(0.01, 0.09, 81)
call_values = [model.value(x, strike, t_exp, df, call_option_type) for x in forwards]
put_values = [model.value(x, strike, t_exp, df, put_option_type) for x in forwards]

plt.figure()
plt.plot(forwards * 100.0, call_values, label="Call")
plt.plot(forwards * 100.0, put_values, label="Put")
plt.axvline(strike * 100.0, linestyle="--", label="Strike")
plt.xlabel("Forward rate (%)")
plt.ylabel("Option value")
plt.title("Black option value versus forward")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
