"""Closed-form European power options under geometric Brownian motion.

Payoff: max(S_T**power - strike, 0) or max(strike - S_T**power, 0).
The carry rate is dividend yield for equity or the foreign risk-free rate
for FX. The discount rate is the domestic risk-free rate for FX.
"""
import numpy as np
from scipy.special import ndtr

from ..utils.error import FinError
from ..utils.global_types import OptionTypes


def _validate(spot, strike, power, maturity, rate, carry, volatility, opt_type):
    if not np.isfinite([spot, strike, power, maturity, rate, carry, volatility]).all():
        raise FinError("Power option inputs must be finite.")
    if spot <= 0.0 or strike < 0.0 or maturity < 0.0 or volatility < 0.0:
        raise FinError("Require spot>0, strike>=0, T>=0 and volatility>=0.")
    if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        raise FinError("Only European power calls/puts are supported.")


def power_option_value(spot, strike, power, maturity, rate, carry, volatility,
                       opt_type):
    """Exact discounted expectation of a European power payoff under GBM."""
    _validate(spot, strike, power, maturity, rate, carry, volatility, opt_type)
    discount = np.exp(-rate * maturity)
    if maturity == 0.0:
        terminal = spot**power
        payoff = (max(terminal - strike, 0.0) if opt_type == OptionTypes.EUROPEAN_CALL
                  else max(strike - terminal, 0.0))
        return discount * payoff

    mean_power = spot**power * np.exp(
        power * (rate - carry) * maturity
        + 0.5 * power * (power - 1.0) * volatility**2 * maturity
    )
    sigma_power = abs(power) * volatility * np.sqrt(maturity)
    if sigma_power == 0.0:
        payoff = (max(mean_power - strike, 0.0) if opt_type == OptionTypes.EUROPEAN_CALL
                  else max(strike - mean_power, 0.0))
        return discount * payoff
    if strike == 0.0:
        return discount * mean_power if opt_type == OptionTypes.EUROPEAN_CALL else 0.0

    d2 = (np.log(spot**power / strike)
          + power * (rate - carry - 0.5 * volatility**2) * maturity) / sigma_power
    d1 = d2 + sigma_power
    if opt_type == OptionTypes.EUROPEAN_CALL:
        return discount * (mean_power * ndtr(d1) - strike * ndtr(d2))
    return discount * (strike * ndtr(-d2) - mean_power * ndtr(-d1))


def power_option_delta(spot, strike, power, maturity, rate, carry, volatility,
                       opt_type):
    """Analytical spot derivative; at a deterministic kink returns 0 by convention."""
    _validate(spot, strike, power, maturity, rate, carry, volatility, opt_type)
    if power == 0.0:
        return 0.0

    discount = np.exp(-rate * maturity)
    mean_power = spot**power * np.exp(
        power * (rate - carry) * maturity
        + 0.5 * power * (power - 1.0) * volatility**2 * maturity
    )
    derivative = discount * power * mean_power / spot
    sigma_power = abs(power) * volatility * np.sqrt(maturity)
    if sigma_power == 0.0:
        if opt_type == OptionTypes.EUROPEAN_CALL:
            return derivative if mean_power > strike else 0.0
        return -derivative if mean_power < strike else 0.0
    if strike == 0:
        return derivative if opt_type == OptionTypes.EUROPEAN_CALL else 0.0

    d1 = ((np.log(spot**power / strike)
          + power * (rate - carry - 0.5 * volatility**2) * maturity)
          / sigma_power + sigma_power)
    return derivative * (ndtr(d1) if opt_type == OptionTypes.EUROPEAN_CALL
                         else -ndtr(-d1))


class PowerOption:
    """Reusable European power option for equity and FX under GBM."""

    def __init__(self, strike, power, opt_type):
        _validate(1.0, strike, power, 1.0, 0.0, 0.0, 0.0, opt_type)
        self.strike, self.power, self.opt_type = strike, power, opt_type

    def value(self, spot, maturity, rate, carry, volatility):
        return power_option_value(
            spot, self.strike, self.power, maturity, rate, carry,
            volatility, self.opt_type
        )

    def delta(self, spot, maturity, rate, carry, volatility):
        return power_option_delta(
            spot, self.strike, self.power, maturity, rate, carry,
            volatility, self.opt_type
        )

