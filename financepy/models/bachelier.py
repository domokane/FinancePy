# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from scipy.stats import norm

from ..utils.error import FinError
from ..utils.global_types import OptionTypes
from ..utils.helpers import label_to_string

# NOTE: Need to convert option types to use enums.
# NOTE: Perhaps just turn this into a function rather than a class.

########################################################################################


class Bachelier:
    """Bachelier's Model which prices call and put options in the forward
    measure assuming the underlying rate follows a normal process.
    """

    ####################################################################################

    def __init__(self, volatility: float) -> None:
        """Create FinModel black using parameters."""

        if volatility <= 0.0:
            raise FinError("Volatility must be positive")

        self.volatility = volatility

    ####################################################################################

    def value(
        self,
        forward_rate: float,  # Forward rate F
        strike_rate: float,  # Strike Rate K
        time_to_expiry: float,  # Time to Expiry (years)
        df: float,  # Discount Factor to expiry date
        call_or_put: OptionTypes,  # Call or put
    ) -> float:
        """Price a call or put option using Bachelier's model."""

        # Vectorisations
        f = np.asarray(forward_rate, dtype=float)
        k = np.asarray(strike_rate, dtype=float)
        t = np.asarray(time_to_expiry, dtype=float)
        df = np.asarray(df, dtype=float)

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be non-negative")

        if np.any(df <= 0.0):
            raise FinError("Discount factor must be positive")

        v = self.volatility
        root_t = np.sqrt(np.maximum(t, 0.0))
        std = v * root_t
        d = np.where(std > 0.0, (f - k) / std, 0.0)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            values = df * ((f - k) * norm.cdf(d) + std * norm.pdf(d))
            intrinsic = df * np.maximum(f - k, 0.0)
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            values = df * ((k - f) * norm.cdf(-d) + std * norm.pdf(d))
            intrinsic = df * np.maximum(k - f, 0.0)
        else:
            raise FinError("Option type must be call or put")

        values = np.where(t <= 0.0, intrinsic, values)
        return values

    ####################################################################################

    def __repr__(self) -> str:

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        return s
