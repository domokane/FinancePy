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
        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        root_t = np.sqrt(t)
        v = self.volatility
        d = (f - k) / (v * root_t)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * ((f - k) * norm.cdf(d) + v * root_t * norm.pdf(d))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * ((k - f) * norm.cdf(-d) + v * root_t * norm.pdf(d))
        else:
            raise FinError("Option type must be a European Call(C) or Put(P)")

    ####################################################################################

    def __repr__(self) -> str:

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        return s

