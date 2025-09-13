# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes

from ..utils.math import normcdf

# NOTE: Keeping this separate from SABR for the moment.

########################################################################################


class BlackShifted:
    """Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. This model also allows
    the distribution to be shifted to the negative in order to allow for
    negative interest rates."""

    ####################################################################################

    def __init__(self, volatility: float, shift: float, implementation: int = 0) -> None:
        """Create FinModel black using parameters."""
        self.volatility = volatility
        self.shift = shift
        self.implementation = 0
        self.num_steps = 0
        self.seed = 0
        self.param1 = 0
        self.param2 = 0

    ####################################################################################

    def value(
        self,
        forward_rate: float,  # Forward rate
        strike_rate: float,  # Strike Rate
        time_to_expiry: float,  # time to expiry in years
        df: float,  # Discount Factor to expiry date
        call_or_put: OptionTypes,  # Call or put
    ) -> float:
        """Price a derivative using Black's model which values in the forward
        measure following a change of measure. The sign of the shift is the
        same as Matlab."""

        s = self.shift
        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        sqrt_t = np.sqrt(t)
        vol = self.volatility

        d1 = np.log((f + s) / (k + s)) + vol * vol * t / 2
        d1 = d1 / (vol * sqrt_t)
        d2 = d1 - vol * sqrt_t

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * ((f + s) * normcdf(d1) - (k + s) * normcdf(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * ((k + s) * normcdf(-d2) - (f + s) * normcdf(-d1))

        raise ValueError("Option type must be European Call(C) or Put(P)")

    ####################################################################################

    def __repr__(self) -> str:

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        s += label_to_string("SHIFT", self.shift)
        s += label_to_string("IMPLEMENTATION", self.implementation)
        s += label_to_string("NUMSTEPS", self.num_steps)
        return s
