##############################################################################
# Copyright (C) 2026 Dominic O'Kane
##############################################################################

import numpy as np
from ..utils.global_vars import G_SMALL
from ..utils.error import FinError
from ..utils.global_types import OptionTypes
from ..utils.global_types import DigitalOptionTypes
from ..utils.math import normcdf_vect


def bs_digital_option_value(s0, t, K, r, q, vol, call_put, cash_asset):

    ln_s0_k = np.log(s0 / K)
    sqrt_t = np.sqrt(t)

    if abs(vol) < G_SMALL:
        vol = G_SMALL

    d1 = ln_s0_k + (r - q + vol * vol / 2.0) * t
    d1 = d1 / vol / sqrt_t
    d2 = d1 - vol * sqrt_t

    v = 0.0

    if cash_asset == DigitalOptionTypes.CASH_OR_NOTHING.value:

        if call_put == OptionTypes.EUROPEAN_CALL.value:
            v = np.exp(-r * t) * normcdf_vect(d2)
        elif call_put == OptionTypes.EUROPEAN_PUT.value:
            v = np.exp(-r * t) * normcdf_vect(-d2)

    elif cash_asset == DigitalOptionTypes.ASSET_OR_NOTHING.value:

        if call_put == OptionTypes.EUROPEAN_CALL.value:
            v = s0 * np.exp(-q * t) * normcdf_vect(d1)
        elif call_put == OptionTypes.EUROPEAN_PUT.value:
            v = s0 * np.exp(-q * t) * normcdf_vect(-d1)
    else:
        raise FinError("Unknown digital option type.")

    return v
