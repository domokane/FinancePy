# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union

import numpy as np

from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.error import FinError
from ...utils.helpers import check_argument_types, label_to_string
from .discount_curve import DiscountCurve
from .interpolator import InterpTypes


########################################################################################
class FXDiscountCurve(DiscountCurve):
    """Discount curve implied from FX spot, FX forwards, and a domestic curve."""

    ####################################################################################
    def __init__(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        forward_dts: list,
        forward_fx_rates: Union[list, np.ndarray],
        domestic_curve: DiscountCurve,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):
        """Create an implied foreign discount curve from FX forward quotes.

        The FX rate convention is domestic currency per unit of foreign currency.
        Given spot S, forward F(T), and domestic discount factor DF_dom(T), the
        implied foreign discount factor is:

            DF_for(T) = F(T) / S * DF_dom(T)
        """
        check_argument_types(self.__init__, locals())

        if spot_fx_rate <= 0.0:
            raise FinError("Spot FX rate must be positive.")

        if len(forward_dts) == 0:
            raise FinError("At least one FX forward date is required.")

        if len(forward_dts) != len(forward_fx_rates):
            raise FinError("Forward dates and FX forward rates must have same length.")

        if domestic_curve.value_dt != value_dt:
            raise FinError("Domestic curve value date must match curve value date.")

        for dt in forward_dts:
            if dt <= value_dt:
                raise FinError("Forward dates must be after the value date.")

        forward_fx_rates = np.asarray(forward_fx_rates, dtype=float)

        if np.any(~np.isfinite(forward_fx_rates)) or np.any(forward_fx_rates <= 0.0):
            raise FinError("FX forward rates must be finite and positive.")

        domestic_dfs = domestic_curve.df(forward_dts)
        foreign_dfs = forward_fx_rates / spot_fx_rate * domestic_dfs

        self.spot_fx_rate = spot_fx_rate
        self.forward_fx_rates = forward_fx_rates
        self.domestic_curve = domestic_curve

        super().__init__(
            value_dt,
            forward_dts,
            foreign_dfs,
            interp_type,
            time_dc_type,
        )

    ####################################################################################
    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SPOT FX RATE", self.spot_fx_rate)
        s += "\n"
        s += super().__repr__()
        return s

    ####################################################################################
    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
