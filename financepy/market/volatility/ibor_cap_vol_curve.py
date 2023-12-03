##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.helpers import label_to_string
from ...utils.global_vars import gDaysInYear
from ...utils.day_count import DayCount, DayCountTypes

##########################################################################
# TODO: Calibration
# TODO: Interpolation
# TODO: Integration kernel for LMM
##########################################################################


class IborCapVolCurve():
    """ Class to manage a term structure of cap (flat) volatilities and to
    do the conversion to caplet (spot) volatilities. This does not manage a
    strike dependency, only a term structure. The cap and caplet volatilies
    are keyed off the cap and caplet maturity dates. However this volatility
    only applies to the evolution of the Ibor rate out to the caplet start
    dates. Note also that this class also handles floor vols."""

    def __init__(self,
                 curve_dt,  # Valuation date for cap volatility
                 cap_maturity_dts,  # curve date + maturity dates for caps
                 cap_sigmas,  # Flat cap volatility for cap maturity dates
                 dc_type):
        """ Create a cap/floor volatility curve given a curve date, a list of
        cap maturity dates and a vector of cap volatilities. To avoid confusion
        first date of the capDates must be equal to the curve date and first
        cap volatility for this date must equal zero. The internal times are
        calculated according to the provided day count convention. Note cap and
        floor volatilities are the same for the same strike and tenor, I just
        refer to cap volatilities in the code for code simplicity. """

        num_times = len(cap_maturity_dts)
        num_vols = len(cap_sigmas)

        if num_times != num_vols:
            raise FinError("Date and volatility vectors not same length.")

        if num_times < 2:
            raise FinError("FinCapVolCurve requires at least two dates/vols")

        if curve_dt != cap_maturity_dts[0]:
            raise FinError("FinCapFloorDates must start on curve date")

        self._curve_dt = curve_dt

        if cap_sigmas[0] != 0.0:
            raise FinError("Curve date cap floor volatility must equal zero")

        self._cap_sigmas = np.array(cap_sigmas)
        self._caplet_gammas = []

        # Basic validation of dates
        prev_dt = self._curve_dt
        for dt in cap_maturity_dts[1:]:
            if dt < prev_dt:
                raise FinError("CapFloorLet Dates not in increasing order")

            if dt == prev_dt:
                raise FinError("Two successive dates are identical")

            prev_dt = dt

        self._cap_maturity_dts = cap_maturity_dts

        if isinstance(dc_type, DayCountTypes) is False:
            raise FinError("DayCountType must be of type DayCountTypes.")

        self._dc_type = dc_type

        self.generate_caplet_vols()

###############################################################################

    def generate_caplet_vols(self):
        """ Bootstrap caplet volatilities from cap volatilities using similar
        notation to Hull's book (page 32.11). The first volatility in the
        vector of caplet vols is zero. """

        self._times = []
        self._taus = []

        day_counter = DayCount(self._dc_type)
        prev_dt = self._curve_dt
        numCaps = len(self._cap_maturity_dts)

        for dt in self._cap_maturity_dts:
            t = (dt - self._curve_dt) / gDaysInYear
            self._times.append(t)
            tau = day_counter.year_frac(prev_dt, dt)[0]
            self._taus.append(tau)
            prev_dt = dt

        fwd_rateVol = self._cap_sigmas[0]
        self._caplet_gammas = np.zeros(numCaps)
        self._caplet_gammas[0] = 0.0
        cumIbor2Tau = (fwd_rateVol**2) * self._taus[0]

        sumTau = 0.0
        for i in range(1, len(self._cap_maturity_dts)):
            t = self._times[i]
            tau = self._taus[i]
            sumTau += tau
            volCap = self._cap_sigmas[i]
            volIbor2 = ((volCap**2) * sumTau - cumIbor2Tau) / tau

            if volIbor2 < 0.0:
                raise FinError("Error due to negative caplet variance.")

            volIbor = np.sqrt(volIbor2)
            self._caplet_gammas[i] = volIbor
            cumIbor2Tau += volIbor2 * self._taus[i]

###############################################################################

    def caplet_vol(self, dt):
        """ Return the forward rate caplet/floorlet volatility for a specific
        forward caplet expiry date. The period of the volatility is the
        the intercaplet spacing period used when creating the class object.
        The volatility interpolation is piecewise flat. """

        if isinstance(dt, Date):
            t = (dt - self._curve_dt) / gDaysInYear
        else:
            t = dt

        if t <= self._times[1]:
            return self._caplet_gammas[1]

        if 1 == 0:
            print(self._times)
            print(self._caplet_gammas)
            print(t)

        num_vols = len(self._times)
        vol = self._caplet_gammas[1]

        for i in range(1, num_vols):
            if self._times[i] >= t:
                vol = self._caplet_gammas[i]
                return vol

        return self._caplet_gammas[-1]

###############################################################################

    def cap_vol(self, dt):
        """ Return the cap flat volatility for a specific cap maturity date for
        the last caplet/floorlet in the cap/floor. The volatility interpolation
        is piecewise flat. """

        if isinstance(dt, Date):
            t = (dt - self._curve_dt) / gDaysInYear
        else:
            t = dt

        num_vols = len(self._times)
        vol = self._cap_sigmas[0]

        if 1 == 0:
            print(self._times)
            print(self._caplet_gammas)
            print(t)

        for i in range(1, num_vols):
            if self._times[i] >= t:
                vol = self._cap_sigmas[i]
                return vol

        return self._cap_sigmas[-1]

###############################################################################

    def __repr__(self):
        """ Output the contents of the FinCapVolCurve class object. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        num_times = len(self._times)
        s += " TIME     TAU    CAP VOL    CAPLET VOL"
        for i in range(0, num_times):
            t = self._times[i]
            tau = self._taus[i]
            volCap = self._cap_sigmas[i]
            fwdIborVol = self._capletVols[i]
            s += label_to_string("%7.4f  %6.4f  %9.4f  %9.4f"
                                 % (t, tau, volCap*100.0, fwdIborVol*100.0))

        return s

###############################################################################
