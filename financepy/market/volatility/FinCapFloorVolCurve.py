##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes

##########################################################################
# TODO: This is unfinished
##########################################################################


class FinCapFloorVolCurve():
    ''' Class to manage a term structure of cap (flat) volatilities and to
    do the conversion to caplet (spot) volatilities. This does not manage a
    strike dependency, only a term structure. For that you need the class
    FinCapFloorVolSurface. '''

    def __init__(self,
                 curveDate,  # Valuation date for cap volatility
                 capFloorStartDates,  # curve date + start dates for caplets
                 capFloorVols,  # Flat cap volatility out to c/flet start dates
                 dayCountType):
        ''' Create a cap/floor volatility curve given a curve date, a list of
        cap/floor dates and a vector of cap volatilities. To avoid confusion
        first date of the capDates must be equal to the curve date and first
        cap volatility for this date must equal zero. The internal times are
        calculated according to the provided day count convention.'''

        numTimes = len(capFloorStartDates)
        numVols = len(capFloorVols)

        if numTimes != numVols:
            raise FinError("Date and volatility vectors not same length.")

        if numTimes < 2:
            raise FinError("FinCapVolCurve requires at least two dates/vols")

        if curveDate != capFloorStartDates[0]:
            raise FinError("FinCapFloorDates must start on curve date")

        self._curveDate = curveDate

        if capFloorVols[0] != 0.0:
            raise FinError("Curve date cap floor volatility must equal zero")

        self._capFloorVols = capFloorVols

        # Basic validation of dates
        prevDt = self._curveDate
        for dt in capFloorStartDates[1:]:
            if dt < prevDt:
                raise FinError("CapFloorlet Dates not in increasing order")

            if dt == prevDt:
                raise FinError("Two successive dates are identical")

            prevDt = dt

        self._capFloorStartDates = capFloorStartDates

        if isinstance(dayCountType, FinDayCountTypes) is False:
            raise FinError("DayCountType must be of type FinDayCountTypes.")

        self._dayCountType = dayCountType

        self.generateCapFloorletVols()

###############################################################################

    def generateCapFloorletVols(self):
        ''' Bootstrap caplet volatilities from cap volatilities. '''

        self._times = []
        self._taus = []

        dayCounter = FinDayCount(self._dayCountType)

        prevDt = self._curveDate

        for dt in self._capFloorStartDates:
            t = (dt - self._curveDate) / gDaysInYear
            tau = dayCounter.yearFrac(prevDt, dt)
            self._times.append(t)
            self._taus.append(tau)
            prevDt = dt

        self._capFloorletVols = []
        fwdRateVol = self._capFloorVols[0]
        self._capFloorletVols.append(fwdRateVol)
        cumLibor2Tau = (fwdRateVol**2) * self._taus[0]

        sumTau = 0.0
        for i in range(1, len(self._capFloorStartDates)):
            t = self._times[i]
            tau = self._taus[i]
            sumTau += tau
            volCapFloor = self._capFloorVols[i]
            volLibor2 = ((volCapFloor**2) * sumTau - cumLibor2Tau)/tau
            volLibor = np.sqrt(volLibor2)
            self._capFloorletVols.append(volLibor)
            cumLibor2Tau += volLibor2 * self._taus[i]

        self._capFloorVols = np.array(self._capFloorVols)
        self._capFloorletVols = np.array(self._capFloorletVols)

###############################################################################

    def capFloorletVol(self, t):
        ''' Return the forward rate caplet/floorlet volatility for a specific
        forward caplet start date. The period of the volatility is implicitly
        the intercaplet spacing period used when creating the class object. The
        volatility interpolation is piecewise flat. '''

        numVols = len(self._times)
        vol = self._capFloorletVols[0]

        for i in range(0, numVols):
            if self._times[i] > t:
                vol = self._capFloorletVols[i]
                return vol

        return self._capFloorletVols[-1]

###############################################################################

    def capFloorVol(self, t):
        ''' Return the cap flat volatility to a specific forward start date for
        the last caplet/floorlet in the cap/floor. The volatility interpolation
        is piecewise flat. '''

        numVols = len(self._times)
        vol = self._capFloorVols[0]

        for i in range(0, numVols):
            if self._times[i] >= t:
                vol = self._capFloorVols[i]
                return vol

        return self._capFloorVols[-1]

###############################################################################

    def __repr__(self):
        ''' Output the contents of the FinCapVolCurve class object. '''

        numTimes = len(self._times)
        s = " TIME     TAU    CAP VOL    CAPLET VOL"
        for i in range(0, numTimes):
            t = self._times[i]
            tau = self._taus[i]
            volCap = self._capVols[i]
            fwdLiborVol = self._capletVols[i]
            s += labelToString("%7.4f  %6.4f  %9.4f  %9.4f"
                               % (t, tau, volCap*100.0, fwdLiborVol*100.0))

        return s

###############################################################################
