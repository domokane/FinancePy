##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinError import FinError
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes

##########################################################################
# TODO: This is unfinished


class FinCapVolCurve():
    ''' Class to manage a term structure of cap (flat) volatilities and to 
    do the conversion to caplet (spot) volatilities. '''

    def __init__(self,
                 curveDate,
                 capDates,
                 capVols,  # Flat vols out to the cap maturity 
                 dayCountType):

        numTimes = len(capDates)
        numVols = len(capVols)

        if numTimes != numVols:
            raise FinError("Date and volatility vectors not same length.")

        # Basic validation of dates
        prevDt = curveDate
        for dt in capDates:
            if dt < prevDt:
                raise FinError("Caplet Dates not in increasing order")

            if dt == prevDt:
                raise FinError("Two successive dates are identical")

            prevDt = dt

        if isinstance(dayCountType, FinDayCountTypes) is False:
            raise FinError("DayCountType must be of type FinDayCountTypes.")

        # Begin work
        self._curveDate = curveDate
        self._times = []
        self._taus = []
        self._capVols = capVols

        dayCounter = FinDayCount(dayCountType)

        prevDt = curveDate
        for dt in capDates:
            t = (dt - curveDate) / gDaysInYear
            self._times.append(t)
            tau = dayCounter.yearFrac(prevDt, dt)
            self._taus.append(tau)
            prevDt = dt

        self._capletVols = []
        fwdRateVol = capVols[0]
        self._capletVols.append(fwdRateVol)
        cumLibor2Tau = (fwdRateVol**2) * self._taus[0]

        for i in range(1, numVols):
            t = self._times[i]
            tau = self._taus[i]
            volCap = self._capVols[i]
            volLibor2 = ((volCap**2) * t - cumLibor2Tau)/tau
            volLibor = np.sqrt(volLibor2)
            self._capletVols.append(volLibor)
            cumLibor2Tau += volLibor2 * self._taus[i]

        self._capVols = np.array(self._capletVols)
        self._capletVols = np.array(self._capletVols)

###############################################################################

    def capletVol(self, t):
        ''' Return the forward rate vol at a specific forward. '''

        numVols = len(self._times)
        vol = self._capletVols[0]

        for i in range(0, numVols):
            if self._times[i] > t:
                vol = self._capletVols[i]
                return vol

        return self._capletVols[-1]

###############################################################################

    def capVol(self, t):
        ''' Return the forward rate vol at a specific forward. '''

        numVols = len(self._times)
        vol = self._capVols[0]

        for i in range(0, numVols):
            if self._times[i] >= t:
                vol = self._capVols[i]
                return vol

        return self._capVols[-1]

###############################################################################

    def dump(self):
        numTimes = len(self._times)
        print(" Time     tau    CapletVol    FwdVol")
        for i in range(0, numTimes):
            t = self._times[i]
            tau = self._taus[i]
            volCap = self._capVols[i]
            fwdLiborVol = self._capletVols[i]
            print("%7.4f  %6.4f  %9.4f  %9.4f"
                  %(t, tau, volCap*100.0, fwdLiborVol*100.0))

###############################################################################
