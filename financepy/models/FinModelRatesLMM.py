##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, jit

from ..finutils.FinError import FinError
from ..finutils.FinMath import accruedInterpolator
from ..market.curves.FinInterpolate import FinInterpMethods, uinterpolate
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinDayCount import FinDayCount, FinDayCountTypes
from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorTypes
from ..finutils.FinGlobalVariables import gDaysInYear

interp = FinInterpMethods.FLAT_FORWARDS.value

small = 1e-10

##########################################################################
# LIBOR MARKET MODEL
##########################################################################


@njit(fastmath=True, cache=True)
def simulateLibors(taus, initialLibors, sigmas, numPaths):
    ''' Simulate Libor rates in the Spot Measure '''

    # Number of points we evolve on the forward curve
    numFwdRates = len(taus)

    # We assume the same time steps as the forward rates
    numTimeSteps = numFwdRates

    if len(sigmas) != numFwdRates:
        raise FinError("Number of fwd libor vols not same as num libor fwds")

    # 3D matrix which stores the Libor forwards
    F = np.zeros((numPaths, numFwdRates, numTimeSteps))

    # Generate independent Gaussian RVs
    Z = np.random.rand(numPaths, numTimeSteps)

    # Initialise forward curves with today's values
    for iPath in range(0, numPaths):
        for iFwd in range(0, numFwdRates):
            F[iPath, iFwd, 0] = initialLibors[iFwd]

    sqrtdt = np.zeros(numTimeSteps)
    for iTime in range(0, numTimeSteps):
        sqrtdt[iTime] = np.sqrt(taus[iTime])

    # Simulate in the spot measure
    for iPath in range(0, numPaths):
        for iTime in range(0, numTimeSteps-1):
            for iFwd in range(iTime+1, numFwdRates):
                sigmai = sigmas[iFwd-iTime-1]
                mu = 0.0
                for kFwd in range(iTime+1, iFwd+1):  # i is kFwd
                    tau = taus[kFwd]
                    fwd = F[iPath, kFwd, iTime]
                    sigmak = sigmas[kFwd-iTime-1]
                    mu += sigmak * tau * fwd / (1.0 + tau * fwd)

                dt = taus[iTime]
                S = sigmai * (mu - (1/2) * (sigmai)) * dt
                S = S + sigmai * Z[iPath, iTime] * sqrtdt[iTime]
                S = np.exp(S)
                F[iPath, iFwd, iTime+1] = F[iPath, iFwd, iTime] * S

    return F

###############################################################################
###############################################################################
# LIBOR MARKET MODEL
###############################################################################
###############################################################################


class FinModelRatesLMM():

    def __init__(self,
                 valuationDate,
                 fwdCurveDates,
                 dayCountType):
        ''' Constructs the Libor Market model. '''

        self._valuationDate = valuationDate
        self._numPaths = 0
        self._dayCountType = dayCountType

        prevDt = fwdCurveDates[0]
        for dt in fwdCurveDates[1:]:
            if dt <= prevDt:
                raise FinError("Fwd Curve Dates need to be increasing.")
            prevDt = dt

        self._fwdCurveDates = fwdCurveDates

        dayCounter = FinDayCount(dayCountType)

        self._accrualFactors = []
        self._fwdEndTimes = []
        self._fwdStartTimes = []

        prevDt = fwdCurveDates[0]

        for dt in fwdCurveDates[1:]:

            accFac = dayCounter.yearFrac(prevDt, dt)
            self._accrualFactors.append(accFac)

            prevTime = (prevDt - self._valuationDate) / gDaysInYear
            nextTime = (dt - self._valuationDate) / gDaysInYear

            self._fwdStartTimes.append(prevTime)
            self._fwdEndTimes.append(nextTime)

            prevDt = dt

        self._accrualFactors = np.array(self._accrualFactors)
        self._fwdStartTimes = np.array(self._fwdStartTimes)
        self._fwdEndTimes = np.array(self._fwdEndTimes)
        self._numForwards = len(self._accrualFactors)

###############################################################################

    def simulateLibors(self,
                       liborCurve,
                       capletVolatilityCurve,
                       numPaths,
                       seed):

        initialLiborRates = []
        prevDt = self._valuationDate

        for dt in self._fwdCurveDates[1:]:
            liborRate = liborCurve.fwdRate(prevDt, dt, self._dayCountType)
            initialLiborRates.append(liborRate)
            prevDt = dt

        initialLiborRates = np.array(initialLiborRates)

        np.random.seed(seed)

        self._numPaths = numPaths

        fwdVols = capletVolatilityCurve._fwdLiborVols

        self._libors = simulateLibors(self._accrualFactors,
                                      initialLiborRates,
                                      fwdVols,
                                      numPaths)

        return self._libors

###############################################################################

    def bondPrice(self, tmat):

        if tmat < 0.0:
            raise FinError("Bond price input tmat is negative.")

        if tmat > self._fwdEndTimes[-1]:
            raise FinError("Bond price input tmat is beyond end of time grid.")

        numForwards = len(self._fwdEndTimes)

        numPeriods = -999

        for iMaturity in range(0, numForwards):
            if tmat < self._fwdEndTimes[iMaturity]:
                numPeriods = iMaturity
                break

        if numPeriods == -999:
            numPeriods = numForwards

        v = 0.0

        # In the spot measure we calculate the rolled up MM account
        for iPath in range(0, self._numPaths):
            N = 1.0
            for iTime in range(0, numPeriods):
                acc = self._accrualFactors[iTime]
                libor = self._libors[iPath, iTime, iTime]
                N = N * (1.0 + acc * libor)
            v += 1.0/N

        v /= self._numPaths

        return v

###############################################################################

    def valueCapFloor(self, strike, optionType, notional):

        capFloorValue = 0.0

        if optionType == FinLiborCapFloorTypes.CAP:

            for p in range(0, self._numPaths):
                df = 1.0
                for i in range(0, len(self._accFactors)):
                    acc = self._accFactors[i]
                    libor = self._libors[p, i, i]
                    df = df / (1.0 + acc * libor)
                    capletValue = acc * max(0.0, libor - strike) * df
                    print(p, i, acc, libor, df, capletValue)

                    capFloorValue += capletValue
        elif optionType == FinLiborCapFloorTypes.FLOOR:
            pass
        else:
            raise FinError("Unknown CapFloorType")

        capFloorValue /= self._numPaths
        capFloorValue *= notional
        return capFloorValue

###############################################################################

    def dump(self):
        for iPath in range(0, self._numPaths):
            for iTime in range(0, self._numForwards):
                s = ""
                for iFwd in range(0, self._numForwards):
                    s += ("%9.5f " % self._libors[iPath, iFwd, iTime])

                print("iPath:", iPath, "iTime:", iTime, "FWDS:", s)
            print("")

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = "LMM Model\n"
        s += labelToString("Sigma", self._sigma)
        s += labelToString("a", self._a)
        return s

###############################################################################

# # @njit(fastmath=True, cache=True)
# def simulateLiborsHull(taus, initialLibors, capletVolCurve, numPaths):

#     if len(taus) != len(initialLibors):
#         raise FinError("Taus not same length as initial Libors")

#     v = capletVolCurve._fwdGammaVols

#     if len(taus) != len(v):
#         raise FinError("Taus not same as sigmas")

#     R = [[1.000],
#          [0.924, 1.000],
#          [0.707, 0.924, 1.000],
#          [0.557, 0.833, 0.981, 1.000],
#          [0.454, 0.760, 0.951, 0.997, 1.000],
#          [0.760, 0.951, 0.997, 0.963, 0.924, 1.000],
#          [0.843, 0.985, 0.976, 0.916, 0.862, 0.990, 1.000],
#          [0.837, 0.983, 0.979, 0.921, 0.867, 0.992, 1.000, 1.00],
#          [0.837, 0.983, 0.979, 0.920, 0.867, 0.992, 1.000, 1.000, 1.00],
#          [0.920, 1.000, 0.928, 0.838, 0.767, 0.954, 0.986, 0.985, 0.985, 1.000]
#          ]

#     RR = np.zeros((10,10))
#     for i in range(0,10):
#         for j in range(0,i+1):
#             print("i:", i, "j:", j, "==>", RR[i,j])
#             RR[i,j] = R[i][j]
#             RR[j,i] = R[i][j]

#  #   print(RR)

#     numeraire = 0

#     numMaturities = len(taus)
#     numTimeSteps = len(taus)

#     ###########################################################################
#     # LIBORs is a 3D matrix with indices
#     # [iPath, iTime, iMaturity]
#     # The index iPath is the path
#     # The index iTime is the current time of the forward rate
#     # The index iMaturity is the point on the forward curve
#     ###########################################################################

#     F = np.zeros((numPaths, numMaturities, numMaturities))
#     Z = np.random.rand(numPaths, numMaturities, numMaturities)

#     # initialise Libors to initial term structure
#     for iPath in range(0, numPaths):
#         for i in range(0, numMaturities):
#             F[iPath, i, 0] = initialLibors[i]

#     sqrtTaus = np.zeros(numMaturities)
#     for j in range(0, numMaturities):
#         sqrtTaus[j] = np.sqrt(taus[j])

#     M = numPaths
#     N = numTimeSteps
#     m = numMaturities
#     dt = 1.0

#     # Loop over simulation paths
#     for p in range(0, numPaths):
#         drift = 0
#         mu = 0.0
#         mu1 = 0.0
#         # Loop over time steps
#         for l in range(0, N):
#             # Loop over curve points
#             for k in range(1, m+1):
#                 z = np.random.rand()

#                 if k < numeraire:

#                     for j in range(k+1, numeraire):
#                         mu += taus[j]*v[k]*F[j-1]/(1+taus[j]*F[j-1])*dt
#                 else:

#                     for j in range(numeraire, k+1):
#                         mu1 += taus[j]*v[k]*F[j-1]/(1+taus[j]*F[j-1])*dt

#             drift = -mu + mu1
#             logf = np.log(v[k]*drift - 0.5 * v[k]*v[k]*dt + v[k] * z * sqrtTaus[j])
#             F[p,l,k] = F[p,l,k-1] * np.exp(logf)

#     return F