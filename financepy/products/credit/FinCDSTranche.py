##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Add __repr__ method

import numpy as np
from math import sqrt


from ...models.FinModelGaussianCopula1F import trSurvProbGaussian
from ...models.FinModelGaussianCopula1F import trSurvProbAdjBinomial
from ...models.FinModelGaussianCopula1F import trSurvProbRecursion
from ...models.FinModelGaussianCopulaLHP import trSurvProbLHP

from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes

from ...products.credit.FinCDS import FinCDS
from ...products.credit.FinCDSCurve import FinCDSCurve

from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...market.curves.FinInterpolator import FinInterpTypes, interpolate
from ...finutils.FinError import FinError

from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDate import FinDate

###############################################################################

from enum import Enum


class FinLossDistributionBuilder(Enum):
    RECURSION = 1
    ADJUSTED_BINOMIAL = 2
    GAUSSIAN = 3
    LHP = 4

###############################################################################


class FinCDSTranche(object):

    def __init__(self,
                 stepInDate: FinDate,
                 maturityDate: FinDate,
                 k1: float,
                 k2: float,
                 notional: float = ONE_MILLION,
                 runningCoupon: float = 0.0,
                 longProtection: bool = True,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):

        checkArgumentTypes(self.__init__, locals())

        if k1 >= k2:
            raise FinError("K1 must be less than K2")

        self._k1 = k1
        self._k2 = k2

        self._stepInDate = stepInDate
        self._maturityDate = maturityDate
        self._notional = notional
        self._runningCoupon = runningCoupon
        self._longProtection = longProtection
        self._dayCountType = dayCountType
        self._dateGenRuleType = dateGenRuleType
        self._calendarType = calendarType
        self._freqType = freqType
        self._busDayAdjustType = busDayAdjustType

        notional = 1.0

        self._cdsContract = FinCDS(self._stepInDate,
                                   self._maturityDate,
                                   self._runningCoupon,
                                   notional,
                                   self._longProtection,
                                   self._freqType,
                                   self._dayCountType,
                                   self._calendarType,
                                   self._busDayAdjustType,
                                   self._dateGenRuleType)

###############################################################################

    def valueBC(self,
                valuationDate,
                issuerCurves,
                upfront,
                runningCoupon,
                corr1,
                corr2,
                numPoints=50,
                model=FinLossDistributionBuilder.RECURSION):

        numCredits = len(issuerCurves)
        k1 = self._k1
        k2 = self._k2
        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        if tmat < 0.0:
            raise FinError("Value date is after maturity date")

        if abs(k1 - k2) < 0.00000001:
            output = np.zeros(4)
            output[0] = 0.0
            output[1] = 0.0
            output[2] = 0.0
            output[3] = 0.0
            return output

        if k1 > k2:
            raise FinError("K1 > K2")

        kappa = k2 / (k2 - k1)

        recoveryRates = np.zeros(numCredits)

        paymentDates = self._cdsContract._adjustedDates
        numTimes = len(paymentDates)

        beta1 = sqrt(corr1)
        beta2 = sqrt(corr2)
        betaVector1 = np.zeros(numCredits)
        for bb in range(0, numCredits):
            betaVector1[bb] = beta1

        betaVector2 = np.zeros(numCredits)
        for bb in range(0, numCredits):
            betaVector2[bb] = beta2

        qVector = np.zeros(numCredits)
        qt1 = np.zeros(numTimes)
        qt2 = np.zeros(numTimes)
        trancheTimes = np.zeros(numTimes)
        trancheSurvivalCurve = np.zeros(numTimes)

        trancheTimes[0] = 0
        trancheSurvivalCurve[0] = 1.0
        qt1[0] = 1.0
        qt2[0] = 1.0

        for i in range(1, numTimes):

            t = (paymentDates[i] - valuationDate) / gDaysInYear

            for j in range(0, numCredits):

                issuerCurve = issuerCurves[j]
                vTimes = issuerCurve._times
                qRow = issuerCurve._values
                recoveryRates[j] = issuerCurve._recoveryRate
                qVector[j] = interpolate(
                    t, vTimes, qRow, FinInterpTypes.FLAT_FWD_RATES.value)

            if model == FinLossDistributionBuilder.RECURSION:
                qt1[i] = trSurvProbRecursion(
                    0.0, k1, numCredits, qVector, recoveryRates,
                    betaVector1, numPoints)
                qt2[i] = trSurvProbRecursion(
                    0.0, k2, numCredits, qVector, recoveryRates,
                    betaVector2, numPoints)
            elif model == FinLossDistributionBuilder.ADJUSTED_BINOMIAL:
                qt1[i] = trSurvProbAdjBinomial(
                    0.0, k1, numCredits, qVector, recoveryRates,
                    betaVector1, numPoints)
                qt2[i] = trSurvProbAdjBinomial(
                    0.0, k2, numCredits, qVector, recoveryRates,
                    betaVector2, numPoints)
            elif model == FinLossDistributionBuilder.GAUSSIAN:
                qt1[i] = trSurvProbGaussian(
                    0.0,
                    k1,
                    numCredits,
                    qVector,
                    recoveryRates,
                    betaVector1,
                    numPoints)
                qt2[i] = trSurvProbGaussian(
                    0.0,
                    k2,
                    numCredits,
                    qVector,
                    recoveryRates,
                    betaVector2,
                    numPoints)
            elif model == FinLossDistributionBuilder.LHP:
                qt1[i] = trSurvProbLHP(
                    0.0, k1, numCredits, qVector, recoveryRates, beta1)
                qt2[i] = trSurvProbLHP(
                    0.0, k2, numCredits, qVector, recoveryRates, beta2)
            else:
                raise FinError(
                    "Unknown model type only full and AdjBinomial allowed")

            if qt1[i] > qt1[i - 1]:
                raise FinError(
                    "Tranche K1 survival probabilities not decreasing.")

            if qt2[i] > qt2[i - 1]:
                raise FinError(
                    "Tranche K2 survival probabilities not decreasing.")

            trancheSurvivalCurve[i] = kappa * qt2[i] + (1.0 - kappa) * qt1[i]
            trancheTimes[i] = t

        curveRecovery = 0.0  # For tranches only
        liborCurve = issuerCurves[0]._liborCurve
        trancheCurve = FinCDSCurve(
            valuationDate, [], liborCurve, curveRecovery)
        trancheCurve._times = trancheTimes
        trancheCurve._values = trancheSurvivalCurve

        protLegPV = self._cdsContract.protectionLegPV(
            valuationDate, trancheCurve, curveRecovery)
        riskyPV01 = self._cdsContract.riskyPV01(valuationDate, trancheCurve)['clean_rpv01']

        mtm = self._notional * (protLegPV - upfront - riskyPV01 * runningCoupon)

        if not self._longProtection:
            mtm *= -1.0

        trancheOutput = np.zeros(4)
        trancheOutput[0] = mtm
        trancheOutput[1] = riskyPV01 * self._notional * runningCoupon
        trancheOutput[2] = protLegPV * self._notional
        trancheOutput[3] = protLegPV / riskyPV01

        return trancheOutput

###############################################################################
