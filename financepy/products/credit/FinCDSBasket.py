##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: There are several speed ups for the Monte-Carlo including calculating
# all default baskets at the same time.

import numpy as np

from ...finutils.FinError import FinError

from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes

from ...products.credit.FinCDS import FinCDS

from ...models.FinModelGaussianCopula1F import homogeneousBasketLossDbn
from ...models.FinModelGaussianCopula import defaultTimesGC
from ...models.FinModelStudentTCopula import FinModelStudentTCopula

from ...products.credit.FinCDSCurve import FinCDSCurve

from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...market.curves.FinInterpolator import interpolate, FinInterpTypes

from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString

###############################################################################
# TODO: Convert functions to use NUMBA!!
###############################################################################


class FinCDSBasket(object):

    ''' Class to deal with n-to-default CDS baskets. '''

    def __init__(self,
                 stepInDate: FinDate,
                 maturityDate: FinDate,
                 notional: float = ONE_MILLION,
                 runningCoupon: float = 0.0,
                 longProtection: bool = True,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):

        checkArgumentTypes(self.__init__, locals())

        self._stepInDate = stepInDate
        self._maturityDate = maturityDate
        self._notional = notional
        self._runningCoupon = runningCoupon / 10000.0
        self._longProtection = longProtection
        self._dayCountType = dayCountType
        self._dateGenRuleType = dateGenRuleType
        self._calendarType = calendarType
        self._freqType = freqType
        self._busDayAdjustType = busDayAdjustType

        self._cdsContract = FinCDS(self._stepInDate,
                                   self._maturityDate,
                                   self._runningCoupon,
                                   1.0,
                                   self._longProtection,
                                   self._freqType,
                                   self._dayCountType,
                                   self._calendarType,
                                   self._busDayAdjustType,
                                   self._dateGenRuleType)

###############################################################################

    def valueLegs_MC(self,
                     valuationDate,
                     nToDefault,
                     defaultTimes,
                     issuerCurves,
                     liborCurve):
        ''' Value the legs of the default basket using Monte Carlo. The default
        times are an input so this valuation is not model dependent. '''

        numCredits = defaultTimes.shape[0]
        numTrials = defaultTimes.shape[1]

        adjustedDates = self._cdsContract._adjustedDates
        numFlows = len(adjustedDates)
        dayCount = FinDayCount(self._dayCountType)

        averageAccrualFactor = 0.0

        rpv01ToTimes = np.zeros(numFlows)
        for iTime in range(1, numFlows):
            t = (adjustedDates[iTime] - valuationDate) / gDaysInYear
            dt0 = adjustedDates[iTime - 1]
            dt1 = adjustedDates[iTime]
            accrualFactor = dayCount.yearFrac(dt0, dt1)[0]
            averageAccrualFactor += accrualFactor
            rpv01ToTimes[iTime] = rpv01ToTimes[iTime - 1] + \
                accrualFactor * liborCurve._df(t)

        averageAccrualFactor /= numFlows

        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        rpv01 = 0.0
        prot = 0.0

        assetTau = np.zeros(numCredits)

        for iTrial in range(0, numTrials):
            for iCredit in range(0, numCredits):
                assetTau[iCredit] = defaultTimes[iCredit, iTrial]

            # ORDER THE DEFAULT TIMES
            assetTau.sort()

            # GET THE Nth DEFAULT TIME
            minTau = assetTau[nToDefault - 1]

            if minTau < tmat:
                numPaymentsIndex = int(minTau / averageAccrualFactor)
                rpv01Trial = rpv01ToTimes[numPaymentsIndex]
                rpv01Trial += (minTau - numPaymentsIndex *
                               averageAccrualFactor)

                # DETERMINE IDENTITY OF N-TO-DEFAULT CREDIT IF BASKET NOT HOMO
                assetIndex = 0
                for iCredit in range(0, numCredits):
                    if minTau == defaultTimes[iCredit, iTrial]:
                        assetIndex = iCredit
                        break

                protTrial = (1.0 - issuerCurves[assetIndex]._recoveryRate)
                protTrial *= liborCurve._df(minTau)

            else:

                numPaymentsIndex = int(tmat / averageAccrualFactor)
                rpv01Trial = rpv01ToTimes[numPaymentsIndex]
                protTrial = 0.0

            rpv01 += rpv01Trial
            prot += protTrial

        rpv01 = rpv01 / numTrials
        prot = prot / numTrials
        return (rpv01, prot)

###############################################################################

    def valueGaussian_MC(self,
                         valuationDate,
                         nToDefault,
                         issuerCurves,
                         correlationMatrix,
                         liborCurve,
                         numTrials,
                         seed):
        ''' Value the default basket using a Gaussian copula model. This
        depends on the issuer curves and correlation matrix. '''

        numCredits = len(issuerCurves)

        if nToDefault > numCredits or nToDefault < 1:
            raise FinError("nToDefault must be 1 to numCredits")

        defaultTimes = defaultTimesGC(issuerCurves,
                                      correlationMatrix,
                                      numTrials,
                                      seed)

        rpv01, protPV = self.valueLegs_MC(valuationDate,
                                          nToDefault,
                                          defaultTimes,
                                          issuerCurves,
                                          liborCurve)

        spd = protPV / rpv01
        value = self._notional * (protPV - self._runningCoupon * rpv01)

        if not self._longProtection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def valueStudentT_MC(self,
                         valuationDate,
                         nToDefault,
                         issuerCurves,
                         correlationMatrix,
                         degreesOfFreedom,
                         liborCurve,
                         numTrials,
                         seed):
        ''' Value the default basket using the Student-T copula. '''

        numCredits = len(issuerCurves)

        if nToDefault > numCredits or nToDefault < 1:
            raise FinError("nToDefault must be 1 to numCredits")

        model = FinModelStudentTCopula()

        defaultTimes = model.defaultTimes(issuerCurves,
                                          correlationMatrix,
                                          degreesOfFreedom,
                                          numTrials,
                                          seed)

        rpv01, protPV = self.valueLegs_MC(valuationDate,
                                          nToDefault,
                                          defaultTimes,
                                          issuerCurves,
                                          liborCurve)

        spd = protPV / rpv01
        value = self._notional * (protPV - self._runningCoupon * rpv01)

        if not self._longProtection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def value1FGaussian_Homo(self,
                             valuationDate,
                             nToDefault,
                             issuerCurves,
                             betaVector,
                             liborCurve,
                             numPoints=50):
        ''' Value default basket using 1 factor Gaussian copula and analytical
        approach which is only exact when all recovery rates are the same. '''

        numCredits = len(issuerCurves)

        if numCredits == 0:
            raise FinError("Num Credits is zero")

        if nToDefault < 1 or nToDefault > numCredits:
            raise FinError("NToDefault must be 1 to numCredits")

        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        if tmat < 0.0:
            raise FinError("Value date is after maturity date")

        paymentDates = self._cdsContract._adjustedDates
        numTimes = len(paymentDates)

        issuerSurvivalProbabilities = np.zeros(numCredits)
        recoveryRates = np.zeros(numCredits)
        basketTimes = np.zeros(numTimes)
        basketSurvivalCurve = np.zeros(numTimes)

        basketTimes[0] = 0.0
        basketSurvivalCurve[0] = 1.0

        for iTime in range(1, numTimes):

            t = (paymentDates[iTime] - valuationDate) / gDaysInYear

            for iCredit in range(0, numCredits):
                issuerCurve = issuerCurves[iCredit]
                recoveryRates[iCredit] = issuerCurve._recoveryRate
                issuerSurvivalProbabilities[iCredit] = interpolate(
                    t, issuerCurve._times, issuerCurve._values,
                    FinInterpTypes.FLAT_FWD_RATES.value)

            lossDbn = homogeneousBasketLossDbn(issuerSurvivalProbabilities,
                                               recoveryRates,
                                               betaVector,
                                               numPoints)

            basketSurvivalCurve[iTime] = 1.0
            for iToDefault in range(nToDefault, numCredits + 1):
                basketSurvivalCurve[iTime] -= lossDbn[iToDefault]

            basketTimes[iTime] = t

        curveRecovery = recoveryRates[0]
        liborCurve = issuerCurves[0]._liborCurve
        basketCurve = FinCDSCurve(valuationDate, [], liborCurve, curveRecovery)
        basketCurve._times = basketTimes
        basketCurve._values = basketSurvivalCurve

        protLegPV = self._cdsContract.protectionLegPV(
            valuationDate, basketCurve, curveRecovery)
        riskyPV01 = self._cdsContract.riskyPV01(valuationDate, basketCurve)['clean_rpv01']

        # Long protection
        mtm = self._notional * (protLegPV - riskyPV01 * self._runningCoupon)

        if not self._longProtection:
            mtm *= -1.0

        basketOutput = np.zeros(4)
        basketOutput[0] = mtm
        basketOutput[1] = riskyPV01 * self._notional * self._runningCoupon
        basketOutput[2] = protLegPV * self._notional
        basketOutput[3] = protLegPV / riskyPV01

        return basketOutput

###############################################################################

    def __repr__(self):
        ''' print out details of the CDS contract and all of the calculated
        cashflows '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("STEP-IN DATE", self._stepInDate)
        s += labelToString("MATURITY", self._maturityDate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("RUNNING COUPON", self._runningCoupon*10000, "bp\n")
        s += labelToString("DAYCOUNT", self._dayCountType)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUSDAYRULE", self._busDayAdjustType)
        s += labelToString("DATEGENRULE", self._dateGenRuleType)

#        header = "PAYMENT_DATE, YEAR_FRAC, FLOW"
#        valueTable = [self._adjustedDates, self._accrualFactors, self._flows]
#        precision = "12.6f"
#        s += tableToString(header, valueTable, precision)

        return s

###############################################################################
