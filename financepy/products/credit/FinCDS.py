##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from numba import njit, float64, int64
from math import exp, log
from copy import deepcopy

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinCalendar import FinCalendar, FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinHelperFunctions import labelToString, tableToString
from ...market.curves.FinInterpolator import FinInterpTypes, _uinterpolate

from ...finutils.FinHelperFunctions import checkArgumentTypes

useFlatHazardRateIntegral = True
standardRecovery = 0.40

###############################################################################
# TODO: Perform protection leg pv analytically using fact that hazard rate and
#       interest rates are flat between their combined node points. Right now I
#       do not find the protection leg PV calculations to be a bottleneck,
#       especially given the speedup benefits of using NUMBA.
###############################################################################


@njit(float64[:](float64, float64, float64[:], float64[:], float64[:],
                 float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def _riskyPV01_NUMBA(teff,
                     accrualFactorPCDToNow,
                     paymentTimes,
                     yearFracs,
                     npIborTimes,
                     npIborValues,
                     npSurvTimes,
                     npSurvValues,
                     pv01Method):
    ''' Fast calculation of the risky PV01 of a CDS using NUMBA.
    The output is a numpy array of the full and clean risky PV01.'''

    method = FinInterpTypes.FLAT_FWD_RATES.value

    couponAccruedIndicator = 1

    # Method 0 : This is the market standard which assumes that the coupon
    # accrued is treated as though on average default occurs roughly midway
    # through a coupon period.

    tncd = paymentTimes[1]

    # The first coupon is a special case which needs to be handled carefully
    # taking into account what coupon has already accrued and what has not
    qeff = _uinterpolate(teff, npSurvTimes, npSurvValues, method)
    q1 = _uinterpolate(tncd, npSurvTimes, npSurvValues, method)
    z1 = _uinterpolate(tncd, npIborTimes, npIborValues, method)

    # this is the part of the coupon accrued from previous coupon date to now
    # accrualFactorPCDToNow = dayCount.yearFrac(pcd,teff)

    # reference credit survives to the premium payment date
    fullRPV01 = q1 * z1 * yearFracs[1]

    # coupon accrued from previous coupon to today paid in full at default
    # before coupon payment
    fullRPV01 = fullRPV01 + z1 * \
        (qeff - q1) * accrualFactorPCDToNow * couponAccruedIndicator

    # future accrued from now to coupon payment date assuming default roughly
    # midway
    fullRPV01 += 0.5 * z1 * \
        (qeff - q1) * (yearFracs[1] - accrualFactorPCDToNow) \
        * couponAccruedIndicator

    for it in range(2, len(paymentTimes)):

        t2 = paymentTimes[it]

        q2 = _uinterpolate(t2, npSurvTimes, npSurvValues, method)
        z2 = _uinterpolate(t2, npIborTimes, npIborValues, method)

        accrualFactor = yearFracs[it]

        # full coupon is paid at the end of the current period if survives to
        # payment date
        fullRPV01 += q2 * z2 * accrualFactor

        #######################################################################

        if couponAccruedIndicator == 1:

            if useFlatHazardRateIntegral:
                # This needs to be updated to handle small h+r
                tau = accrualFactor
                h12 = -log(q2 / q1) / tau
                r12 = -log(z2 / z1) / tau
                alpha = h12 + r12
                expTerm = 1.0 - exp(-alpha * tau) - alpha * \
                    tau * exp(-alpha * tau)
                dfullRPV01 = q1 * z1 * h12 * \
                    expTerm / abs(alpha * alpha + 1e-20)
            else:
                dfullRPV01 = 0.50 * (q1 - q2) * z2 * accrualFactor

            fullRPV01 = fullRPV01 + dfullRPV01

        #######################################################################

        q1 = q2

    cleanRPV01 = fullRPV01 - accrualFactorPCDToNow

    return np.array([fullRPV01, cleanRPV01])

###############################################################################

@njit(float64(float64, float64, float64[:], float64[:], float64[:], float64[:],
              float64, int64, int64), fastmath=True, cache=True)
def _protectionLegPV_NUMBA(teff,
                           tmat,
                           npIborTimes,
                           npIborValues,
                           npSurvTimes,
                           npSurvValues,
                           contractRecovery,
                           numStepsPerYear,
                           protMethod):
    ''' Fast calculation of the CDS protection leg PV using NUMBA to speed up
    the numerical integration over time. '''

    method = FinInterpTypes.FLAT_FWD_RATES.value
    dt = (tmat - teff) / numStepsPerYear
    t = teff
    z1 = _uinterpolate(t, npIborTimes, npIborValues, method)
    q1 = _uinterpolate(t, npSurvTimes, npSurvValues, method)

    protPV = 0.0
    small = 1e-8

    if useFlatHazardRateIntegral is True:

        for _ in range(0, numStepsPerYear):

            t = t + dt
            z2 = _uinterpolate(t, npIborTimes, npIborValues, method)
            q2 = _uinterpolate(t, npSurvTimes, npSurvValues, method)
            # This needs to be updated to handle small h+r
            h12 = -log(q2 / q1) / dt
            r12 = -log(z2 / z1) / dt
            expTerm = exp(-(r12 + h12) * dt)
            dprotPV = h12 * (1.0 - expTerm) * q1 * z1 / \
                (abs(h12 + r12) + small)
            protPV += dprotPV
            q1 = q2
            z1 = z2

    else:

        for _ in range(0, numStepsPerYear):

            t += dt
            z2 = _uinterpolate(t, npIborTimes, npIborValues, method)
            q2 = _uinterpolate(t, npSurvTimes, npSurvValues, method)
            dq = q1 - q2
            dprotPV = 0.5 * (z1 + z2) * dq
            protPV += dprotPV
            q1 = q2
            z1 = z2

    protPV = protPV * (1.0 - contractRecovery)
    return protPV

###############################################################################
###############################################################################
###############################################################################


class FinCDS(object):
    ''' A class which manages a Credit Default Swap. It performs schedule
    generation and the valuation and risk management of CDS. '''

    def __init__(self,
                 stepInDate: FinDate,  # Date protection starts
                 maturityDateOrTenor: (FinDate, str),  # FinDate or tenor
                 runningCoupon: float,  # Annualised coupon on premium fee leg
                 notional: float = ONE_MILLION,
                 longProtection: bool = True,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a CDS from the step-in date, maturity date and coupon '''

        checkArgumentTypes(self.__init__, locals())

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            # To get the next CDS date we move on by the tenor and then roll to
            # the next CDS date after that. We do not holiday adjust it. That
            # is handled in the schedule generation.
            maturityDate = stepInDate.addTenor(maturityDateOrTenor)
            maturityDate = maturityDate.nextCDSDate()

        if stepInDate > maturityDate:
            raise FinError("Step in date after maturity date")

        self._stepInDate = stepInDate
        self._maturityDate = maturityDate
        self._runningCoupon = runningCoupon
        self._notional = notional
        self._longProtection = longProtection
        self._dayCountType = dayCountType
        self._dateGenRuleType = dateGenRuleType
        self._calendarType = calendarType
        self._freqType = freqType
        self._busDayAdjustType = busDayAdjustType

        self._generateAdjustedCDSPaymentDates()
        self._calcFlows()

###############################################################################

    def _generateAdjustedCDSPaymentDates(self):
        ''' Generate CDS payment dates which have been holiday adjusted.'''
        frequency = FinFrequency(self._freqType)
        calendar = FinCalendar(self._calendarType)
        startDate = self._stepInDate
        endDate = self._maturityDate

        self._adjustedDates = []
        numMonths = int(12.0 / frequency)

        unadjustedScheduleDates = []

        if self._dateGenRuleType == FinDateGenRuleTypes.BACKWARD:

            nextDate = endDate
            flowNum = 0

            while nextDate > startDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(-numMonths)
                flowNum += 1

            # Add on the Previous Coupon Date
            unadjustedScheduleDates.append(nextDate)
            flowNum += 1

            # reverse order
            for i in range(0, flowNum):
                dt = unadjustedScheduleDates[flowNum - i - 1]
                self._adjustedDates.append(dt)

            # holiday adjust dates except last one
            for i in range(0, flowNum - 1):

                dt = calendar.adjust(self._adjustedDates[i],
                                     self._busDayAdjustType)

                self._adjustedDates[i] = dt

            finalDate = self._adjustedDates[flowNum - 1]

            # Final date is moved forward by one day
            self._adjustedDates[flowNum - 1] = finalDate.addDays(1)

        elif self._dateGenRuleType == FinDateGenRuleTypes.FORWARD:

            nextDate = startDate
            flowNum = 0

            unadjustedScheduleDates.append(nextDate)
            flowNum = 1

            while nextDate < endDate:
                unadjustedScheduleDates.append(nextDate)
                nextDate = nextDate.addMonths(numMonths)
                flowNum = flowNum + 1

            for i in range(1, flowNum):

                dt = calendar.adjust(unadjustedScheduleDates[i],
                                     self._busDayAdjustType)

                self._adjustedDates.append(dt)

            finalDate = endDate.addDays(1)
            self._adjustedDates.append(finalDate)

        else:
            raise FinError("Unknown FinDateGenRuleType:" +
                           str(self._dateGenRuleType))

###############################################################################

    def _calcFlows(self):
        ''' Calculate cash flow amounts on premium leg. '''
        paymentDates = self._adjustedDates
        dayCount = FinDayCount(self._dayCountType)

        self._accrualFactors = []
        self._flows = []

        self._accrualFactors.append(0.0)
        self._flows.append(0.0)

        numFlows = len(paymentDates)

        for it in range(1, numFlows):
            t0 = paymentDates[it - 1]
            t1 = paymentDates[it]
            accrualFactor = dayCount.yearFrac(t0, t1)[0]
            flow = accrualFactor * self._runningCoupon * self._notional

            self._accrualFactors.append(accrualFactor)
            self._flows.append(flow)

###############################################################################

    def value(self,
              valuationDate,
              issuerCurve,
              contractRecovery=standardRecovery,
              pv01Method=0,
              prot_method=0,
              numStepsPerYear=25):
        ''' Valuation of a CDS contract on a specific valuation date given
        an issuer curve and a contract recovery rate.'''

        rpv01 = self.riskyPV01(valuationDate,
                               issuerCurve,
                               pv01Method)

        fullRPV01 = rpv01['full_rpv01']
        cleanRPV01 = rpv01['clean_rpv01']

        protPV = self.protectionLegPV(valuationDate,
                                      issuerCurve,
                                      contractRecovery,
                                      numStepsPerYear,
                                      prot_method)

        fwdDf = 1.0

        if self._longProtection:
            longProt = +1
        else:
            longProt = -1

        fullPV = fwdDf * longProt * \
            (protPV - self._runningCoupon * fullRPV01 * self._notional)
        cleanPV = fwdDf * longProt * \
            (protPV - self._runningCoupon * cleanRPV01 * self._notional)

#        print("protLeg", protPV, "cleanRPV01", cleanRPV01, "value", cleanPV)

        return {'full_pv': fullPV, 'clean_pv': cleanPV}

###############################################################################

    def creditDV01(self,
                   valuationDate,
                   issuerCurve,
                   contractRecovery=standardRecovery,
                   pv01Method=0,
                   prot_method=0,
                   numStepsPerYear=25):
        ''' Calculation of the change in the value of the CDS contract for a
        one basis point change in the level of the CDS curve.'''

        v0 = self.value(valuationDate,
                        issuerCurve,
                        contractRecovery,
                        pv01Method,
                        prot_method,
                        numStepsPerYear)

        bump = 0.0001  # 1 basis point

        # we create a deep copy to avoid state issues
        bumpedIssuerCurve = deepcopy(issuerCurve)
        for cds in bumpedIssuerCurve._cdsContracts:
            cds._runningCoupon += bump

        bumpedIssuerCurve._buildCurve()

        v1 = self.value(valuationDate,
                        bumpedIssuerCurve,
                        contractRecovery,
                        pv01Method,
                        prot_method,
                        numStepsPerYear)

        creditDV01 = (v1['full_pv'] - v0['full_pv'])
        return creditDV01

###############################################################################

    def interestDV01(self,
                     valuationDate: FinDate,
                     issuerCurve,
                     contractRecovery=standardRecovery,
                     pv01Method: int = 0,
                     prot_method: int = 0,
                     numStepsPerYear: int = 25):
        ''' Calculation of the interest DV01 based on a simple bump of
        the discount factors and reconstruction of the CDS curve. '''

        v0 = self.value(valuationDate,
                        issuerCurve,
                        contractRecovery,
                        pv01Method,
                        prot_method,
                        numStepsPerYear)

        # we create a deep copy to avoid state issues
        newIssuerCurve = deepcopy(issuerCurve)

        bump = 0.0001  # 1 basis point

        for depo in newIssuerCurve._liborCurve._usedDeposits:
            depo._depositRate += bump
        for fra in newIssuerCurve._liborCurve._usedFRAs:
            fra._fraRate += bump
        for swap in newIssuerCurve._liborCurve._usedSwaps:
            swap._fixedLeg._coupon += bump

        newIssuerCurve._liborCurve._buildCurve()

        newIssuerCurve._buildCurve()

        v1 = self.value(valuationDate,
                        newIssuerCurve,
                        contractRecovery,
                        pv01Method,
                        prot_method,
                        numStepsPerYear)

        interestDV01 = (v1['full_pv'] - v0['full_pv'])
        return interestDV01

###############################################################################

    def cashSettlementAmount(self,
                             valuationDate,
                             settlementDate,
                             issuerCurve,
                             contractRecovery=standardRecovery,
                             pv01Method=0,
                             prot_method=0,
                             numStepsPerYear=25):
        ''' Value of the contract on the settlement date including accrued
        interest. '''

        v = self.value(valuationDate,
                       issuerCurve,
                       contractRecovery,
                       pv01Method,
                       prot_method,
                       numStepsPerYear)

        liborCurve = issuerCurve._liborCurve
        df = liborCurve.df(settlementDate)
        v = v / df
        return v

###############################################################################

    def cleanPrice(self,
                   valuationDate,
                   issuerCurve,
                   contractRecovery=standardRecovery,
                   pv01Method=0,
                   prot_method=0,
                   numStepsPerYear=52):
        ''' Value of the CDS contract excluding accrued interest. '''

        riskyPV01 = self.riskyPV01(valuationDate, issuerCurve, pv01Method)

        cleanRPV01 = riskyPV01['clean_rpv01']

        protPV = self.protectionLegPV(valuationDate,
                                      issuerCurve,
                                      contractRecovery,
                                      numStepsPerYear,
                                      prot_method)

        fwdDf = 1.0

        cleanPV = fwdDf * (protPV - self._runningCoupon * cleanRPV01
                           * self._notional)
        cleanPrice = (self._notional - cleanPV) / self._notional * 100.0
        return cleanPrice

###############################################################################

    def riskyPV01OLD(self,
                      valuationDate,
                      issuerCurve,
                      pv01Method=0):
        ''' RiskyPV01 of the contract using the OLD method. '''

        paymentDates = self._adjustedDates
        dayCount = FinDayCount(self._dayCountType)

        couponAccruedIndicator = 1

        # Method 0 : This is the market standard which assumes that the coupon
        # accrued is treated as though on average default occurs roughly midway
        # through a coupon period.

        teff = self._stepInDate
        pcd = paymentDates[0]  # PCD
        ncd = paymentDates[1]  # NCD

        # The first coupon is a special case which must be handled carefully
        # taking into account what coupon has already accrued and what has not
        qeff = issuerCurve.survivalProbability(teff)
        q1 = issuerCurve.survivalProbability(ncd)
        z1 = issuerCurve.df(ncd)

        # this is the part of the coupon accrued from the previous coupon date
        # to now
        accrualFactorPCDToNow = dayCount.yearFrac(pcd, teff)[0]

        # full first coupon is paid at the end of the current period if the
        yearFrac = dayCount.yearFrac(pcd, ncd)[0]

        # reference credit survives to the premium payment date
        fullRPV01 = q1 * z1 * yearFrac

        # coupon accrued from previous coupon to today paid in full at default
        # before coupon payment
        fullRPV01 = fullRPV01 + z1 * \
            (qeff - q1) * accrualFactorPCDToNow * couponAccruedIndicator

        # future accrued from now to coupon payment date assuming default
        # roughly midway
        fullRPV01 = fullRPV01 + 0.5 * z1 * \
            (qeff - q1) * (yearFrac - accrualFactorPCDToNow) \
            * couponAccruedIndicator

        for it in range(2, len(paymentDates)):

            t1 = paymentDates[it - 1]
            t2 = paymentDates[it]
            q2 = issuerCurve.survivalProbability(t2)
            z2 = issuerCurve.df(t2)

            accrualFactor = dayCount.yearFrac(t1, t2)[0]

            # full coupon is paid at the end of the current period if survives
            # to payment date
            fullRPV01 += q2 * z2 * accrualFactor

            ###################################################################

            if couponAccruedIndicator == 1:

                if useFlatHazardRateIntegral:
                    # This needs to be updated to handle small h+r
                    tau = accrualFactor
                    h12 = -log(q2 / q1) / tau
                    r12 = -log(z2 / z1) / tau
                    alpha = h12 + r12
                    expTerm = 1.0 - exp(-alpha * tau) - \
                        alpha * tau * exp(-alpha * tau)
                    dfullRPV01 = q1 * z1 * h12 * \
                        expTerm / abs(alpha * alpha + 1e-20)
                else:
                    dfullRPV01 = 0.50 * (q1 - q2) * z2 * accrualFactor

                fullRPV01 = fullRPV01 + dfullRPV01

            ###################################################################

            q1 = q2

        cleanRPV01 = fullRPV01 - accrualFactorPCDToNow

#        print("OLD PV01",fullRPV01, cleanRPV01)

        return {'full_rpv01': fullRPV01, 'clean_rpv01': cleanRPV01}

###############################################################################

    def accruedDays(self):
        ''' Number of days between the previous coupon and the currrent step
        in date. '''

        # I assume accrued runs to the effective date
        paymentDates = self._adjustedDates
        pcd = paymentDates[0]
        accruedDays = (self._stepInDate - pcd)
        return accruedDays

###############################################################################

    def accruedInterest(self):
        ''' Calculate the amount of accrued interest that has accrued from the
        previous coupon date (PCD) to the stepInDate of the CDS contract. '''

        dayCount = FinDayCount(self._dayCountType)
        paymentDates = self._adjustedDates
        pcd = paymentDates[0]
        accrualFactor = dayCount.yearFrac(pcd, self._stepInDate)[0]
        accruedInterest = accrualFactor * self._notional * self._runningCoupon

        if self._longProtection:
            accruedInterest *= -1.0

        return accruedInterest

###############################################################################

    def protectionLegPV(self,
                        valuationDate,
                        issuerCurve,
                        contractRecovery=standardRecovery,
                        numStepsPerYear=25,
                        protMethod=0):
        ''' Calculates the protection leg PV of the CDS by calling into the
        fast NUMBA code that has been defined above. '''

        teff = (self._stepInDate - valuationDate) / gDaysInYear
        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        liborCurve = issuerCurve._liborCurve

        v = _protectionLegPV_NUMBA(teff,
                                   tmat,
                                   liborCurve._times,
                                   liborCurve._dfs,
                                   issuerCurve._times,
                                   issuerCurve._values,
                                   contractRecovery,
                                   numStepsPerYear,
                                   protMethod)

        return v * self._notional

###############################################################################

    def riskyPV01(self,
                  valuationDate,
                  issuerCurve,
                  pv01Method=0):
        ''' The riskyPV01 is the present value of a risky one dollar paid on
        the premium leg of a CDS contract. '''

        liborCurve = issuerCurve._liborCurve

        paymentTimes = []
        for it in range(0, len(self._adjustedDates)):
            t = (self._adjustedDates[it] - valuationDate) / gDaysInYear
            paymentTimes.append(t)

        # this is the part of the coupon accrued from the previous coupon date
        # to now
        pcd = self._adjustedDates[0]
        eff = self._stepInDate
        dayCount = FinDayCount(self._dayCountType)

        accrualFactorPCDToNow = dayCount.yearFrac(pcd, eff)[0]

        yearFracs = self._accrualFactors
        teff = (eff - valuationDate) / gDaysInYear

        valueRPV01 = _riskyPV01_NUMBA(teff,
                                      accrualFactorPCDToNow,
                                      np.array(paymentTimes),
                                      np.array(yearFracs),
                                      liborCurve._times,
                                      liborCurve._dfs,
                                      issuerCurve._times,
                                      issuerCurve._values,
                                      pv01Method)

        fullRPV01 = valueRPV01[0]
        cleanRPV01 = valueRPV01[1]

#        print("NEW PV01",fullRPV01, cleanRPV01)
        return {'full_rpv01': fullRPV01, 'clean_rpv01': cleanRPV01}

###############################################################################

    def premiumLegPV(self,
                     valuationDate,
                     issuerCurve,
                     pv01Method=0):
        ''' Value of the premium leg of a CDS. '''

        fullRPV01 = self.riskyPV01(valuationDate,
                                   issuerCurve,
                                   pv01Method)['full_rpv01']

        v = fullRPV01 * self._notional * self._runningCoupon
        return v

###############################################################################

    def parSpread(self,
                  valuationDate,
                  issuerCurve,
                  contractRecovery=standardRecovery,
                  numStepsPerYear=25,
                  pv01Method=0,
                  protMethod=0):
        ''' Breakeven CDS coupon that would make the value of the CDS contract
        equal to zero. '''

        cleanRPV01 = self.riskyPV01(valuationDate,
                                    issuerCurve,
                                    pv01Method)['clean_rpv01']

        prot = self.protectionLegPV(valuationDate,
                                    issuerCurve,
                                    contractRecovery,
                                    numStepsPerYear,
                                    protMethod)

        # By convention this is calculated using the clean RPV01
        spd = prot / cleanRPV01 / self._notional
        return spd

###############################################################################

    def valueFastApprox(self,
                        valuationDate,
                        flatContinuousInterestRate,
                        flatCDSCurveSpread,
                        curveRecovery=standardRecovery,
                        contractRecovery=standardRecovery):
        ''' Implementation of fast valuation of the CDS contract using an
        accurate approximation that avoids curve building. '''

        if type(valuationDate) is not FinDate:
            raise FinError("Valuation date must be a FinDate and not " +
                           str(valuationDate))

        t_mat = (self._maturityDate - valuationDate) / gDaysInYear
        t_eff = (self._stepInDate - valuationDate) / gDaysInYear

        h = flatCDSCurveSpread / (1.0 - curveRecovery)
        r = flatContinuousInterestRate
        fwdDf = 1.0
        bumpSize = 0.0001

        if self._longProtection:
            longProtection = +1
        else:
            longProtection = -1

        # The sign of he accrued has already been sign adjusted for direction
        accrued = self.accruedInterest()

        # This is the clean RPV01 as it treats the PV01 stream as though it
        # pays just the accrued for the time between 0 and the maturity
        # It therefore omits the part that has accrued

        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        protPV = h * (1.0 - contractRecovery) * (z / w) * self._notional
        cleanPV = fwdDf * longProtection * \
            (protPV - self._runningCoupon * cleanRPV01 * self._notional)
        fullPV = cleanPV + fwdDf * accrued

        #######################################################################
        # bump CDS spread and calculate
        #######################################################################

        h = (flatCDSCurveSpread + bumpSize) / (1.0 - contractRecovery)
        r = flatContinuousInterestRate
        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        protPV = h * (1.0 - contractRecovery) * (z / w) * self._notional
        cleanPV_credit_bumped = fwdDf * longProtection * \
            (protPV - self._runningCoupon * cleanRPV01 * self._notional)
        fullPV_credit_bumped = cleanPV_credit_bumped \
            + fwdDf * longProtection * accrued
        credit01 = fullPV_credit_bumped - fullPV

        #######################################################################
        # bump Rate and calculate
        #######################################################################

        h = flatCDSCurveSpread / (1.0 - contractRecovery)
        r = flatContinuousInterestRate + bumpSize

        w = r + h
        z = np.exp(-w * t_eff) - np.exp(-w * t_mat)
        cleanRPV01 = (z / w) * 365.0 / 360.0
        protPV = h * (1.0 - contractRecovery) * (z/w) * self._notional
        cleanPV_ir_bumped = fwdDf * longProtection * \
            (protPV - self._runningCoupon * cleanRPV01 * self._notional)
        fullPV_ir_bumped = cleanPV_ir_bumped + fwdDf * longProtection * accrued
        ir01 = fullPV_ir_bumped - fullPV

        return (fullPV, cleanPV, credit01, ir01)

###############################################################################

    def printFlows(self, issuerCurve):

        numFlows = len(self._adjustedDates)

        print("PAYMENT_DATE      YEAR_FRAC      FLOW           DF       SURV_PROB      NPV")

        for it in range(1, numFlows):
            dt = self._adjustedDates[it]
            accFactor = self._accrualFactors[it]
            flow = self._flows[it]
            z = issuerCurve.df(dt)
            q = issuerCurve.survProb(dt)
            print("%15s %10.6f %12.2f %12.6f %12.6f %12.2f" %
                  (dt, accFactor, flow, z, q, flow * z * q))

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
        s += labelToString("ACCRUED DAYS", self.accruedDays())

        header = "PAYMENT_DATE, YEAR_FRAC, FLOW"
        valueTable = [self._adjustedDates, self._accrualFactors, self._flows]
        precision = "12.6f"

        s += tableToString(header, valueTable, precision)

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
