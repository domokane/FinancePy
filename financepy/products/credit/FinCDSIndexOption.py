##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt


from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION, INVROOT2PI, N
from ...finutils.FinError import FinError
from ...products.credit.FinCDSCurve import FinCDSCurve
from ...products.credit.FinCDS import FinCDS
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString

RPV01_INDEX = 1  # 0 is FULL, 1 is CLEAN

###############################################################################


class FinCDSIndexOption(object):

    ''' Class to manage the pricing and risk management of an option to enter
    into a CDS index. Different pricing algorithms are presented.'''

    def __init__(self,
                 expiryDate: FinDate,
                 maturityDate: FinDate,
                 indexCoupon: float,
                 strikeCoupon: float,
                 notional: float = ONE_MILLION,
                 longProtection: bool = True,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Initialisation of the class object. Note that a large number of the
        inputs are set to default values in line with the standard contract.'''

        checkArgumentTypes(self.__init__, locals())

        if expiryDate > maturityDate:
            raise FinError("Expiry date after end date")

        if indexCoupon < 0.0:
            raise FinError("Index coupon is negative")

        if strikeCoupon < 0.0:
            raise FinError("Index Option strike coupon is negative")

        self._expiryDate = expiryDate
        self._maturityDate = maturityDate
        self._indexCoupon = indexCoupon
        self._strikeCoupon = strikeCoupon
        self._notional = notional
        self._longProtection = longProtection

        self._dayCountType = dayCountType
        self._dateGenRuleType = dateGenRuleType
        self._calendarType = calendarType
        self._freqType = freqType
        self._busDayAdjustType = busDayAdjustType

        self._cdsContract = FinCDS(self._expiryDate,
                                   self._maturityDate,
                                   self._indexCoupon,
                                   1.0,
                                   self._longProtection,
                                   self._freqType,
                                   self._dayCountType,
                                   self._calendarType,
                                   self._busDayAdjustType,
                                   self._dateGenRuleType)

###############################################################################

    def valueAdjustedBlack(self,
                           valuationDate,
                           indexCurve,
                           indexRecovery,
                           liborCurve,
                           sigma):
        ''' This approach uses two adjustments to Black's option pricing
        model to value an option on a CDS index. '''

        k = self._strikeCoupon
        c = self._indexCoupon
        timeToExpiry = (self._expiryDate - valuationDate) / gDaysInYear
        df = liborCurve.df(self._expiryDate)
        qExpiryIndex = indexCurve.survProb(timeToExpiry)

        cds = FinCDS(valuationDate, self._maturityDate, k)
        strikeCurve = FinCDSCurve(
            valuationDate, [cds], liborCurve, indexRecovery)
#        qExpiryStrike = strikeCurve.survivalProbability(timeToExpiry)

        strikeRPV01 = self._cdsContract.riskyPV01(
            valuationDate, strikeCurve)['clean_rpv01']
        indexRPV01 = self._cdsContract.riskyPV01(
            valuationDate, indexCurve)['clean_rpv01']

        s = self._cdsContract.parSpread(valuationDate, indexCurve)

        fep = df * (1.0 - qExpiryIndex) * (1.0 - indexRecovery)
        adjFwd = s + fep / indexRPV01
        adjStrike = c + (k - c) * strikeRPV01 / indexRPV01 / qExpiryIndex

        denom = sigma * sqrt(timeToExpiry)
        d1 = log(adjFwd / adjStrike) + 0.5 * sigma * sigma * timeToExpiry
        d2 = log(adjFwd / adjStrike) - 0.5 * sigma * sigma * timeToExpiry
        d1 /= denom
        d2 /= denom

        v_pay = (adjFwd * N(d1) - adjStrike * N(d2)) * indexRPV01
        v_rec = (adjStrike * N(-d2) - adjFwd * N(-d1)) * indexRPV01

        v_pay *= self._notional
        v_rec *= self._notional

        return (v_pay, v_rec)

###############################################################################

    def valueAnderson(self,
                      valuationDate,
                      issuerCurves,
                      indexRecovery,
                      sigma):
        ''' This function values a CDS index option following approach by
        Anderson (2006). This ensures that a no-arbitrage relationship between
        the consituent CDS contract and the CDS index is enforced. It models
        the forward spread as a log-normally distributed quantity and uses the
        credit triangle to compute the forward RPV01. '''

        numCredits = len(issuerCurves)
        timeToExpiry = (self._expiryDate - valuationDate) / gDaysInYear
#        timeToMaturity = (self._maturityDate - valuationDate) / gDaysInYear
        dfToExpiry = issuerCurves[0].df(timeToExpiry)
        liborCurve = issuerCurves[0]._liborCurve

        k = self._strikeCoupon
        c = self._indexCoupon

        strikeCDS = FinCDS(
            self._expiryDate,
            self._maturityDate,
            self._strikeCoupon,
            1.0)
        strikeCurve = FinCDSCurve(valuationDate, [strikeCDS], liborCurve)
        strikeRPV01s = strikeCDS.riskyPV01(valuationDate, strikeCurve)
        qToExpiry = strikeCurve.survProb(timeToExpiry)
        strikeValue = (k - c) * strikeRPV01s['clean_rpv01']
        strikeValue /= (dfToExpiry * qToExpiry)

        expH = 0.0
        h1 = 0.0
        h2 = 0.0

        for iCredit in range(0, numCredits):

            issuerCurve = issuerCurves[iCredit]
            q = issuerCurve.survProb(timeToExpiry)
            dh1 = (1.0 - issuerCurve._recoveryRate) * (1.0 - q)

            s = self._cdsContract.parSpread(valuationDate, issuerCurve)
            rpv01 = self._cdsContract.riskyPV01(valuationDate, issuerCurve)
            dh2 = (s - c) * rpv01['clean_rpv01'] / (dfToExpiry * qToExpiry)

            h1 = h1 + dh1
            h2 = h2 + dh2

        expH = (h1 + h2) / numCredits

        x = self._solveForX(valuationDate,
                            sigma,
                            c,
                            indexRecovery,
                            liborCurve,
                            expH)

        v = self._calcIndexPayerOptionPrice(valuationDate,
                                            x,
                                            sigma,
                                            c,
                                            strikeValue,
                                            liborCurve,
                                            indexRecovery)

        v = v[1]
        v_pay = v * self._notional
        v_rec = v_pay + (strikeValue - expH) * dfToExpiry * self._notional
        strikeValue *= 10000.0
        x *= 10000.0
        expH *= 10000.0
        return v_pay, v_rec, strikeValue, x, expH

###############################################################################

    def _solveForX(self,
                   valuationDate,
                   sigma,
                   indexCoupon,
                   indexRecovery,
                   liborCurve,
                   expH):
        ''' Function to solve for the arbitrage free '''
        x1 = 0.0
        x2 = 0.9999
        ftol = 1e-8
        jmax = 40
        xacc = 0.000000001
        rtb = 999999

        f = self._calcObjFunc(x1, valuationDate, sigma, indexCoupon,
                              indexRecovery, liborCurve) - expH

        fmid = self._calcObjFunc(x2, valuationDate, sigma, indexCoupon,
                                 indexRecovery, liborCurve) - expH

        if f * fmid >= 0.0:
            raise FinError("Solution not bracketed.")

        if f < 0.0:
            rtb = x1
            dx = x2 - x1
        else:
            rtb = x2
            dx = x1 - x2

        for _ in range(0, jmax):
            dx = dx * 0.5
            xmid = rtb + dx
            fmid = self._calcObjFunc(xmid, valuationDate, sigma, indexCoupon,
                                     indexRecovery, liborCurve) - expH
            if fmid <= 0.0:
                rtb = xmid
            if abs(dx) < xacc or abs(fmid) < ftol:
                return rtb

        return rtb

###############################################################################

    def _calcObjFunc(self,
                     x,
                     valuationDate,
                     sigma,
                     indexCoupon,
                     indexRecovery,
                     liborCurve):
        ''' An internal function used in the Anderson valuation. '''

        # The strike value is not relevant here as we want the zeroth element
        # of the return value
        strikeValue = 0.0

        values = self._calcIndexPayerOptionPrice(valuationDate,
                                                 x,
                                                 sigma,
                                                 self._indexCoupon,
                                                 strikeValue,
                                                 liborCurve,
                                                 indexRecovery)

        return values[0]

###############################################################################

    def _calcIndexPayerOptionPrice(self,
                                   valuationDate,
                                   x,
                                   sigma,
                                   indexCoupon,
                                   strikeValue,
                                   liborCurve,
                                   indexRecovery):
        ''' Calculates the intrinsic value of the index payer swap and the
        value of the index payer option which are both returned in an array.
        '''

        z = -6.0
        dz = 0.2
        numZSteps = int(2.0 * abs(z) / dz)

        flowDates = self._cdsContract._adjustedDates
        numFlows = len(flowDates)
        texp = (self._expiryDate - valuationDate) / gDaysInYear
        dfToExpiry = liborCurve.df(self._expiryDate)
        lgd = 1.0 - indexRecovery

        fwdDfs = [1.0] * (numFlows)
        expiryToFlowTimes = [1.0] * (numFlows)

        for iFlow in range(1, numFlows):
            expiryToFlowTimes[iFlow] = (
                flowDates[iFlow] - self._expiryDate) / gDaysInYear
            fwdDfs[iFlow] = liborCurve.df(flowDates[iFlow]) / dfToExpiry

        intH = 0.0
        intMaxH = 0.0

        dayCount = FinDayCount(self._dayCountType)
        pcd = flowDates[0]  # PCD
        eff = self._expiryDate
        accrualFactorPCDToExpiry = dayCount.yearFrac(pcd, eff)[0]

        s0 = exp(-0.5 * sigma * sigma * texp)

        for _ in range(0, numZSteps):
            s = x * s0 * exp(sigma * sqrt(texp) * z)
            pdf = exp(-(z**2) / 2.0)
            z = z + dz

            fwdRPV01 = 0.0
            for iFlow in range(1, numFlows):
                accFactor = self._cdsContract._accrualFactors[iFlow]
                survivalProbability = exp(-s * expiryToFlowTimes[iFlow] / lgd)
                fwdRPV01 += accFactor * survivalProbability * fwdDfs[iFlow]

            fwdRPV01 += -accrualFactorPCDToExpiry
            h = (s - indexCoupon) * fwdRPV01
            maxh = max(h - strikeValue, 0.0)

            intH += h * pdf
            intMaxH += maxh * pdf

        intH *= INVROOT2PI * dz
        intMaxH *= INVROOT2PI * dz * dfToExpiry
        return intH, intMaxH

###############################################################################

    def __repr__(self):
        ''' print out details of the CDS contract and all of the calculated
        cashflows '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("INDEX COUPON", self._indexCoupon*10000, "bp\n")
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("LONG PROTECTION", self._longProtection)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("DAYCOUNT", self._dayCountType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUSDAYRULE", self._busDayAdjustType)
        s += labelToString("DATEGENRULE", self._dateGenRuleType)
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
