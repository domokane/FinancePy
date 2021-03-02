##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt


from ...utils.Calendar import FinCalendarTypes
from ...utils.Calendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...utils.DayCount import DayCount, FinDayCountTypes
from ...utils.Frequency import FinFrequencyTypes
from ...utils.FinGlobalVariables import gDaysInYear
from ...utils.Math import ONE_MILLION, INVROOT2PI, N
from ...utils.FinError import FinError
from ...products.credit.FinCDSCurve import FinCDSCurve
from ...products.credit.FinCDS import FinCDS
from ...utils.FinHelperFunctions import checkArgumentTypes
from ...utils.Date import Date
from ...utils.FinHelperFunctions import labelToString

RPV01_INDEX = 1  # 0 is FULL, 1 is CLEAN

###############################################################################


class FinCDSIndexOption(object):

    """ Class to manage the pricing and risk management of an option to enter
    into a CDS index. Different pricing algorithms are presented."""

    def __init__(self,
                 expiry_date: Date,
                 maturity_date: Date,
                 indexCoupon: float,
                 strikeCoupon: float,
                 notional: float = ONE_MILLION,
                 long_protection: bool = True,
                 freq_type: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 day_count_type: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 bus_day_adjust_type: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        """ Initialisation of the class object. Note that a large number of the
        inputs are set to default values in line with the standard contract."""

        checkArgumentTypes(self.__init__, locals())

        if expiry_date > maturity_date:
            raise FinError("Expiry date after end date")

        if indexCoupon < 0.0:
            raise FinError("Index coupon is negative")

        if strikeCoupon < 0.0:
            raise FinError("Index Option strike coupon is negative")

        self._expiry_date = expiry_date
        self._maturity_date = maturity_date
        self._indexCoupon = indexCoupon
        self._strikeCoupon = strikeCoupon
        self._notional = notional
        self._long_protection = long_protection

        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

        self._cdsContract = FinCDS(self._expiry_date,
                                   self._maturity_date,
                                   self._indexCoupon,
                                   1.0,
                                   self._long_protection,
                                   self._freq_type,
                                   self._day_count_type,
                                   self._calendar_type,
                                   self._bus_day_adjust_type,
                                   self._date_gen_rule_type)

###############################################################################

    def valueAdjustedBlack(self,
                           valuation_date,
                           index_curve,
                           indexRecovery,
                           libor_curve,
                           sigma):
        """ This approach uses two adjustments to Black's option pricing
        model to value an option on a CDS index. """

        k = self._strikeCoupon
        c = self._indexCoupon
        timeToExpiry = (self._expiry_date - valuation_date) / gDaysInYear
        df = libor_curve.df(self._expiry_date)
        qExpiryIndex = index_curve.survProb(timeToExpiry)

        cds = FinCDS(valuation_date, self._maturity_date, k)
        strikeCurve = FinCDSCurve(
            valuation_date, [cds], libor_curve, indexRecovery)
#        qExpiryStrike = strikeCurve.survivalProbability(timeToExpiry)

        strikeRPV01 = self._cdsContract.riskyPV01(
            valuation_date, strikeCurve)['clean_rpv01']
        indexRPV01 = self._cdsContract.riskyPV01(
            valuation_date, index_curve)['clean_rpv01']

        s = self._cdsContract.parSpread(valuation_date, index_curve)

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
                      valuation_date,
                      issuer_curves,
                      indexRecovery,
                      sigma):
        """ This function values a CDS index option following approach by
        Anderson (2006). This ensures that a no-arbitrage relationship between
        the consituent CDS contract and the CDS index is enforced. It models
        the forward spread as a log-normally distributed quantity and uses the
        credit triangle to compute the forward RPV01. """

        numCredits = len(issuer_curves)
        timeToExpiry = (self._expiry_date - valuation_date) / gDaysInYear
#        timeToMaturity = (self._maturity_date - valuation_date) / gDaysInYear
        dfToExpiry = issuer_curves[0].df(timeToExpiry)
        libor_curve = issuer_curves[0]._libor_curve

        k = self._strikeCoupon
        c = self._indexCoupon

        strikeCDS = FinCDS(
            self._expiry_date,
            self._maturity_date,
            self._strikeCoupon,
            1.0)
        strikeCurve = FinCDSCurve(valuation_date, [strikeCDS], libor_curve)
        strikeRPV01s = strikeCDS.riskyPV01(valuation_date, strikeCurve)
        qToExpiry = strikeCurve.survProb(timeToExpiry)
        strikeValue = (k - c) * strikeRPV01s['clean_rpv01']
        strikeValue /= (dfToExpiry * qToExpiry)

        expH = 0.0
        h1 = 0.0
        h2 = 0.0

        for iCredit in range(0, numCredits):

            issuer_curve = issuer_curves[iCredit]
            q = issuer_curve.survProb(timeToExpiry)
            dh1 = (1.0 - issuer_curve._recovery_rate) * (1.0 - q)

            s = self._cdsContract.parSpread(valuation_date, issuer_curve)
            rpv01 = self._cdsContract.riskyPV01(valuation_date, issuer_curve)
            dh2 = (s - c) * rpv01['clean_rpv01'] / (dfToExpiry * qToExpiry)

            h1 = h1 + dh1
            h2 = h2 + dh2

        expH = (h1 + h2) / numCredits

        x = self._solveForX(valuation_date,
                            sigma,
                            c,
                            indexRecovery,
                            libor_curve,
                            expH)

        v = self._calcIndexPayerOptionPrice(valuation_date,
                                            x,
                                            sigma,
                                            c,
                                            strikeValue,
                                            libor_curve,
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
                   valuation_date,
                   sigma,
                   indexCoupon,
                   indexRecovery,
                   libor_curve,
                   expH):
        """ Function to solve for the arbitrage free """
        x1 = 0.0
        x2 = 0.9999
        ftol = 1e-8
        jmax = 40
        xacc = 0.000000001
        rtb = 999999

        f = self._calcObjFunc(x1, valuation_date, sigma, indexCoupon,
                              indexRecovery, libor_curve) - expH

        fmid = self._calcObjFunc(x2, valuation_date, sigma, indexCoupon,
                                 indexRecovery, libor_curve) - expH

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
            fmid = self._calcObjFunc(xmid, valuation_date, sigma, indexCoupon,
                                     indexRecovery, libor_curve) - expH
            if fmid <= 0.0:
                rtb = xmid
            if abs(dx) < xacc or abs(fmid) < ftol:
                return rtb

        return rtb

###############################################################################

    def _calcObjFunc(self,
                     x,
                     valuation_date,
                     sigma,
                     indexCoupon,
                     indexRecovery,
                     libor_curve):
        """ An internal function used in the Anderson valuation. """

        # The strike value is not relevant here as we want the zeroth element
        # of the return value
        strikeValue = 0.0

        values = self._calcIndexPayerOptionPrice(valuation_date,
                                                 x,
                                                 sigma,
                                                 self._indexCoupon,
                                                 strikeValue,
                                                 libor_curve,
                                                 indexRecovery)

        return values[0]

###############################################################################

    def _calcIndexPayerOptionPrice(self,
                                   valuation_date,
                                   x,
                                   sigma,
                                   indexCoupon,
                                   strikeValue,
                                   libor_curve,
                                   indexRecovery):
        """ Calculates the intrinsic value of the index payer swap and the
        value of the index payer option which are both returned in an array.
        """

        z = -6.0
        dz = 0.2
        numZSteps = int(2.0 * abs(z) / dz)

        flow_dates = self._cdsContract._adjustedDates
        numFlows = len(flow_dates)
        texp = (self._expiry_date - valuation_date) / gDaysInYear
        dfToExpiry = libor_curve.df(self._expiry_date)
        lgd = 1.0 - indexRecovery

        fwdDfs = [1.0] * (numFlows)
        expiryToFlowTimes = [1.0] * (numFlows)

        for iFlow in range(1, numFlows):
            expiryToFlowTimes[iFlow] = (
                flow_dates[iFlow] - self._expiry_date) / gDaysInYear
            fwdDfs[iFlow] = libor_curve.df(flow_dates[iFlow]) / dfToExpiry

        intH = 0.0
        intMaxH = 0.0

        dayCount = DayCount(self._day_count_type)
        pcd = flow_dates[0]  # PCD
        eff = self._expiry_date
        accrualFactorPCDToExpiry = dayCount.year_frac(pcd, eff)[0]

        s0 = exp(-0.5 * sigma * sigma * texp)

        for _ in range(0, numZSteps):
            s = x * s0 * exp(sigma * sqrt(texp) * z)
            pdf = exp(-(z**2) / 2.0)
            z = z + dz

            fwdRPV01 = 0.0
            for iFlow in range(1, numFlows):
                acc_factor = self._cdsContract._accrualFactors[iFlow]
                survivalProbability = exp(-s * expiryToFlowTimes[iFlow] / lgd)
                fwdRPV01 += acc_factor * survivalProbability * fwdDfs[iFlow]

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
        """ print out details of the CDS contract and all of the calculated
        cashflows """
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("INDEX COUPON", self._indexCoupon*10000, "bp\n")
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("LONG PROTECTION", self._long_protection)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("DAYCOUNT", self._day_count_type)
        s += labelToString("CALENDAR", self._calendar_type)
        s += labelToString("BUSDAYRULE", self._bus_day_adjust_type)
        s += labelToString("DATEGENRULE", self._date_gen_rule_type)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
