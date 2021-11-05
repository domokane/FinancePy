##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt


from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION, INVROOT2PI, N
from ...utils.error import FinError
from ...products.credit.cds_curve import CDSCurve
from ...products.credit.cds import CDS
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import label_to_string

RPV01_INDEX = 1  # 0 is FULL, 1 is CLEAN

###############################################################################


class CDSIndexOption:

    """ Class to manage the pricing and risk management of an option to enter
    into a CDS index. Different pricing algorithms are presented."""

    def __init__(self,
                 expiry_date: Date,
                 maturity_date: Date,
                 index_coupon: float,
                 strike_coupon: float,
                 notional: float = ONE_MILLION,
                 long_protection: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Initialisation of the class object. Note that a large number of the
        inputs are set to default values in line with the standard contract."""

        check_argument_types(self.__init__, locals())

        if expiry_date > maturity_date:
            raise FinError("Expiry date after end date")

        if index_coupon < 0.0:
            raise FinError("Index coupon is negative")

        if strike_coupon < 0.0:
            raise FinError("Index Option strike coupon is negative")

        self._expiry_date = expiry_date
        self._maturity_date = maturity_date
        self._index_coupon = index_coupon
        self._strike_coupon = strike_coupon
        self._notional = notional
        self._long_protection = long_protection

        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

        self._cds_contract = CDS(self._expiry_date,
                                 self._maturity_date,
                                 self._index_coupon,
                                 1.0,
                                 self._long_protection,
                                 self._freq_type,
                                 self._day_count_type,
                                 self._calendar_type,
                                 self._bus_day_adjust_type,
                                 self._date_gen_rule_type)

###############################################################################

    def value_adjusted_black(self,
                             valuation_date,
                             index_curve,
                             indexRecovery,
                             libor_curve,
                             sigma):
        """ This approach uses two adjustments to Black's option pricing
        model to value an option on a CDS index. """

        k = self._strike_coupon
        c = self._index_coupon
        time_to_expiry = (self._expiry_date - valuation_date) / gDaysInYear
        df = libor_curve.df(self._expiry_date)
        qExpiryIndex = index_curve.survival_prob(time_to_expiry)

        cds = CDS(valuation_date, self._maturity_date, k)
        strikeCurve = CDSCurve(
            valuation_date, [cds], libor_curve, indexRecovery)
#        qExpiryStrike = strikeCurve.survivalProbability(time_to_expiry)

        strikeRPV01 = self._cds_contract.risky_pv01(
            valuation_date, strikeCurve)['clean_rpv01']
        indexRPV01 = self._cds_contract.risky_pv01(
            valuation_date, index_curve)['clean_rpv01']

        s = self._cds_contract.par_spread(valuation_date, index_curve)

        fep = df * (1.0 - qExpiryIndex) * (1.0 - indexRecovery)
        adjFwd = s + fep / indexRPV01
        adjStrike = c + (k - c) * strikeRPV01 / indexRPV01 / qExpiryIndex

        denom = sigma * sqrt(time_to_expiry)
        d1 = log(adjFwd / adjStrike) + 0.5 * sigma * sigma * time_to_expiry
        d2 = log(adjFwd / adjStrike) - 0.5 * sigma * sigma * time_to_expiry
        d1 /= denom
        d2 /= denom

        v_pay = (adjFwd * N(d1) - adjStrike * N(d2)) * indexRPV01
        v_rec = (adjStrike * N(-d2) - adjFwd * N(-d1)) * indexRPV01

        v_pay *= self._notional
        v_rec *= self._notional

        return (v_pay, v_rec)

###############################################################################

    def value_anderson(self,
                       valuation_date,
                       issuer_curves,
                       indexRecovery,
                       sigma):
        """ This function values a CDS index option following approach by
        Anderson (2006). This ensures that a no-arbitrage relationship between
        the constituent CDS contract and the CDS index is enforced. It models
        the forward spread as a log-normally distributed quantity and uses the
        credit triangle to compute the forward RPV01. """

        num_credits = len(issuer_curves)
        time_to_expiry = (self._expiry_date - valuation_date) / gDaysInYear
#        timeToMaturity = (self._maturity_date - valuation_date) / gDaysInYear
        dfToExpiry = issuer_curves[0].df(time_to_expiry)
        libor_curve = issuer_curves[0]._libor_curve

        k = self._strike_coupon
        c = self._index_coupon

        strikeCDS = CDS(
            self._expiry_date,
            self._maturity_date,
            self._strike_coupon,
            1.0)
        strikeCurve = CDSCurve(valuation_date, [strikeCDS], libor_curve)
        strikeRPV01s = strikeCDS.risky_pv01(valuation_date, strikeCurve)
        qToExpiry = strikeCurve.survival_prob(time_to_expiry)
        strikeValue = (k - c) * strikeRPV01s['clean_rpv01']
        strikeValue /= (dfToExpiry * qToExpiry)

        expH = 0.0
        h1 = 0.0
        h2 = 0.0

        for iCredit in range(0, num_credits):

            issuer_curve = issuer_curves[iCredit]
            q = issuer_curve.survival_prob(time_to_expiry)
            dh1 = (1.0 - issuer_curve._recovery_rate) * (1.0 - q)

            s = self._cds_contract.par_spread(valuation_date, issuer_curve)
            rpv01 = self._cds_contract.risky_pv01(valuation_date, issuer_curve)
            dh2 = (s - c) * rpv01['clean_rpv01'] / (dfToExpiry * qToExpiry)

            h1 = h1 + dh1
            h2 = h2 + dh2

        expH = (h1 + h2) / num_credits

        x = self._solve_for_x(valuation_date,
                              sigma,
                              c,
                              indexRecovery,
                              libor_curve,
                              expH)

        v = self._calc_index_payer_option_price(valuation_date,
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

    def _solve_for_x(self,
                     valuation_date,
                     sigma,
                     index_coupon,
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

        f = self._calc_obj_func(x1, valuation_date, sigma, index_coupon,
                                indexRecovery, libor_curve) - expH

        fmid = self._calc_obj_func(x2, valuation_date, sigma, index_coupon,
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
            fmid = self._calc_obj_func(xmid, valuation_date, sigma,
                                       index_coupon,
                                       indexRecovery, libor_curve) - expH
            if fmid <= 0.0:
                rtb = xmid
            if abs(dx) < xacc or abs(fmid) < ftol:
                return rtb

        return rtb

###############################################################################

    def _calc_obj_func(self,
                       x,
                       valuation_date,
                       sigma,
                       index_coupon,
                       indexRecovery,
                       libor_curve):
        """ An internal function used in the Anderson valuation. """

        # The strike value is not relevant here as we want the zeroth element
        # of the return value
        strikeValue = 0.0

        values = self._calc_index_payer_option_price(valuation_date,
                                                     x,
                                                     sigma,
                                                     self._index_coupon,
                                                     strikeValue,
                                                     libor_curve,
                                                     indexRecovery)

        return values[0]

###############################################################################

    def _calc_index_payer_option_price(self,
                                       valuation_date,
                                       x,
                                       sigma,
                                       index_coupon,
                                       strikeValue,
                                       libor_curve,
                                       indexRecovery):
        """ Calculates the intrinsic value of the index payer swap and the
        value of the index payer option which are both returned in an array.
        """

        z = -6.0
        dz = 0.2
        numZSteps = int(2.0 * abs(z) / dz)

        flow_dates = self._cds_contract._adjusted_dates
        num_flows = len(flow_dates)
        texp = (self._expiry_date - valuation_date) / gDaysInYear
        dfToExpiry = libor_curve.df(self._expiry_date)
        lgd = 1.0 - indexRecovery

        fwdDfs = [1.0] * (num_flows)
        expiryToFlowTimes = [1.0] * (num_flows)

        for iFlow in range(1, num_flows):
            expiryToFlowTimes[iFlow] = (
                flow_dates[iFlow] - self._expiry_date) / gDaysInYear
            fwdDfs[iFlow] = libor_curve.df(flow_dates[iFlow]) / dfToExpiry

        intH = 0.0
        intMaxH = 0.0

        day_count = DayCount(self._day_count_type)
        pcd = flow_dates[0]  # PCD
        eff = self._expiry_date
        accrual_factorPCDToExpiry = day_count.year_frac(pcd, eff)[0]

        s0 = exp(-0.5 * sigma * sigma * texp)

        for _ in range(0, numZSteps):
            s = x * s0 * exp(sigma * sqrt(texp) * z)
            pdf = exp(-(z**2) / 2.0)
            z = z + dz

            fwdRPV01 = 0.0
            for iFlow in range(1, num_flows):
                acc_factor = self._cds_contract._accrual_factors[iFlow]
                survivalProbability = exp(-s * expiryToFlowTimes[iFlow] / lgd)
                fwdRPV01 += acc_factor * survivalProbability * fwdDfs[iFlow]

            fwdRPV01 += -accrual_factorPCDToExpiry
            h = (s - index_coupon) * fwdRPV01
            maxh = max(h - strikeValue, 0.0)

            intH += h * pdf
            intMaxH += maxh * pdf

        intH *= INVROOT2PI * dz
        intMaxH *= INVROOT2PI * dz * dfToExpiry
        return intH, intMaxH

###############################################################################

    def __repr__(self):
        """ print out details of the CDS contract and all of the calculated
        cash flows """
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("INDEX COUPON", self._index_coupon*10000, "bp\n")
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("LONG PROTECTION", self._long_protection)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAYCOUNT", self._day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUSDAYRULE", self._bus_day_adjust_type)
        s += label_to_string("DATEGENRULE", self._date_gen_rule_type)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
