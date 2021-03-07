##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: There are several speed ups for the Monte-Carlo including calculating
# all default baskets at the same time.

import numpy as np

from ...utils.FinError import FinError

from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes

from ...products.credit.cds import FinCDS

from ...models.credit_gaussian_copula_onefactor import homogeneousBasketLossDbn
from ...models.credit_gaussian_copula import default_timesGC
from ...models.credit_student_t_copula import FinModelStudentTCopula

from ...products.credit.cds_curve import FinCDSCurve

from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...market.discount.interpolator import interpolate, FinInterpTypes

from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import labelToString

###############################################################################
# TODO: Convert functions to use NUMBA!!
###############################################################################


class FinCDSBasket(object):

    """ Class to deal with n-to-default CDS baskets. """

    def __init__(self,
                 step_in_date: Date,
                 maturity_date: Date,
                 notional: float = ONE_MILLION,
                 running_coupon: float = 0.0,
                 long_protection: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):

        check_argument_types(self.__init__, locals())

        self._step_in_date = step_in_date
        self._maturity_date = maturity_date
        self._notional = notional
        self._running_coupon = running_coupon / 10000.0
        self._long_protection = long_protection
        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

        self._cds_contract = FinCDS(self._step_in_date,
                                   self._maturity_date,
                                   self._running_coupon,
                                   1.0,
                                   self._long_protection,
                                   self._freq_type,
                                   self._day_count_type,
                                   self._calendar_type,
                                   self._bus_day_adjust_type,
                                   self._date_gen_rule_type)

###############################################################################

    def valueLegs_MC(self,
                     valuation_date,
                     nToDefault,
                     default_times,
                     issuer_curves,
                     libor_curve):
        """ Value the legs of the default basket using Monte Carlo. The default
        times are an input so this valuation is not model dependent. """

        num_credits = default_times.shape[0]
        num_trials = default_times.shape[1]

        adjusted_dates = self._cds_contract._adjusted_dates
        num_flows = len(adjusted_dates)
        dayCount = DayCount(self._day_count_type)

        averageAccrualFactor = 0.0

        rpv01ToTimes = np.zeros(num_flows)
        for iTime in range(1, num_flows):
            t = (adjusted_dates[iTime] - valuation_date) / gDaysInYear
            dt0 = adjusted_dates[iTime - 1]
            dt1 = adjusted_dates[iTime]
            accrual_factor = dayCount.year_frac(dt0, dt1)[0]
            averageAccrualFactor += accrual_factor
            rpv01ToTimes[iTime] = rpv01ToTimes[iTime - 1] + \
                accrual_factor * libor_curve._df(t)

        averageAccrualFactor /= num_flows

        tmat = (self._maturity_date - valuation_date) / gDaysInYear

        rpv01 = 0.0
        prot = 0.0

        assetTau = np.zeros(num_credits)

        for iTrial in range(0, num_trials):
            for iCredit in range(0, num_credits):
                assetTau[iCredit] = default_times[iCredit, iTrial]

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
                for iCredit in range(0, num_credits):
                    if minTau == default_times[iCredit, iTrial]:
                        assetIndex = iCredit
                        break

                protTrial = (1.0 - issuer_curves[assetIndex]._recovery_rate)
                protTrial *= libor_curve._df(minTau)

            else:

                numPaymentsIndex = int(tmat / averageAccrualFactor)
                rpv01Trial = rpv01ToTimes[numPaymentsIndex]
                protTrial = 0.0

            rpv01 += rpv01Trial
            prot += protTrial

        rpv01 = rpv01 / num_trials
        prot = prot / num_trials
        return (rpv01, prot)

###############################################################################

    def valueGaussian_MC(self,
                         valuation_date,
                         nToDefault,
                         issuer_curves,
                         correlationMatrix,
                         libor_curve,
                         num_trials,
                         seed):
        """ Value the default basket using a Gaussian copula model. This
        depends on the issuer discount and correlation matrix. """

        num_credits = len(issuer_curves)

        if nToDefault > num_credits or nToDefault < 1:
            raise FinError("nToDefault must be 1 to num_credits")

        default_times = default_timesGC(issuer_curves,
                                      correlationMatrix,
                                      num_trials,
                                      seed)

        rpv01, prot_pv = self.valueLegs_MC(valuation_date,
                                          nToDefault,
                                          default_times,
                                          issuer_curves,
                                          libor_curve)

        spd = prot_pv / rpv01
        value = self._notional * (prot_pv - self._running_coupon * rpv01)

        if not self._long_protection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def valueStudentT_MC(self,
                         valuation_date,
                         nToDefault,
                         issuer_curves,
                         correlationMatrix,
                         degreesOfFreedom,
                         libor_curve,
                         num_trials,
                         seed):
        """ Value the default basket using the Student-T copula. """

        num_credits = len(issuer_curves)

        if nToDefault > num_credits or nToDefault < 1:
            raise FinError("nToDefault must be 1 to num_credits")

        model = FinModelStudentTCopula()

        default_times = model.default_times(issuer_curves,
                                          correlationMatrix,
                                          degreesOfFreedom,
                                          num_trials,
                                          seed)

        rpv01, prot_pv = self.valueLegs_MC(valuation_date,
                                          nToDefault,
                                          default_times,
                                          issuer_curves,
                                          libor_curve)

        spd = prot_pv / rpv01
        value = self._notional * (prot_pv - self._running_coupon * rpv01)

        if not self._long_protection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def value1FGaussian_Homo(self,
                             valuation_date,
                             nToDefault,
                             issuer_curves,
                             betaVector,
                             libor_curve,
                             num_points=50):
        """ Value default basket using 1 factor Gaussian copula and analytical
        approach which is only exact when all recovery rates are the same. """

        num_credits = len(issuer_curves)

        if num_credits == 0:
            raise FinError("Num Credits is zero")

        if nToDefault < 1 or nToDefault > num_credits:
            raise FinError("NToDefault must be 1 to num_credits")

        tmat = (self._maturity_date - valuation_date) / gDaysInYear

        if tmat < 0.0:
            raise FinError("Value date is after maturity date")

        payment_dates = self._cds_contract._adjusted_dates
        numTimes = len(payment_dates)

        issuerSurvivalProbabilities = np.zeros(num_credits)
        recovery_rates = np.zeros(num_credits)
        basketTimes = np.zeros(numTimes)
        basketSurvivalCurve = np.zeros(numTimes)

        basketTimes[0] = 0.0
        basketSurvivalCurve[0] = 1.0

        for iTime in range(1, numTimes):

            t = (payment_dates[iTime] - valuation_date) / gDaysInYear

            for iCredit in range(0, num_credits):
                issuer_curve = issuer_curves[iCredit]
                recovery_rates[iCredit] = issuer_curve._recovery_rate
                issuerSurvivalProbabilities[iCredit] = interpolate(
                    t, issuer_curve._times, issuer_curve._values,
                    FinInterpTypes.FLAT_FWD_RATES.value)

            lossDbn = homogeneousBasketLossDbn(issuerSurvivalProbabilities,
                                               recovery_rates,
                                               betaVector,
                                               num_points)

            basketSurvivalCurve[iTime] = 1.0
            for iToDefault in range(nToDefault, num_credits + 1):
                basketSurvivalCurve[iTime] -= lossDbn[iToDefault]

            basketTimes[iTime] = t

        curveRecovery = recovery_rates[0]
        libor_curve = issuer_curves[0]._libor_curve
        basketCurve = FinCDSCurve(valuation_date, [], libor_curve, curveRecovery)
        basketCurve._times = basketTimes
        basketCurve._values = basketSurvivalCurve

        protLegPV = self._cds_contract.protectionLegPV(
            valuation_date, basketCurve, curveRecovery)
        riskyPV01 = self._cds_contract.riskyPV01(valuation_date, basketCurve)['clean_rpv01']

        # Long protection
        mtm = self._notional * (protLegPV - riskyPV01 * self._running_coupon)

        if not self._long_protection:
            mtm *= -1.0

        basketOutput = np.zeros(4)
        basketOutput[0] = mtm
        basketOutput[1] = riskyPV01 * self._notional * self._running_coupon
        basketOutput[2] = protLegPV * self._notional
        basketOutput[3] = protLegPV / riskyPV01

        return basketOutput

###############################################################################

    def __repr__(self):
        """ print out details of the CDS contract and all of the calculated
        cash flows """
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("STEP-IN DATE", self._step_in_date)
        s += labelToString("MATURITY", self._maturity_date)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("RUNNING COUPON", self._running_coupon*10000, "bp\n")
        s += labelToString("DAYCOUNT", self._day_count_type)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("CALENDAR", self._calendar_type)
        s += labelToString("BUSDAYRULE", self._bus_day_adjust_type)
        s += labelToString("DATEGENRULE", self._date_gen_rule_type)

#        header = "PAYMENT_DATE, YEAR_FRAC, FLOW"
#        valueTable = [self._adjusted_dates, self._accrual_factors, self._flows]
#        precision = "12.6f"
#        s += tableToString(header, valueTable, precision)

        return s

###############################################################################
