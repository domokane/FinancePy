##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: There are several speed ups for the Monte-Carlo including calculating
# all default baskets at the same time.

import numpy as np

from ...utils.error import FinError

from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes

from ...products.credit.cds import CDS

from ...models.gauss_copula_onefactor import homog_basket_loss_dbn
from ...models.gauss_copula import default_times_gc
from ...models.student_t_copula import StudentTCopula

from ...products.credit.cds_curve import CDSCurve

from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...market.curves.interpolator import interpolate, InterpTypes

from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import label_to_string

###############################################################################
# TODO: Convert functions to use NUMBA!!
###############################################################################


class CDSBasket:

    """ Class to deal with n-to-default CDS baskets. """

    def __init__(self,
                 step_in_date: Date,
                 maturity_date: Date,
                 notional: float = ONE_MILLION,
                 running_coupon: float = 0.0,
                 long_protection: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 dc_type: DayCountTypes = DayCountTypes.ACT_360,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):

        check_argument_types(self.__init__, locals())

        self._step_in_date = step_in_date
        self._maturity_date = maturity_date
        self._notional = notional
        self._running_coupon = running_coupon / 10000.0
        self._long_protection = long_protection
        self._dc_type = dc_type
        self._dg_type = dg_type
        self._cal_type = cal_type
        self._freq_type = freq_type
        self._bd_type = bd_type

        self._cds_contract = CDS(self._step_in_date,
                                 self._maturity_date,
                                 self._running_coupon,
                                 1.0,
                                 self._long_protection,
                                 self._freq_type,
                                 self._dc_type,
                                 self._cal_type,
                                 self._bd_type,
                                 self._dg_type)

###############################################################################

    def value_legs_mc(self,
                      value_date,
                      n_to_default,
                      default_times,
                      issuer_curves,
                      libor_curve):
        """ Value the legs of the default basket using Monte Carlo. The default
        times are an input so this valuation is not model dependent. """

        num_credits = default_times.shape[0]
        num_trials = default_times.shape[1]

        payment_dates = self._cds_contract._payment_dates
        num_payments = len(payment_dates)
        day_count = DayCount(self._dc_type)

        averageAccrualFactor = 0.0

        rpv01ToTimes = np.zeros(num_payments)

        for iTime in range(1, num_payments):

            t = (payment_dates[iTime] - value_date) / gDaysInYear
            dt0 = payment_dates[iTime - 1]
            dt1 = payment_dates[iTime]
            accrual_factor = day_count.year_frac(dt0, dt1)[0]
            averageAccrualFactor += accrual_factor
            rpv01ToTimes[iTime] = rpv01ToTimes[iTime - 1] + \
                accrual_factor * libor_curve._df(t)

        averageAccrualFactor /= num_payments

        tmat = (self._maturity_date - value_date) / gDaysInYear

        rpv01 = 0.0
        prot = 0.0

        assetTau = np.zeros(num_credits)

        for i_trial in range(0, num_trials):

            for i_credit in range(0, num_credits):

                assetTau[i_credit] = default_times[i_credit, i_trial]

            # ORDER THE DEFAULT TIMES
            assetTau.sort()

            # GET THE Nth DEFAULT TIME
            minTau = assetTau[n_to_default - 1]

            if minTau < tmat:

                numPaymentsIndex = int(minTau / averageAccrualFactor)
                rpv01Trial = rpv01ToTimes[numPaymentsIndex]
                rpv01Trial += (minTau - numPaymentsIndex *
                               averageAccrualFactor)

                # DETERMINE IDENTITY OF N-TO-DEFAULT CREDIT IF BASKET NOT HOMO
                assetIndex = 0
                for i_credit in range(0, num_credits):
                    if minTau == default_times[i_credit, i_trial]:
                        assetIndex = i_credit
                        break

                protTrial = (1.0 - issuer_curves[assetIndex]._recovery_rate)
                protTrial *= libor_curve._df(minTau)

            else:

                numPaymentsIndex = int(tmat / averageAccrualFactor)
                rpv01Trial = rpv01ToTimes[-1]
                protTrial = 0.0

            rpv01 += rpv01Trial
            prot += protTrial

        rpv01 = rpv01 / num_trials
        prot = prot / num_trials
        return (rpv01, prot)

###############################################################################

    def value_gaussian_mc(self,
                          value_date,
                          n_to_default,
                          issuer_curves,
                          correlation_matrix,
                          libor_curve,
                          num_trials,
                          seed):
        """ Value the default basket using a Gaussian copula model. This
        depends on the issuer discount and correlation matrix. """

        num_credits = len(issuer_curves)

        if n_to_default > num_credits or n_to_default < 1:
            raise FinError("n_to_default must be 1 to num_credits")

        default_times = default_times_gc(issuer_curves,
                                         correlation_matrix,
                                         num_trials,
                                         seed)

        rpv01, prot_pv = self.value_legs_mc(value_date,
                                            n_to_default,
                                            default_times,
                                            issuer_curves,
                                            libor_curve)

        spd = prot_pv / rpv01
        value = self._notional * (prot_pv - self._running_coupon * rpv01)

        if not self._long_protection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def value_student_t_mc(self,
                           value_date,
                           n_to_default,
                           issuer_curves,
                           correlation_matrix,
                           degreesOfFreedom,
                           libor_curve,
                           num_trials,
                           seed):
        """ Value the default basket using the Student-T copula. """

        num_credits = len(issuer_curves)

        if n_to_default > num_credits or n_to_default < 1:
            raise FinError("n_to_default must be 1 to num_credits")

        model = StudentTCopula()

        default_times = model.default_times(issuer_curves,
                                            correlation_matrix,
                                            degreesOfFreedom,
                                            num_trials,
                                            seed)

        rpv01, prot_pv = self.value_legs_mc(value_date,
                                            n_to_default,
                                            default_times,
                                            issuer_curves,
                                            libor_curve)

        spd = prot_pv / rpv01
        value = self._notional * (prot_pv - self._running_coupon * rpv01)

        if not self._long_protection:
            value = value * -1.0

        return (value, rpv01, spd)

###############################################################################

    def value_1f_gaussian_homo(self,
                               value_date,
                               n_to_default,
                               issuer_curves,
                               beta_vector,
                               libor_curve,
                               num_points=50):
        """ Value default basket using 1 factor Gaussian copula and analytical
        approach which is only exact when all recovery rates are the same. """

        num_credits = len(issuer_curves)

        if num_credits == 0:
            raise FinError("Num Credits is zero")

        if n_to_default < 1 or n_to_default > num_credits:
            raise FinError("n_to_default must be 1 to num_credits")

        tmat = (self._maturity_date - value_date) / gDaysInYear

        if tmat < 0.0:
            raise FinError("Value date is after maturity date")

        payment_dates = self._cds_contract._payment_dates
        num_times = len(payment_dates)

        issuerSurvivalProbabilities = np.zeros(num_credits)
        recovery_rates = np.zeros(num_credits)
        basketTimes = np.zeros(num_times)
        basketSurvivalCurve = np.zeros(num_times)

        basketTimes[0] = 0.0
        basketSurvivalCurve[0] = 1.0

        for iTime in range(0, num_times):

            t = (payment_dates[iTime] - value_date) / gDaysInYear

            for i_credit in range(0, num_credits):
                issuer_curve = issuer_curves[i_credit]
                recovery_rates[i_credit] = issuer_curve._recovery_rate
                issuerSurvivalProbabilities[i_credit] = interpolate(
                    t, issuer_curve._times, issuer_curve._values,
                    InterpTypes.FLAT_FWD_RATES.value)

            lossDbn = homog_basket_loss_dbn(issuerSurvivalProbabilities,
                                            recovery_rates,
                                            beta_vector,
                                            num_points)

            basketSurvivalCurve[iTime] = 1.0
            for iToDefault in range(n_to_default, num_credits + 1):
                basketSurvivalCurve[iTime] -= lossDbn[iToDefault]

            basketTimes[iTime] = t

        curveRecovery = recovery_rates[0]
        libor_curve = issuer_curves[0]._libor_curve
        basketCurve = CDSCurve(value_date, [], libor_curve, curveRecovery)
        basketCurve._times = basketTimes
        basketCurve._values = basketSurvivalCurve

        protLegPV = self._cds_contract.protection_leg_pv(
            value_date, basketCurve, curveRecovery)
        risky_pv01 = self._cds_contract.risky_pv01(
            value_date, basketCurve)['clean_rpv01']

        # Long protection
        mtm = self._notional * (protLegPV - risky_pv01 * self._running_coupon)

        if not self._long_protection:
            mtm *= -1.0

        basketOutput = np.zeros(4)
        basketOutput[0] = mtm
        basketOutput[1] = risky_pv01 * self._notional * self._running_coupon
        basketOutput[2] = protLegPV * self._notional
        basketOutput[3] = protLegPV / risky_pv01

        return basketOutput

###############################################################################

    def __repr__(self):
        """ print out details of the CDS contract and all of the calculated
        cash flows """
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("STEP-IN DATE", self._step_in_date)
        s += label_to_string("MATURITY", self._maturity_date)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("RUNNING COUPON",
                             self._running_coupon*10000, "bp\n")
        s += label_to_string("DAYCOUNT", self._dc_type)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUSDAYRULE", self._bd_type)
        s += label_to_string("DATEGENRULE", self._dg_type)

#       header = "PAYMENT_DATE, YEAR_FRAC, FLOW"
#       valueTable = [self._payment_dates, self._accrual_factors, self._flows]
#       precision = "12.6f"
#       s += tableToString(header, valueTable, precision)

        return s

###############################################################################
