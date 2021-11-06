##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Add __repr__ method

import numpy as np
from math import sqrt

from ...models.gauss_copula_onefactor import tranch_surv_prob_gaussian
from ...models.gauss_copula_onefactor import tranche_surv_prob_adj_binomial
from ...models.gauss_copula_onefactor import tranche_surv_prob_recursion
from ...models.gauss_copula_lhp import tr_surv_prob_lhp

from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes

from ...products.credit.cds import CDS
from ...products.credit.cds_curve import CDSCurve

from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...market.curves.interpolator import InterpTypes, interpolate
from ...utils.error import FinError

from ...utils.helpers import check_argument_types
from ...utils.date import Date

###############################################################################

from enum import Enum


class FinLossDistributionBuilder(Enum):
    RECURSION = 1
    ADJUSTED_BINOMIAL = 2
    GAUSSIAN = 3
    LHP = 4


###############################################################################


class CDSTranche:

    def __init__(self,
                 step_in_date: Date,
                 maturity_date: Date,
                 k1: float,
                 k2: float,
                 notional: float = ONE_MILLION,
                 running_coupon: float = 0.0,
                 long_protection: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):

        check_argument_types(self.__init__, locals())

        if k1 >= k2:
            raise FinError("K1 must be less than K2")

        self._k1 = k1
        self._k2 = k2

        self._step_in_date = step_in_date
        self._maturity_date = maturity_date
        self._notional = notional
        self._running_coupon = running_coupon
        self._long_protection = long_protection
        self._day_count_type = day_count_type
        self._date_gen_rule_type = date_gen_rule_type
        self._calendar_type = calendar_type
        self._freq_type = freq_type
        self._bus_day_adjust_type = bus_day_adjust_type

        notional = 1.0

        self._cds_contract = CDS(self._step_in_date,
                                 self._maturity_date,
                                 self._running_coupon,
                                 notional,
                                 self._long_protection,
                                 self._freq_type,
                                 self._day_count_type,
                                 self._calendar_type,
                                 self._bus_day_adjust_type,
                                 self._date_gen_rule_type)

    ###############################################################################

    def value_bc(self,
                 valuation_date,
                 issuer_curves,
                 upfront,
                 running_coupon,
                 corr1,
                 corr2,
                 num_points=50,
                 model=FinLossDistributionBuilder.RECURSION):

        num_credits = len(issuer_curves)
        k1 = self._k1
        k2 = self._k2
        tmat = (self._maturity_date - valuation_date) / gDaysInYear

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

        recovery_rates = np.zeros(num_credits)

        payment_dates = self._cds_contract._adjusted_dates
        num_times = len(payment_dates)

        beta1 = sqrt(corr1)
        beta2 = sqrt(corr2)
        beta_vector1 = np.zeros(num_credits)
        for bb in range(0, num_credits):
            beta_vector1[bb] = beta1

        beta_vector2 = np.zeros(num_credits)
        for bb in range(0, num_credits):
            beta_vector2[bb] = beta2

        qVector = np.zeros(num_credits)
        qt1 = np.zeros(num_times)
        qt2 = np.zeros(num_times)
        trancheTimes = np.zeros(num_times)
        trancheSurvivalCurve = np.zeros(num_times)

        trancheTimes[0] = 0
        trancheSurvivalCurve[0] = 1.0
        qt1[0] = 1.0
        qt2[0] = 1.0

        for i in range(1, num_times):

            t = (payment_dates[i] - valuation_date) / gDaysInYear

            for j in range(0, num_credits):
                issuer_curve = issuer_curves[j]
                vTimes = issuer_curve._times
                qRow = issuer_curve._values
                recovery_rates[j] = issuer_curve._recovery_rate
                qVector[j] = interpolate(
                    t, vTimes, qRow, InterpTypes.FLAT_FWD_RATES.value)

            if model == FinLossDistributionBuilder.RECURSION:
                qt1[i] = tranche_surv_prob_recursion(
                    0.0, k1, num_credits, qVector, recovery_rates,
                    beta_vector1, num_points)
                qt2[i] = tranche_surv_prob_recursion(
                    0.0, k2, num_credits, qVector, recovery_rates,
                    beta_vector2, num_points)
            elif model == FinLossDistributionBuilder.ADJUSTED_BINOMIAL:
                qt1[i] = tranche_surv_prob_adj_binomial(
                    0.0, k1, num_credits, qVector, recovery_rates,
                    beta_vector1, num_points)
                qt2[i] = tranche_surv_prob_adj_binomial(
                    0.0, k2, num_credits, qVector, recovery_rates,
                    beta_vector2, num_points)
            elif model == FinLossDistributionBuilder.GAUSSIAN:
                qt1[i] = tranch_surv_prob_gaussian(
                    0.0,
                    k1,
                    num_credits,
                    qVector,
                    recovery_rates,
                    beta_vector1,
                    num_points)
                qt2[i] = tranch_surv_prob_gaussian(
                    0.0,
                    k2,
                    num_credits,
                    qVector,
                    recovery_rates,
                    beta_vector2,
                    num_points)
            elif model == FinLossDistributionBuilder.LHP:
                qt1[i] = tr_surv_prob_lhp(
                    0.0, k1, num_credits, qVector, recovery_rates, beta1)
                qt2[i] = tr_surv_prob_lhp(
                    0.0, k2, num_credits, qVector, recovery_rates, beta2)
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
        libor_curve = issuer_curves[0]._libor_curve
        trancheCurve = CDSCurve(
            valuation_date, [], libor_curve, curveRecovery)
        trancheCurve._times = trancheTimes
        trancheCurve._values = trancheSurvivalCurve

        protLegPV = self._cds_contract.protection_leg_pv(
            valuation_date, trancheCurve, curveRecovery)
        risky_pv01 = self._cds_contract.risky_pv01(
            valuation_date, trancheCurve)['clean_rpv01']

        mtm = self._notional * (protLegPV - upfront -
                                risky_pv01 * running_coupon)

        if not self._long_protection:
            mtm *= -1.0

        trancheOutput = np.zeros(4)
        trancheOutput[0] = mtm
        trancheOutput[1] = risky_pv01 * self._notional * running_coupon
        trancheOutput[2] = protLegPV * self._notional
        trancheOutput[3] = protLegPV / risky_pv01

        return trancheOutput

###############################################################################
