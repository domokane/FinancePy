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
from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import label_to_string

from ...models.gauss_copula_onefactor import homog_basket_loss_dbn
from ...models.gauss_copula import default_times_gc
from ...models.student_t_copula import StudentTCopula

from ...market.curves.interpolator import interpolate, InterpTypes

from ...products.credit.cds_curve import CDSCurve
from ...products.credit.cds import CDS

###############################################################################
# TODO: Convert functions to use NUMBA!!
###############################################################################


class CDSBasket:
    """Class to deal with n-to-default CDS baskets."""

    def __init__(
        self,
        step_in_dt: Date,
        maturity_dt: Date,
        notional: float = ONE_MILLION,
        running_cpn: float = 0.0,
        long_protect: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):

        check_argument_types(self.__init__, locals())

        self.step_in_dt = step_in_dt
        self.maturity_dt = maturity_dt
        self.notional = notional
        self.running_cpn = running_cpn / 10000.0
        self.long_protect = long_protect
        self.dc_type = dc_type
        self.dg_type = dg_type
        self.cal_type = cal_type
        self.freq_type = freq_type
        self.bd_type = bd_type

        self.cds_contract = CDS(
            self.step_in_dt,
            self.maturity_dt,
            self.running_cpn,
            1.0,
            self.long_protect,
            self.freq_type,
            self.dc_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

    ###########################################################################

    def value_legs_mc(
        self, value_dt, n_to_default, default_times, issuer_curves, libor_curve
    ):
        """Value the legs of the default basket using Monte Carlo. The default
        times are an input so this valuation is not model dependent."""

        num_credits = default_times.shape[0]
        num_trials = default_times.shape[1]

        payment_dts = self.cds_contract.payment_dts
        num_payments = len(payment_dts)
        day_count = DayCount(self.dc_type)

        avg_acc_factor = 0.0

        rpv01_to_times = np.zeros(num_payments)

        for i_time in range(1, num_payments):

            t = (payment_dts[i_time] - value_dt) / g_days_in_year
            dt0 = payment_dts[i_time - 1]
            dt1 = payment_dts[i_time]
            accrual_factor = day_count.year_frac(dt0, dt1)[0]
            avg_acc_factor += accrual_factor
            rpv01_to_times[i_time] = rpv01_to_times[
                i_time - 1
            ] + accrual_factor * libor_curve.df_t(t)

        avg_acc_factor /= num_payments

        t_mat = (self.maturity_dt - value_dt) / g_days_in_year

        rpv01 = 0.0
        prot = 0.0

        asset_tau = np.zeros(num_credits)

        for i_trial in range(0, num_trials):

            for i_credit in range(0, num_credits):

                asset_tau[i_credit] = default_times[i_credit, i_trial]

            # ORDER THE DEFAULT TIMES
            asset_tau.sort()

            # GET THE Nth DEFAULT TIME
            min_tau = asset_tau[n_to_default - 1]

            if min_tau < t_mat:

                num_payment_amounts_index = int(min_tau / avg_acc_factor)
                rpv01_trial = rpv01_to_times[num_payment_amounts_index]
                rpv01_trial += (
                    min_tau - num_payment_amounts_index * avg_acc_factor
                )

                # DETERMINE IDENTITY OF N-TO-DEFAULT CREDIT IF BASKET NOT HOMO
                asset_index = 0
                for i_credit in range(0, num_credits):
                    if min_tau == default_times[i_credit, i_trial]:
                        asset_index = i_credit
                        break

                prot_trial = 1.0 - issuer_curves[asset_index].recovery_rate
                prot_trial *= libor_curve.df_t(min_tau)

            else:

                num_payment_amounts_index = int(t_mat / avg_acc_factor)
                rpv01_trial = rpv01_to_times[-1]
                prot_trial = 0.0

            rpv01 += rpv01_trial
            prot += prot_trial

        rpv01 = rpv01 / num_trials
        prot = prot / num_trials
        return (rpv01, prot)

    ###########################################################################

    def value_gaussian_mc(
        self,
        value_dt,
        n_to_default,
        issuer_curves,
        corr_matrix,
        libor_curve,
        num_trials,
        seed,
    ):
        """Value the default basket using a Gaussian copula model. This
        depends on the issuer discount and correlation matrix."""

        num_credits = len(issuer_curves)

        if n_to_default > num_credits or n_to_default < 1:
            raise FinError("n_to_default must be 1 to num_credits")

        default_times = default_times_gc(
            issuer_curves, corr_matrix, num_trials, seed
        )

        rpv01, prot_pv = self.value_legs_mc(
            value_dt, n_to_default, default_times, issuer_curves, libor_curve
        )

        spd = prot_pv / rpv01
        value = self.notional * (prot_pv - self.running_cpn * rpv01)

        if not self.long_protect:
            value = value * -1.0

        return (value, rpv01, spd)

    ###########################################################################

    def value_student_t_mc(
        self,
        value_dt,
        n_to_default,
        issuer_curves,
        corr_matrix,
        degrees_of_freedom,
        libor_curve,
        num_trials,
        seed,
    ):
        """Value the default basket using the Student-T copula."""

        num_credits = len(issuer_curves)

        if n_to_default > num_credits or n_to_default < 1:
            raise FinError("n_to_default must be 1 to num_credits")

        model = StudentTCopula()

        default_times = model.default_times(
            issuer_curves, corr_matrix, degrees_of_freedom, num_trials, seed
        )

        rpv01, prot_pv = self.value_legs_mc(
            value_dt, n_to_default, default_times, issuer_curves, libor_curve
        )

        spd = prot_pv / rpv01
        value = self.notional * (prot_pv - self.running_cpn * rpv01)

        if not self.long_protect:
            value = value * -1.0

        return (value, rpv01, spd)

    ###########################################################################

    def value_1f_gaussian_homo(
        self,
        value_dt,
        n_to_default,
        issuer_curves,
        beta_vector,
        libor_curve,
        num_points=50,
    ):
        """Value default basket using 1 factor Gaussian copula and analytical
        approach which is only exact when all recovery rates are the same."""

        num_credits = len(issuer_curves)

        if num_credits == 0:
            raise FinError("Num Credits is zero")

        if n_to_default < 1 or n_to_default > num_credits:
            raise FinError("n_to_default must be 1 to num_credits")

        t_mat = (self.maturity_dt - value_dt) / g_days_in_year

        if t_mat < 0.0:
            raise FinError("Value date is after maturity date")

        payment_dts = self.cds_contract.payment_dts
        num_times = len(payment_dts)

        issuer_surv_probs = np.zeros(num_credits)
        recovery_rates = np.zeros(num_credits)
        basket_times = np.zeros(num_times)
        basket_surv_curve = np.zeros(num_times)

        basket_times[0] = 0.0
        basket_surv_curve[0] = 1.0

        for i_time in range(0, num_times):

            t = (payment_dts[i_time] - value_dt) / g_days_in_year

            for i_credit in range(0, num_credits):
                issuer_curve = issuer_curves[i_credit]
                recovery_rates[i_credit] = issuer_curve.recovery_rate
                issuer_surv_probs[i_credit] = interpolate(
                    t,
                    issuer_curve._times,
                    issuer_curve._values,
                    InterpTypes.FLAT_FWD_RATES.value,
                )

            loss_dbn = homog_basket_loss_dbn(
                issuer_surv_probs, recovery_rates, beta_vector, num_points
            )

            basket_surv_curve[i_time] = 1.0
            for i_to_default in range(n_to_default, num_credits + 1):
                basket_surv_curve[i_time] -= loss_dbn[i_to_default]

            basket_times[i_time] = t

        curve_recovery = recovery_rates[0]
        libor_curve = issuer_curves[0].libor_curve
        basket_curve = CDSCurve(value_dt, [], libor_curve, curve_recovery)
        basket_curve._times = basket_times
        basket_curve._values = basket_surv_curve

        prot_leg_pv = self.cds_contract.prot_leg_pv(
            value_dt, basket_curve, curve_recovery
        )
        risky_pv01 = self.cds_contract.risky_pv01(value_dt, basket_curve)[
            "clean_rpv01"
        ]

        # Long protection
        mtm = self.notional * (prot_leg_pv - risky_pv01 * self.running_cpn)

        if not self.long_protect:
            mtm *= -1.0

        basket_output = np.zeros(4)
        basket_output[0] = mtm
        basket_output[1] = risky_pv01 * self.notional * self.running_cpn
        basket_output[2] = prot_leg_pv * self.notional
        basket_output[3] = prot_leg_pv / risky_pv01

        return basket_output

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all of the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("STEP-IN DATE", self.step_in_dt)
        s += label_to_string("MATURITY", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string(
            "RUNNING COUPON", self.running_cpn * 10000, "bp\n"
        )
        s += label_to_string("DAYCOUNT", self.dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUSDAYRULE", self.bd_type)
        s += label_to_string("DATEGENRULE", self.dg_type)

        #  header = "PAYMENT_dt, YEAR_FRAC, FLOW"
        #  value_table = [self.payment_dts, self.accrual_factors, self.flows]
        #  precision = "12.6f"
        #  s += tableToString(header, value_table, precision)

        return s


###############################################################################
