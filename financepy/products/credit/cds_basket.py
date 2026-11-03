##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: There are several speed ups for the Monte-Carlo including calculating
# all default baskets at the same time.

from copy import deepcopy

import numpy as np

from ...utils.error import FinError
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.global_vars import ONE_MILLION
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import label_to_string

from ...models.gauss_copula_onefactor import gauss_homog_basket_loss_dbn_adaptive
# from ...models.gauss_copula_onefactor import homog_basket_loss_dbn

from ...models.gauss_copula import gauss_default_times
from ...models.student_t_copula import student_t_default_times

from ...market.curves.interpolator import interpolate, InterpTypes
from ...market.curves.cds_curve import CDSCurve
from ...products.credit.cds import CDS

from ...utils.check_values import check_curve_dt
from ...utils.global_vars import CLEAN, DIRTY

########################################################################################
# TODO: Implement Student-t recursion
########################################################################################


class CDSBasket:
    """Class to deal with n-to-default CDS baskets."""

    def __init__(
        self,
        step_in_dt: Date,
        maturity_dt: Date,
        num_credits: int,
        n_to_default: int,
        notional: float = ONE_MILLION,
        running_cpn: float = 0.0,
        long_protect: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes | list | tuple = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ) -> None:

        check_argument_types(self.__init__, locals())

        if num_credits < 1:
            raise FinError("Num Credits is zero")

        if n_to_default < 1:
            raise FinError("N-to-Default is zero")

        if n_to_default > num_credits:
            raise FinError("N-to-Default is greater than num_credits")

        if step_in_dt > maturity_dt:
            raise FinError("Step-in date after maturity date")

        self.step_in_dt = step_in_dt
        self.maturity_dt = maturity_dt
        self.num_credits = num_credits
        self.n_to_default = n_to_default
        self.notional = notional
        self.running_cpn = running_cpn / 10000.0
        self.long_protect = long_protect
        self.accrual_dc_type = accrual_dc_type
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
            self.accrual_dc_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

    ####################################################################################

    def value_gaussian_mc(
        self,
        value_dt: Date,
        issuer_curves: list[CDSCurve],
        corr_matrix: np.ndarray,
        num_trials: int,
        seed: int,
        use_moment_matching: bool = False,
        use_sobol: bool = False
    ):
        """Value the default basket using a Gaussian copula model. This
        depends on the issuer discount and correlation matrix."""

        num_credits = len(issuer_curves)

        if num_credits != self.num_credits:
            raise FinError("Num curves not same as num credits")

        if num_trials <= 0:
            raise FinError("Number of Monte Carlo trials must be positive.")

        libor_curve = issuer_curves[0].libor_curve

        check_curve_dt(value_dt, libor_curve)
        check_curve_dt(value_dt, *issuer_curves)

        default_times = gauss_default_times(issuer_curves,
                                            corr_matrix,
                                            num_trials,
                                            seed,
                                            use_moment_matching,
                                            use_sobol)

        dirty_rpv01, prot_leg_pv = self.value_legs_mc(value_dt,
                                                      default_times,
                                                      issuer_curves)

        # Need to calculate Clean RPV01
        day_count = DayCount(self.accrual_dc_type)

        accrual_start_dt = self.cds_contract.accrual_start_dts[0]

        accrued = day_count.year_frac(
            accrual_start_dt,
            value_dt,
        )[0]

        clean_rpv01 = dirty_rpv01 - accrued

        risky_pv01 = np.array([dirty_rpv01, clean_rpv01])

        if clean_rpv01 <= 0.0:
            raise FinError("Cannot compute basket spread: risky PV01 is non-positive."
                           )
        mtm = self.notional * (prot_leg_pv - self.running_cpn * risky_pv01)

        if not self.long_protect:
            mtm = mtm * -1.0

        basket_output = np.zeros(5)
        basket_output[DIRTY] = mtm[DIRTY]
        basket_output[CLEAN] = mtm[CLEAN]
        basket_output[2] = risky_pv01[DIRTY] * self.notional * self.running_cpn
        basket_output[3] = prot_leg_pv * self.notional
        basket_output[4] = prot_leg_pv / risky_pv01[CLEAN]  # This is the fair spread

        return basket_output

    ###########################################################################

    def value_student_t_mc(
        self,
        value_dt: Date,
        issuer_curves: list[CDSCurve],
        corr_matrix: np.ndarray,
        degrees_of_freedom: int,
        num_trials: int,
        seed: int,
        use_moment_matching: bool = False,
        use_sobol: bool = False
    ):
        """Value the default basket using the Student-T copula."""

        num_credits = len(issuer_curves)

        if num_credits != self.num_credits:
            raise FinError("Num curves not same as num credits")

        if num_trials <= 0:
            raise FinError("Number of Monte Carlo trials must be positive.")

        libor_curve = issuer_curves[0].libor_curve

        check_curve_dt(value_dt, libor_curve)
        check_curve_dt(value_dt, *issuer_curves)

        default_times = student_t_default_times(issuer_curves,
                                                corr_matrix,
                                                degrees_of_freedom,
                                                num_trials,
                                                seed,
                                                use_moment_matching,
                                                use_sobol)

        dirty_rpv01, prot_leg_pv = self.value_legs_mc(value_dt,
                                                      default_times,
                                                      issuer_curves)

        # Need to calculate Clean RPV01
        day_count = DayCount(self.accrual_dc_type)

        accrual_start_dt = self.cds_contract.accrual_start_dts[0]

        accrued = day_count.year_frac(
            accrual_start_dt,
            value_dt,
        )[0]

        clean_rpv01 = dirty_rpv01 - accrued

        if clean_rpv01 <= 0.0:
            raise FinError("Cannot compute basket spread: risky PV01 is non-positive.")

        risky_pv01 = np.array([dirty_rpv01, clean_rpv01])

        mtm = self.notional * (prot_leg_pv - self.running_cpn * risky_pv01)

        if not self.long_protect:
            mtm = mtm * -1.0

        basket_output = np.zeros(5)
        basket_output[DIRTY] = mtm[DIRTY]
        basket_output[CLEAN] = mtm[CLEAN]
        basket_output[2] = risky_pv01[DIRTY] * self.notional * self.running_cpn
        basket_output[3] = prot_leg_pv * self.notional
        basket_output[4] = prot_leg_pv / risky_pv01[CLEAN]  # This is the fair spread

        return basket_output

    ###########################################################################

    def value_1f_gaussian_homo(
        self,
        value_dt: Date,
        issuer_curves: list[CDSCurve],
        beta_vector: np.ndarray,
        num_points: int = 25,
    ):
        """Value default basket using 1 factor Gaussian copula and analytical
        approach which is only exact when all recovery rates are the same."""

        num_credits = len(issuer_curves)

        if num_credits != self.num_credits:
            raise FinError("Num curves not same as num credits")

        libor_curve = issuer_curves[0].libor_curve

        check_curve_dt(value_dt, libor_curve)
        check_curve_dt(value_dt, *issuer_curves)

        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEAR

        if t_mat < 0.0:
            raise FinError("Value date is after maturity date")

        payment_dts = self.cds_contract.payment_dts
        num_times = len(payment_dts)

        libor_curve = issuer_curves[0].libor_curve

        issuer_surv_probs = np.zeros(num_credits)
        recovery_rates = np.zeros(num_credits)
        basket_times = np.zeros(num_times)
        basket_surv_curve = np.zeros(num_times)

        basket_times[0] = 0.0
        basket_surv_curve[0] = 1.0

        for i_credit in range(0, num_credits):
            issuer_curve = issuer_curves[i_credit]
            recovery_rates[i_credit] = issuer_curve.recovery_rate

        if not np.allclose(
            recovery_rates,
            recovery_rates[0],
            atol=1e-10,
            rtol=0.0
        ):
            raise FinError(
                "Analytical basket valuation requires "
                "homogeneous recovery rates."
            )

        for i_time in range(0, num_times):

            t = (payment_dts[i_time] - value_dt) / G_DAYS_IN_YEAR

            if t <= 1e-12:
                continue

            for i_credit in range(0, num_credits):

                issuer_curve = issuer_curves[i_credit]
                issuer_surv_probs[i_credit] = interpolate(
                    t,
                    issuer_curve._times,
                    issuer_curve._qs,
                    InterpTypes.FLAT_FWD_RATES.value,
                )

            loss_dbn = gauss_homog_basket_loss_dbn_adaptive(issuer_surv_probs,
                                                            recovery_rates,
                                                            beta_vector,
                                                            num_points)

            basket_surv_curve[i_time] = 1.0
            for i_to_default in range(self.n_to_default, num_credits + 1):
                basket_surv_curve[i_time] -= loss_dbn[i_to_default]

            basket_times[i_time] = t

        curve_recovery = recovery_rates[0]
        basket_curve = CDSCurve(value_dt, [], libor_curve, curve_recovery)
        basket_curve._times = basket_times
        basket_curve._qs = basket_surv_curve

        prot_leg_pv = self.cds_contract.prot_leg_pv(value_dt,
                                                    basket_curve,
                                                    curve_recovery)

        # Dirty and clean risky PV01
        risky_pv01 = np.asarray(
            self.cds_contract.rpv01(value_dt, basket_curve)
        )

        if np.any(risky_pv01 <= 0.0):
            raise FinError("Cannot compute basket spread: risky PV01 is non-positive.")

        # Dirty and clean MTMs
        mtm = self.notional * (
            prot_leg_pv - risky_pv01 * self.running_cpn
        )

        if not self.long_protect:
            mtm *= -1.0

        basket_output = np.zeros(5)
        basket_output[DIRTY] = mtm[DIRTY]
        basket_output[CLEAN] = mtm[CLEAN]
        basket_output[2] = risky_pv01[DIRTY] * self.notional * self.running_cpn
        basket_output[3] = prot_leg_pv * self.notional
        basket_output[4] = prot_leg_pv / risky_pv01[CLEAN]

        return basket_output

    ###########################################################################

    def value_legs_mc(
        self,
        value_dt: Date,
        default_times: np.ndarray,
        issuer_curves: list[CDSCurve],
    ):
        """Vectorised MC valuation of n-to-default basket legs."""

        num_credits, num_trials = default_times.shape

        if num_credits != self.num_credits:
            raise FinError("Num curves not same as num credits")

        if not 1 <= self.n_to_default <= num_credits:
            raise FinError("Invalid n_to_default")

        if self.num_credits != len(issuer_curves):
            raise FinError("Default-time dimension does not match issuer curves.")

        if num_trials == 0:
            raise FinError("Number of Monte Carlo trials must be positive.")

        payment_dts = self.cds_contract.payment_dts
        day_count = DayCount(self.accrual_dc_type)

        # Schedule times
        pay_times = np.array([
            (dt - value_dt) / G_DAYS_IN_YEAR
            for dt in payment_dts
        ])

        accrual_start_dt = self.cds_contract.accrual_start_dts[0]
        accrual_start_time = (
            accrual_start_dt - value_dt
        ) / G_DAYS_IN_YEAR

        accrual_factors = np.empty(len(payment_dts))

        accrual_factors[0] = day_count.year_frac(
            accrual_start_dt, payment_dts[0]
        )[0]

        for i in range(1, len(payment_dts)):
            accrual_factors[i] = day_count.year_frac(
                payment_dts[i - 1], payment_dts[i]
            )[0]

        libor_curve = issuer_curves[0].libor_curve
        df_pay = libor_curve.df_t(pay_times)

        cum_pv01 = np.cumsum(accrual_factors * df_pay)
        full_risky_pv01 = cum_pv01[-1]

        t_mat = (
            self.maturity_dt - value_dt
        ) / G_DAYS_IN_YEAR

        recoveries = np.array([
            curve.recovery_rate for curve in issuer_curves
        ])

        # Find nth defaulter in every simulation.
        path_indices = np.arange(num_trials)

        if self.n_to_default == 1:
            nth_indices = np.argmin(default_times, axis=0)

        elif self.n_to_default == num_credits:
            nth_indices = np.argmax(default_times, axis=0)

        else:
            nth_indices = np.argpartition(
                default_times,
                self.n_to_default - 1,
                axis=0
            )[self.n_to_default - 1, :]

        tau = default_times[nth_indices, path_indices]

        # Start with the full premium annuity and zero protection.
        risky_pv01_paths = np.full(num_trials, full_risky_pv01)
        prot_paths = np.zeros(num_trials)

        # Only paths triggered before maturity need further work.
        triggered = tau < t_mat

        if np.any(triggered):

            tau_def = tau[triggered]
            nth_def = nth_indices[triggered]

            # Number of scheduled payment dates before default.
            pos = np.searchsorted(
                pay_times, tau_def, side="left"
            )

            previous = np.maximum(pos - 1, 0)

            scheduled_pv01 = np.where(
                pos > 0,
                cum_pv01[previous],
                0.0
            )

            last_start_time = np.where(
                pos > 0,
                pay_times[previous],
                accrual_start_time
            )

            partial_accrual = np.maximum(
                0.0,
                (tau_def - last_start_time)
                * G_DAYS_IN_YEAR / 360.0
            )

            # vectorised discount factors for all default times
            df_def = libor_curve.df_t(tau_def)

            risky_pv01_paths[triggered] = (
                scheduled_pv01 + partial_accrual * df_def
            )

            prot_paths[triggered] = (
                (1.0 - recoveries[nth_def]) * df_def
            )

        return np.mean(risky_pv01_paths), np.mean(prot_paths)

    ###########################################################################

    def spread_dv01s(self, value_dt, issuer_curves, beta_vector):

        v0 = self.value_1f_gaussian_homo(
            value_dt, issuer_curves, beta_vector
        )

        bump = 0.0001
        spread_dv01_vector = []

        for i_credit in range(self.num_credits):

            bumped_curves = list(issuer_curves)
            bumped_curve = deepcopy(issuer_curves[i_credit])

            for cds in bumped_curve.cds_contracts:
                cds.running_cpn += bump

            bumped_curve.build_curve()
            bumped_curves[i_credit] = bumped_curve

            v1 = self.value_1f_gaussian_homo(
                value_dt, bumped_curves, beta_vector
            )

            spread_dv01_vector.append(v1[0] - v0[0])

        return spread_dv01_vector

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all of the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("STEP-IN DATE", self.step_in_dt)
        s += label_to_string("MATURITY", self.maturity_dt)
        s += label_to_string("NUM CREDITS", self.num_credits)
        s += label_to_string("N-TO-DEFAULT", self.n_to_default)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("RUNNING COUPON", self.running_cpn * 10000, " bp\n")
        s += label_to_string("DAY COUNT", self.accrual_dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST", self.bd_type)
        s += label_to_string("DATE GEN RULE", self.dg_type)

        #  header = "PAYMENT_DT, YEAR_FRAC, PAYMENT"
        #  value_table = [self.payment_dts, self.accrual_factors, self.flows]
        #  precision = "12.6f"
        #  s += tableToString(header, value_table, precision)

        return s


########################################################################################
