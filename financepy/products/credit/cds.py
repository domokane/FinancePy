##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

from copy import deepcopy
from math import exp, log


import numpy as np
from numba import njit, float64, int64

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.calendar import Calendar, CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.math import ONE_MILLION
from ...utils.helpers import label_to_string, table_to_string
from ...market.curves.interpolator import InterpTypes, _uinterpolate

from ...utils.helpers import check_argument_types

USE_FLAT_HAZARD_RATE_INTEGRAL = True
STANDARD_RECOVERY_RATE = 0.40
GLOB_NUM_STEPS_PER_YEAR = 25
ONE_BP = 0.0001
ONE_PCT = 0.01

# Premium accrues on an Act/360 basis while time is measured in
# calendar years, hence the 365/360 adjustment to the clean RPV01.
KAPPA = 365.0 / 360.0

# RETURN ORDER
DIRTY = 0
CLEAN = 1

########################################################################################
# TODO: Perform protection leg pv analytically using fact that hazard rate and
#       interest rates are flat between their combined node points. Right now I
#       do not find the protection leg PV calculations to be a bottleneck,
#       especially given the speedup benefits of using NUMBA.
########################################################################################


@njit(
    float64[:](
        float64,
        float64,
        float64[:],
        float64[:],
        float64[:],
        float64[:],
        float64[:],
        float64[:],
        int64,
    ),
    fastmath=True,
    cache=True,
)
def _risky_pv01_numba(
    t_eff,
    accrual_factor_pcd_to_now,
    payment_times,
    year_fracs,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    pv01_method,
):
    """Fast calculation of the risky PV01 of a CDS using NUMBA.
    The output is a numpy array of the full and clean risky PV01."""

    method = InterpTypes.FLAT_FWD_RATES.value
    debug = False

    if debug:
        print("===================")
        print("t_eff", t_eff)
        print("Acc", accrual_factor_pcd_to_now)
        print("Payments", payment_times)
        print("Alphas", year_fracs)
        print("QTimes", np_surv_times)
        print("QValues", np_surv_values)

    cpn_accd_indicator = 1

    # Method 0 : This is the market standard which assumes that the cpn
    # accrued is treated as though on average default occurs roughly midway
    # through a cpn period.

    t_ncd = payment_times[0]

    # The first cpn is a special case which needs to be handled carefully
    # taking into account what cpn has already accrued and what has not
    qeff = _uinterpolate(t_eff, np_surv_times, np_surv_values, method)
    q1 = _uinterpolate(t_ncd, np_surv_times, np_surv_values, method)
    z1 = _uinterpolate(t_ncd, np_ibor_times, np_ibor_values, method)

    # this is the part of the cpn accrued from previous cpn date to now
    # accrual_factor_pcd_to_now = day_count.year_frac(pcd,teff)
    # reference credit survives to the premium payment date
    dirty_rpv01 = q1 * z1 * year_fracs[1]

    # cpn accrued from previous cpn to today paid in full at default
    # before cpn payment
    dq = qeff - q1
    dirty_rpv01 += z1 * dq * accrual_factor_pcd_to_now * cpn_accd_indicator

    # future accrued from now to cpn payment date assuming default roughly
    # midway
    dirty_rpv01 += (
        0.5
        * z1
        * dq
        * (year_fracs[1] - accrual_factor_pcd_to_now)
        * cpn_accd_indicator
    )

    for it in range(1, len(payment_times)):

        t2 = payment_times[it]
        q2 = _uinterpolate(t2, np_surv_times, np_surv_values, method)
        z2 = _uinterpolate(t2, np_ibor_times, np_ibor_values, method)
        accrual_factor = year_fracs[it]

        # full cpn is paid at the end of the current period if survives
        dirty_rpv01 += q2 * z2 * accrual_factor

        #        print(it, t2, z2, q2, accrual_factor, full_rpv01)

        #######################################################################

        if cpn_accd_indicator == 1:

            if USE_FLAT_HAZARD_RATE_INTEGRAL:
                # This needs to be updated to handle small h+r
                tau = accrual_factor
                h12 = -log(q2 / q1) / tau
                r12 = -log(z2 / z1) / tau
                alpha = h12 + r12
                exp_term = (
                    1.0 - exp(-alpha * tau) - alpha * tau * exp(-alpha * tau)
                )
                d_dirty_rpv01 = (
                    q1 * z1 * h12 * exp_term / abs(alpha * alpha + 1e-20)
                )
            else:
                d_dirty_rpv01 = 0.50 * (q1 - q2) * z2 * accrual_factor

            dirty_rpv01 = dirty_rpv01 + d_dirty_rpv01

        q1 = q2
        z1 = z2

    clean_rpv01 = dirty_rpv01 - accrual_factor_pcd_to_now

    return np.array([dirty_rpv01, clean_rpv01])


########################################################################################


@njit(
    float64(
        float64,
        float64,
        float64[:],
        float64[:],
        float64[:],
        float64[:],
        float64,
        int64,
        int64,
    ),
    fastmath=True,
    cache=True,
)
def _prot_leg_pv_numba(
    t_eff,
    t_mat,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    contract_recovery_rate,
    num_steps_per_year,
    prot_method,
):
    """Fast calculation of the CDS protection leg PV using NUMBA to speed up
    the numerical integration over time."""

    if t_eff < 0.0:
        raise FinError("Error: Protection leg starts in past: t_eff < 0")

    method = InterpTypes.FLAT_FWD_RATES.value
    dt = 1.0 / num_steps_per_year
    num_steps = int((t_mat - t_eff) * num_steps_per_year + 0.50)
    dt = (t_mat - t_eff) / num_steps

    t = t_eff
    z1 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
    q1 = _uinterpolate(t, np_surv_times, np_surv_values, method)

    prot_pv = 0.0
    small = 1e-8

    if USE_FLAT_HAZARD_RATE_INTEGRAL:

        log_z1 = np.log(z1)
        log_q1 = np.log(q1)

        for _ in range(0, num_steps):
            t = t + dt
            z2 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
            q2 = _uinterpolate(t, np_surv_times, np_surv_values, method)

            log_z2 = np.log(z2)
            log_q2 = np.log(q2)

            # This needs to be updated to handle small h+r
            # h12 = -log(q2 / q1) / dt
            # r12 = -log(z2 / z1) / dt

            h12 = -(log_q2 - log_q1) / dt
            r12 = -(log_z2 - log_z1) / dt

            exp_term = exp(-(r12 + h12) * dt)
            dprot_pv = (
                h12 * (1.0 - exp_term) * q1 * z1 / (abs(h12 + r12) + small)
            )
            prot_pv += dprot_pv

            q1, z1 = q2, z2
            log_q1, log_z1 = log_q2, log_z2

    else:

        for _ in range(0, num_steps):
            t += dt
            z2 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
            q2 = _uinterpolate(t, np_surv_times, np_surv_values, method)
            dq = q1 - q2
            dprot_pv = 0.5 * (z1 + z2) * dq
            prot_pv += dprot_pv
            q1 = q2
            z1 = z2

    prot_pv = prot_pv * (1.0 - contract_recovery_rate)
    return prot_pv


########################################################################################
########################################################################################
########################################################################################


class CDS:
    """A class which manages a Credit Default Swap. It performs schedule
    generation and the valuation and risk management of CDS."""

    def __init__(
        self,
        step_in_dt: Date,  # Date protection starts
        maturity_dt_or_tenor: Union[Date, str],  # Date or tenor
        running_cpn: float,  # Annualised cpn on premium fee leg
        notional: float = ONE_MILLION,
        long_protect: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create a CDS from the step-in date, maturity date and cpn"""

        check_argument_types(self.__init__, locals())

        if isinstance(maturity_dt_or_tenor, Date):
            maturity_dt = maturity_dt_or_tenor
        else:
            # To get the next CDS date we move on by the tenor and then roll to
            # the next CDS date after that. We do not holiday adjust it. That
            # is handled in the schedule generation.
            maturity_dt = step_in_dt.add_tenor(maturity_dt_or_tenor)
            maturity_dt = maturity_dt.next_cds_date()

        if step_in_dt > maturity_dt:
            raise FinError("Step in date after maturity date")

        self.step_in_dt = step_in_dt
        self.maturity_dt = maturity_dt
        self.running_cpn = running_cpn
        self.notional = notional
        self.long_protect = long_protect
        self.accrual_dc_type = accrual_dc_type
        self.dg_type = dg_type
        self.cal_type = cal_type
        self.freq_type = freq_type
        self.bd_type = bd_type

        self._generate_adjusted_cds_payment_dts()
        self._calc_flows()

    ####################################################################################

    def _generate_adjusted_cds_payment_dts(self):
        """Generate CDS payment dates which have been holiday adjusted."""

        frequency = annual_frequency(self.freq_type)
        calendar = Calendar(self.cal_type)
        start_dt = self.step_in_dt

        self.payment_dts = []
        self.accrual_start_dts = []
        self.accrual_end_dts = []
        num_months = int(12.0 / frequency)

        # We generate unadjusted dates - not adjusted for weekends or holidays
        unadjusted_schedule_dts = []

        if self.dg_type == DateGenRuleTypes.BACKWARD:

            # We start at end date and step backwards
            next_dt = self.maturity_dt

            unadjusted_schedule_dts.append(next_dt)

            # the unadjusted dates start at end date and end at previous cpn date
            while next_dt > start_dt:
                next_dt = next_dt.add_months(-num_months)
                unadjusted_schedule_dts.append(next_dt)

            # now we adjust for holiday using business day adjustment
            # convention specified
            adjusted_dts = []

            for date in reversed(unadjusted_schedule_dts):
                adjusted = calendar.adjust(date, self.bd_type)
                adjusted_dts.append(adjusted)

        # https://www.cdsmodel.com/assets/cds-model/docs/Standard%20CDS%20Examples.pdf
        # Payment   = [20-MAR-2009, 22-JUN-2009, 21-SEP-2009, 21-DEC-2009, 22-MAR-2010]
        # Acc Start = [22-DEC-2008, 20-MAR-2009, 22-JUN-2009, 21-SEP-2009, 21-DEC-2009]
        # Acc End   = [19-MAR-2009, 21-JUN-2009, 20-SEP-2009, 20-DEC-2009, 20-MAR-2010]

        elif self.dg_type == DateGenRuleTypes.FORWARD:

            # We start at start date and step forwards

            next_dt = start_dt

            # the unadjusted dates start at start date and end at last date
            # before maturity date
            while next_dt < self.maturity_dt:
                unadjusted_schedule_dts.append(next_dt)
                next_dt = next_dt.add_months(num_months)

            # We then append the maturity date
            unadjusted_schedule_dts.append(self.maturity_dt)

            adjusted_dts = []
            for date in unadjusted_schedule_dts:
                adjusted = calendar.adjust(date, self.bd_type)
                adjusted_dts.append(adjusted)

        # eg. Date(20, 2, 2009) to Date(20, 3, 2010) with DateGenRuleTypes.FORWARD
        # Payment       = [20-MAY-2009, 20-AUG-2009, 20-NOV-2009, 22-FEB-2010]
        # Accrual Start = [20-FEB-2009, 20-MAY-2009, 20-AUG-2009, 20-NOV-2009]
        # Accrual End   = [19-MAY-2009, 19-AUG-2009, 19-NOV-2009, 20-MAR-2010]

        else:
            raise FinError("Unknown DATE_GEN_RULEType:" + str(self.dg_type))

        # We only include dates which fall after the CDS start date
        self.payment_dts = adjusted_dts[1:]

        # Accrual start dates run from previous cpn date to penultimate cpn date
        self.accrual_start_dts = adjusted_dts[:-1]

        # Accrual end dates are one day before the start of the next accrual period
        self.accrual_end_dts = [
            date.add_days(-1) for date in self.accrual_start_dts[1:]
        ]

        # Final accrual end date is the maturity date
        self.accrual_end_dts.append(self.maturity_dt)

    ###########################################################################

    def _calc_flows(self):
        """Calculate cash flow amounts on premium leg."""
        day_count = DayCount(self.accrual_dc_type)

        self.accrual_factors = []
        self.flows = []

        # self.accrual_factors.append(0.0)
        # self.flows.append(0.0)

        for t0, t1 in zip(self.accrual_start_dts, self.accrual_end_dts):
            # Adding a day because `year_frac` is non-inclusive
            # eg. 20th to 22nd should be 3 days
            accrual_factor = day_count.year_frac(t0, t1.add_days(1))[0]
            flow = accrual_factor * self.running_cpn * self.notional

            self.accrual_factors.append(accrual_factor)
            self.flows.append(flow)

    ###########################################################################

    def value(
        self,
        value_dt,
        issuer_curve,
        contract_recovery_rate,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Valuation of a CDS contract on a specific valuation date given
        an issuer curve and a contract recovery rate."""

        rpv01 = self.rpv01(value_dt, issuer_curve, pv01_method)

        dirty_rpv01 = rpv01[DIRTY]
        clean_rpv01 = rpv01[CLEAN]

        prot_pv = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            num_steps_per_year,
            prot_method,
        )

        fwd_df = 1.0

        if self.long_protect:
            omega = +1
        else:
            omega = -1

        dirty_pv = (
            fwd_df
            * omega
            * (prot_pv - self.running_cpn * dirty_rpv01 * self.notional)
        )

        clean_pv = (
            fwd_df
            * omega
            * (prot_pv - self.running_cpn * clean_rpv01 * self.notional)
        )

        return (dirty_pv, clean_pv)

    ###########################################################################

    def spread_dv01(
        self,
        value_dt,
        issuer_curve,
        contract_recovery_rate,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Calculation of the change in the value of the CDS contract for a
        one basis point change in the level of the CDS curve."""

        v0 = self.value(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        bump = 0.0001  # 1 basis point

        # we create a deep copy to avoid state issues
        bumped_issuer_curve = deepcopy(issuer_curve)
        for cds in bumped_issuer_curve.cds_contracts:
            cds.running_cpn += bump

        bumped_issuer_curve.build_curve()

        v1 = self.value(
            value_dt,
            bumped_issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        credit_dv01 = v1[DIRTY] - v0[DIRTY]
        return credit_dv01

    ###########################################################################

    def ir_dv01(
        self,
        value_dt: Date,
        issuer_curve,
        contract_recovery_rate,
        pv01_method: int = 0,
        prot_method: int = 0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Calculation of the interest DV01 based on a simple bump of
        the discount factors and reconstruction of the CDS curve."""

        v0 = self.value(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        # we create a deep copy to avoid state issues
        bumped_issuer_curve = deepcopy(issuer_curve)
        original_discount_curve = bumped_issuer_curve.libor_curve

        if not hasattr(original_discount_curve, "bump_parallel"):
            raise FinError(
            f"Discount curve type "
            f"{type(original_discount_curve).__name__} "
            "does not implement bump_parallel()."
        )

        bumped_issuer_curve.libor_curve = (original_discount_curve.bump_parallel(ONE_BP))
        bumped_issuer_curve.build_curve()

        v1 = self.value(
            value_dt,
            bumped_issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        interest_dv01 = v1[DIRTY] - v0[DIRTY]
        return interest_dv01

    ###########################################################################

    def recovery_dv01(
        self,
        value_dt: Date,
        issuer_curve,
        contract_recovery_rate: float,
        pv01_method: int = 0,
        prot_method: int = 0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """PV change when contract and curve recovery both increase by 1%.

        The market CDS calibration quotes and discount curve are held fixed.

        Returns:
            PV(curve recovery + 0.01, contract recovery + 0.01)
            - PV(curve recovery, contract recovery)
        """

        bump = ONE_PCT  # 0.01

        if not 0.0 <= contract_recovery_rate < 1.0:
            raise FinError(
                "Contract recovery rate must lie in [0.0, 1.0)."
            )

        if contract_recovery_rate + bump >= 1.0:
            raise FinError(
                "Bumped contract recovery rate must be less than 1.0."
            )

        if not hasattr(issuer_curve, "recovery_rate"):
            raise FinError(
                f"Issuer curve type {type(issuer_curve).__name__} "
                "does not expose recovery_rate."
            )

        curve_recovery_rate = issuer_curve.recovery_rate

        if not 0.0 <= curve_recovery_rate < 1.0:
            raise FinError(
                "Issuer curve recovery rate must lie in [0.0, 1.0)."
            )

        if curve_recovery_rate + bump >= 1.0:
            raise FinError(
                "Bumped issuer curve recovery rate must be less than 1.0."
            )

        # Base valuation:
        # - original survival curve
        # - original contract recovery
        value_base = self.value(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        # Bump the recovery used to calibrate the survival curve.
        bumped_issuer_curve = deepcopy(issuer_curve)
        bumped_issuer_curve.recovery_rate = curve_recovery_rate + bump

        # Rebuild survival probabilities while holding CDS quotes and
        # the discount curve unchanged.
        bumped_issuer_curve.build_curve()

        # Value using both:
        # - bumped survival curve
        # - bumped contract recovery
        value_bumped = self.value(
            value_dt,
            bumped_issuer_curve,
            contract_recovery_rate + bump,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        return value_bumped[DIRTY] - value_base[DIRTY]

    ###########################################################################

    def upfront(
        self,
        value_dt,
        settle_dt,
        issuer_curve,
        contract_recovery_rate,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Return clean percentage PV expressed on T+3 settlement date."""

        discount_curve = issuer_curve.libor_curve

        if settle_dt < discount_curve.value_dt:
            raise FinError(
                "Settlement date cannot precede the discount curve date: "
                f"settle_dt={settle_dt}, "
                f"curve.value_dt={discount_curve.value_dt}"
            )

        values = self.value(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        settlement_df = discount_curve.df(settle_dt)

        if not np.isfinite(settlement_df) or settlement_df <= 0.0:
            raise FinError(
            "Settlement-date discount factor must be finite and positive: "
            f"df={settlement_df}"
        )

        upfront = values[CLEAN] / settlement_df / self.notional
        return upfront

    ###########################################################################

    def cash_settlement_amount(
        self,
        value_dt,
        settle_dt,
        issuer_curve,
        contract_recovery_rate,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Return dirty amount paid on the T+3 settlement date."""

        discount_curve = issuer_curve.libor_curve

        if settle_dt < discount_curve.value_dt:
            raise FinError(
                "Settlement date cannot precede the discount curve date: "
                f"settle_dt={settle_dt}, "
                f"curve.value_dt={discount_curve.value_dt}"
            )

        values = self.value(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

        settlement_df = discount_curve.df(settle_dt)

        if not np.isfinite(settlement_df) or settlement_df <= 0.0:
            raise FinError(
            "Settlement-date discount factor must be finite and positive: "
            f"df={settlement_df}"
        )

        csa_pv = values[DIRTY] / settlement_df
        return csa_pv

    ###########################################################################

    def clean_price(
        self,
        value_dt,
        issuer_curve,
        contract_recovery_rate,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
    ):
        """Value of the CDS contract excluding accrued interest."""

        risky_pv01 = self.rpv01(value_dt, issuer_curve, pv01_method)

        clean_rpv01 = risky_pv01[CLEAN]

        prot_pv = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            num_steps_per_year,
            prot_method,
        )

        fwd_df = 1.0

        clean_pv = fwd_df * (
            prot_pv - self.running_cpn * clean_rpv01 * self.notional
        )

        clean_price = (self.notional - clean_pv) / self.notional * 100.0

        return clean_price

    ###########################################################################

    def accrued_days(self, settle_dt):
        """Number of days between the previous coupon and the currrent step
        in date."""

        # I assume accrued runs to the effective date
        pcd, _ = self.get_pcd(settle_dt)
        eff_dt = max(settle_dt, self.step_in_dt)
        accrued_days = eff_dt - pcd
        return accrued_days

    ###########################################################################

    def accrued_interest(self, settle_dt):
        """Calculate the amount of accrued interest that has accrued from the
        previous cpn date (PCD) to the step_in_dt of the CDS contract."""

        day_count = DayCount(self.accrual_dc_type)
        pcd, _ = self.get_pcd(settle_dt)
        eff_dt = max(settle_dt, self.step_in_dt)
        accrual_factor = day_count.year_frac(pcd, eff_dt)[0]
        # This is always a positive quantity
        accrued_interest = accrual_factor * self.notional * self.running_cpn
        return accrued_interest

    ###########################################################################

    def prot_leg_pv(
        self,
        value_dt,
        issuer_curve,
        contract_recovery_rate=STANDARD_RECOVERY_RATE,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
        prot_method=0,
    ):
        """Calculates the protection leg PV of the CDS by calling into the
        fast NUMBA code that has been defined above."""

        t_eff = (self.step_in_dt - value_dt) / G_DAYS_IN_YEAR

        # An existing contract may have protection which as of now started in the past
        # We need to ensure that looking forward the protection starts immediately
        if t_eff < 0.0:
            t_eff = 0.0

        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEAR

        libor_curve = issuer_curve.libor_curve

        v = _prot_leg_pv_numba(
            t_eff,
            t_mat,
            libor_curve._times,
            libor_curve._dfs,
            issuer_curve._times,
            issuer_curve._qs,
            contract_recovery_rate,
            num_steps_per_year,
            prot_method,
        )

        return v * self.notional

    ###########################################################################

    def get_pcd(self, value_dt):
        """Get the previous coupon date before the value date"""
        pcd = self.accrual_start_dts[0]
        start_index = 0

        # I turned off this warning for CDS Index options
        #        if value_dt < pcd:
        #            raise FinError("Value date before start of first accrual period.")

        for acc_start_dt in self.accrual_start_dts[1:]:
            if value_dt < acc_start_dt:
                break
            else:
                pcd = acc_start_dt
                start_index = start_index + 1

        #        print("Value date:", value_dt, "PCD:", pcd, "Start Index", start_index)

        return pcd, start_index

    ###########################################################################

    def rpv01(self, value_dt, issuer_curve, pv01_method=0):
        """The risky_pv01 is the present value of a risky one dollar paid on
        the premium leg of a CDS contract."""

        libor_curve = issuer_curve.libor_curve

        payment_times = []
        for date in self.payment_dts:
            t = (date - value_dt) / G_DAYS_IN_YEAR

            if t > 0.0:
                payment_times.append(t)

        # this is the part of the cpn accrued from the previous cpn date to now

        pcd, start_index = self.get_pcd(value_dt)

        day_count = DayCount(self.accrual_dc_type)

        eff_dt = self.step_in_dt

        # A contract whose protection started in the past starts today
        eff_dt = max(eff_dt, value_dt)

        accrual_factor_pcd_to_now = day_count.year_frac(pcd, eff_dt)[0]

        #    print("accrued pcd to now", accrual_factor_pcd_to_now)

        year_fracs = self.accrual_factors[start_index:]

        t_eff = (eff_dt - value_dt) / G_DAYS_IN_YEAR

        #    print("riskyPV01 t_eff", t_eff)

        value_rpv01 = _risky_pv01_numba(
            t_eff,
            accrual_factor_pcd_to_now,
            np.array(payment_times),
            np.array(year_fracs),
            libor_curve._times,
            libor_curve._dfs,
            issuer_curve._times,
            issuer_curve._qs,
            pv01_method,
        )

        dirty_rpv01 = value_rpv01[DIRTY]
        clean_rpv01 = value_rpv01[CLEAN]

        return (dirty_rpv01, clean_rpv01)

    ###########################################################################

    def premium_leg_pv(self, value_dt, issuer_curve, pv01_method=0):
        """Value of the premium leg of a CDS."""

        dirty_rpv01 = self.rpv01(value_dt, issuer_curve, pv01_method)[DIRTY]
        v = dirty_rpv01 * self.notional * self.running_cpn
        return v

    ###########################################################################

    def par_spread(
        self,
        value_dt,
        issuer_curve,
        contract_recovery_rate=STANDARD_RECOVERY_RATE,
        num_steps_per_year=GLOB_NUM_STEPS_PER_YEAR,
        pv01_method=0,
        prot_method=0,
    ):
        """Breakeven CDS cpn that would make the value of the CDS contract
        equal to zero."""

        clean_rpv01 = self.rpv01(value_dt, issuer_curve, pv01_method)[CLEAN]

        #        print("cleanRPV01", clean_rpv01)

        prot = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            contract_recovery_rate,
            num_steps_per_year,
            prot_method,
        )

        #        print("prot", prot)

        # By convention this is calculated using the clean RPV01
        spd = prot / clean_rpv01 / self.notional
        return spd

    ###########################################################################

    def value_fast_approx(
        self,
        value_dt,
        flat_cont_interest_rate,
        flat_cds_curve_spread,
        curve_recovery=STANDARD_RECOVERY_RATE,
        contract_recovery_rate=STANDARD_RECOVERY_RATE,
        bump_size=ONE_BP,
        recovery_bump_size=ONE_PCT,
    ):
        """Fast approximate CDS valuation using flat hazard and discount curves.

        The hazard rate is implied from the flat CDS spread via the credit
        triangle ``h = spread / (1 - curve_recovery)`` and both survival
        and discounting are treated as flat exponentials. The premium leg
        ignores the actual coupon schedule; accrual-on-default is captured
        only heuristically through the 365/360 clean-RPV01 adjustment. All
        PVs are in currency units of the contract notional.

        Args:
            value_dt: Valuation date (a ``Date``). Scalar only.
            flat_cont_interest_rate: Flat continuously compounded discount
                rate. Scalar or array-like.
            flat_cds_curve_spread: Flat CDS par spread (decimal, e.g. 0.01
                for 100 bp). Must be non-negative. Scalar or array-like.
            curve_recovery: Recovery rate used to imply the hazard rate
                from the spread. Must lie in ``[0, 1)``. Scalar or
                array-like.
            contract_recovery_rate: Contractual recovery rate used on the
                protection leg. Must lie in ``[0, 1]``. Scalar or
                array-like.
            bump_size: Absolute bump applied to the spread and to the
                interest rate for the finite-difference sensitivities.
                Scalar or array-like, strictly positive.
            recovery_bump_size: Absolute bump applied to the contractual
                recovery rate for the recovery sensitivity. Scalar or
                array-like, strictly positive.

        Returns:
            A tuple ``(full_pv, clean_pv, spread_dv01, ir_dv01, recovery01)``:
                full_pv: Clean PV plus accrued interest.
                clean_pv: PV excluding accrued interest.
                credit01: PV change for a 1 bp increase in the CDS spread.
                ir01: PV change for a 1 bp increase in the interest rate.
                recovery01: PV change for a 1 percentage-point increase in
                    the contractual recovery.

            Elements are floats if every input was scalar, otherwise ndarrays
            with the broadcast shape of the inputs. The sensitivities are
            one-sided (upward) finite differences, rescaled to the stated
            units whatever bump sizes are supplied. Recovery01 bumps the
            recovery assumption everywhere it enters: the hazard rate is
            re-implied from the unchanged quoted spread under the bumped
            curve recovery, and the contractual loss uses the bumped
            contract recovery. It is therefore small and vanishes for an
            at-market contract.

        Raises:
            FinError: If any input is invalid, the inputs cannot be
                broadcast together, or the contract has matured.
        """
        if not isinstance(value_dt, Date):
            raise FinError(
                "Valuation date must be a Date and not " + str(value_dt)
            )

        rate = np.asarray(flat_cont_interest_rate, dtype=float)
        spread = np.asarray(flat_cds_curve_spread, dtype=float)
        curve_rec = np.asarray(curve_recovery, dtype=float)
        contract_rec = np.asarray(contract_recovery_rate, dtype=float)
        bump = np.asarray(bump_size, dtype=float)
        rec_bump = np.asarray(recovery_bump_size, dtype=float)

        inputs = (rate, spread, curve_rec, contract_rec, bump, rec_bump)
        scalar_output = all(arr.ndim == 0 for arr in inputs)

        try:
            out_shape = np.broadcast_shapes(*(arr.shape for arr in inputs))
        except ValueError as exc:
            raise FinError(
                "Inputs could not be broadcast together: " + str(exc))

        if np.any(curve_rec < 0.0) or np.any(curve_rec >= 1.0):
            raise FinError("Curve recovery must lie in [0.0, 1.0).")

        if np.any(contract_rec < 0.0) or np.any(contract_rec > 1.0):
            raise FinError("Contract recovery rate must lie in [0.0, 1.0].")

        if np.any(bump <= 0.0):
            raise FinError("Bump size must be positive.")

        if np.any(rec_bump <= 0.0):
            raise FinError("Recovery bump size must be positive.")

        if np.any(contract_rec + rec_bump > 1.0):
            raise FinError("Bumped contract recovery rate cannot exceed 1.0.")

        if np.any(curve_rec + rec_bump >= 1.0):
            raise FinError("Bumped curve recovery must be less than 1.0.")

        if np.any(spread < 0.0):
            raise FinError("CDS spread cannot be negative.")

        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEAR

        if t_mat <= 0.0:
            raise FinError(
                "Contract has matured: maturity date is on or before the "
                "valuation date."
            )

        t_eff = max((self.step_in_dt - value_dt) / G_DAYS_IN_YEAR, 0.0)

        if t_mat <= t_eff:
            raise FinError(
                "Maturity date must be after the effective protection "
                "start date."
            )

        if self.long_protect:
            direction = 1.0
        else:
            direction = -1.0

        accrued = self.accrued_interest(value_dt) # positive
        delta = accrued / self.notional / self.running_cpn

        horizon = t_mat - t_eff

        def _dirty_and_clean_pv(spread_, rate_, curve_recovery_, contract_recovery_):
            """Dirty and clean PV under flat spread/rate/recovery inputs.

            All arguments broadcast; returns arrays of the broadcast shape.
            """
            hazard_rate = KAPPA * spread_ / (1.0 - curve_recovery_)
            decay_rate = rate_ + hazard_rate

            # Risky annuity: integral of exp(-decay_rate * t) over
            # [t_eff, t_mat]. Written with expm1 so it is stable as
            # decay_rate -> 0, where it tends to the horizon length. The
            # zero-decay elements are masked out of the division and
            # replaced with the exact limit.
            x = decay_rate * horizon
            is_zero = x == 0.0
            safe_decay = np.where(is_zero, 1.0, decay_rate)
            annuity_per_df = np.where(
                is_zero,
                horizon,
                -np.expm1(-x) / safe_decay,
            )

            risky_annuity = np.exp(-decay_rate * t_eff) * annuity_per_df

            rpv01_dirty = risky_annuity * KAPPA + delta
            rpv01_clean = risky_annuity * KAPPA

            protection_pv = (
                hazard_rate
                * (1.0 - contract_recovery_)
                * risky_annuity
                * self.notional
            )

            premium_pv_dirty = self.running_cpn * rpv01_dirty * self.notional
            premium_pv_clean = self.running_cpn * rpv01_clean * self.notional

            dirty_pv = direction * (protection_pv - premium_pv_dirty)
            clean_pv = direction * (protection_pv - premium_pv_clean)

            return dirty_pv, clean_pv, rpv01_dirty, rpv01_clean

        dirty_pv, clean_pv, rpv01_dirty, rpv01_clean = _dirty_and_clean_pv(
            spread, rate, curve_rec, contract_rec
        )

        dirty_pv_spread_up, _, _, _ = _dirty_and_clean_pv(
            spread + bump, rate, curve_rec, contract_rec
        )

        dirty_pv_rate_up, _ , _, _  = _dirty_and_clean_pv(
            spread, rate + bump, curve_rec, contract_rec
        )

        # The recovery bump moves the recovery assumption everywhere it
        # enters: the hazard rate is re-implied from the (unchanged) quoted
        # spread under the bumped curve recovery, and the contractual loss
        # uses the bumped contract recovery. Bumping the contractual
        # recovery alone while holding the implied hazard fixed would mix
        # two inconsistent recovery assumptions in one valuation.
        dirty_pv_recovery_up, _ , _, _ = _dirty_and_clean_pv(
            spread, rate, curve_rec + rec_bump, contract_rec + rec_bump
        )

        # One-sided differences, rescaled to the conventional units so the
        # reported figures stay "per 1 bp" / "per 1 pp" for any bump size.
        spread_dv01 = (dirty_pv_spread_up - dirty_pv) * (ONE_BP / bump)
        ir_dv01 = (dirty_pv_rate_up - dirty_pv) * (ONE_BP / bump)
        rec_dv01 = (dirty_pv_recovery_up - dirty_pv) * (ONE_PCT / rec_bump)

        results = (dirty_pv, clean_pv, rpv01_dirty, rpv01_clean, spread_dv01, ir_dv01, rec_dv01)

        if scalar_output:
            return tuple(float(arr) for arr in results)

        return tuple(
            np.broadcast_to(arr, out_shape).copy() for arr in results
        )

    ###########################################################################

    def print_payments(self, value_dt, issuer_curve):
        """We only print payments after the current valuation date"""
        num_flows = len(self.payment_dts)

        print(
            "PAYMENT_DT      YEAR_FRAC      PAYMENT           DF       SURV_PROB      NPV"
        )

        for it in range(0, num_flows):
            dt = self.payment_dts[it]

            if dt > value_dt:
                acc_factor = self.accrual_factors[it]
                flow = self.flows[it]
                z = issuer_curve.df(dt)
                q = issuer_curve.survival_prob(dt)
                print(
                    "%15s %10.6f %12.2f %12.6f %12.6f %12.2f"
                    % (dt, acc_factor, flow, z, q, flow * z * q)
                )

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all the calculated
        cash flows"""
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("STEP_IN_DATE", self.step_in_dt)
        s += label_to_string("MATURITY", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("LONG_PROT", self.long_protect)
        s += label_to_string("COUPON", self.running_cpn * 10000, "bp\n")
        s += label_to_string("DAY_COUNT", self.accrual_dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS_DAY_ADJUST", self.bd_type)
        s += label_to_string("DATE_GEN_RULE", self.dg_type)

        header = "PAYMENT_DT, YEAR_FRAC, ACCRUAL_START, ACCRUAL_END, PAYMENT"
        value_table = [
            self.payment_dts,
            self.accrual_factors,
            self.accrual_start_dts,
            self.accrual_end_dts,
            self.flows,
        ]
        precision = "12.6f"

        s += table_to_string(header, value_table, precision)

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
