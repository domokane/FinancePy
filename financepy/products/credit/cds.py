##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

from ...utils.date import Date
from ...utils.valuation_log import log_amount, log_percentage, log_bps, log_float, log_int, log_date
from ...utils.valuation_log import log_value, log_table, log_row, log_section, suppress_logging, logging_enabled
from ...market.curves.cds_curve import CDSCurve
from ...utils.helpers import check_argument_types
from ...utils.check_values import check_curve_dt
from ...utils.helpers import label_to_string, table_to_string
from ...utils.global_vars import ONE_MILLION
from ...utils.global_vars import CLEAN, DIRTY
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.global_types import CDSModelTypes
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.calendar import Calendar, CalendarTypes
from ...utils.error import FinError
from typing import Union

from ...models.cds_fast_model import cds_fast_value
from ...models.cds_okane_turnbull_model import prot_leg_pv_ot_numba
from ...models.cds_okane_turnbull_model import rpv01_ot_numba
from ...models.cds_isda_model import prot_leg_pv_isda_numba
from ...models.cds_isda_model import rpv01_isda_numba

from copy import deepcopy
import numpy as np


STANDARD_RECOVERY_RATE = 0.40
ONE_BP = 0.0001
ONE_PCT = 0.01

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
        cal_type: CalendarTypes | list | tuple = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
        accrue_from_step_in_dt: bool = False
    ) -> None:
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

        if step_in_dt >= maturity_dt:
            raise FinError(
                "Step-in date must precede maturity date."
            )

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
        self.accrue_from_step_in_dt = accrue_from_step_in_dt

        self._rpv01_cache_dt = None
        self._rpv01_cache = None

        self._generate_adjusted_cds_payment_dts()
        self._calc_flows()

    ####################################################################################

    def _generate_adjusted_cds_payment_dts(self):
        """Generate CDS payment dates which have been holiday adjusted."""

        self._rpv01_cache_dt = None
        self._rpv01_cache = None

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
        self.accrual_end_dts = [date.add_days(-1) for date in self.accrual_start_dts[1:]]

        # Final accrual end date is the maturity date
        self.accrual_end_dts.append(self.maturity_dt)

    ###########################################################################

    def _calc_flows(self):
        """Calculate cash flow amounts on premium leg."""

        self._rpv01_cache_dt = None
        self._rpv01_cache = None

        day_count = DayCount(self.accrual_dc_type)

        self.accrual_factors = []
        self.flows = []

        # self.accrual_factors.append(0.0)
        # self.flows.append(0.0)

        for t0, t1 in zip(self.accrual_start_dts, self.accrual_end_dts):

            if self.accrue_from_step_in_dt:
                t0 = max(t0, self.step_in_dt)

            accrual_factor = day_count.year_frac(
                t0, t1.add_days(1)
            )[0]

            flow = (
                accrual_factor
                * self.running_cpn
                * self.notional
            )

            self.accrual_factors.append(accrual_factor)
            self.flows.append(flow)

    ###########################################################################

    def _prepare_rpv01(self, value_dt: Date):
        """Prepare valuation-date-dependent premium leg inputs.

        These inputs are independent of the discount curve,
        survival curve, recovery rate and pricing model.
        """

        if (self._rpv01_cache_dt is not None and self._rpv01_cache_dt == value_dt):
            return self._rpv01_cache

        step_in_dt = self._get_step_in_dt(value_dt)

        pcd, start_index = self.get_pcd(step_in_dt)

        payment_dts = []
        payment_times = []
        accrual_start_times = []
        accrual_end_times = []
        year_fracs = []

        for i in range(start_index, len(self.payment_dts)):

            dt = self.payment_dts[i]

            if dt > value_dt:

                start_dt = self.accrual_start_dts[i]

                if self.accrue_from_step_in_dt:
                    start_dt = max(start_dt, self.step_in_dt)

                payment_dts.append(dt)

                payment_times.append(
                    (dt - value_dt) / G_DAYS_IN_YEAR
                )

                accrual_start_times.append(
                    (start_dt - value_dt) / G_DAYS_IN_YEAR
                )

                accrual_end_times.append(
                    (self.accrual_end_dts[i] - value_dt + 1)
                    / G_DAYS_IN_YEAR
                )

                year_fracs.append(self.accrual_factors[i])

        accrual_start_dt = pcd

        if self.accrue_from_step_in_dt:
            accrual_start_dt = max(pcd, self.step_in_dt)

        day_count = DayCount(self.accrual_dc_type)

        accrued_fraction = day_count.year_frac(
            accrual_start_dt,
            step_in_dt,
        )[0]

        self._rpv01_cache = (
            step_in_dt,
            pcd,
            accrual_start_dt,
            accrued_fraction,
            payment_dts,
            np.array(payment_times),
            np.array(accrual_start_times),
            np.array(accrual_end_times),
            np.array(year_fracs),
        )

        self._rpv01_cache_dt = value_dt

        return self._rpv01_cache

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Valuation of a CDS contract on a specific valuation date
        given an issuer curve and a contract recovery rate."""

        check_curve_dt(value_dt, issuer_curve)

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        log_section("Valuation Details")

        rpv01 = self.rpv01(value_dt,
                           issuer_curve,
                           settle_dt,
                           pay_accrued_on_default)

        dirty_rpv01 = rpv01[DIRTY]
        clean_rpv01 = rpv01[CLEAN]

        log_float("Dirty RPV01", dirty_rpv01, 6, 2)
        log_float("Clean RPV01", clean_rpv01, 6, 2)

        prot_pv = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            recovery_rate,
        )

        log_percentage("Recovery Rate", recovery_rate, 2)
        log_amount("Protection PV", prot_pv, 3)

        omega = +1 if self.long_protect else -1

        dirty_pv = omega * (
            prot_pv
            - self.running_cpn * dirty_rpv01 * self.notional
        )

        clean_pv = omega * (
            prot_pv
            - self.running_cpn * clean_rpv01 * self.notional
        )

        log_amount("Dirty PV", dirty_pv, 2)
        log_amount("Clean PV", clean_pv, 2)

        return (dirty_pv, clean_pv)

    ###########################################################################

    def spread_dv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Calculation of the change in the value of the CDS contract for a
        one basis point change in the level of the CDS curve."""

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        with suppress_logging():

            v0 = self.value(
                value_dt,
                issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
            )

        # we create a deep copy to avoid state issues
        bumped_issuer_curve = deepcopy(issuer_curve)
        for cds in bumped_issuer_curve.cds_contracts:
            cds.running_cpn += ONE_BP

        bumped_issuer_curve.build_curve()

        with suppress_logging():

            v1 = self.value(
                value_dt,
                bumped_issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
            )

        credit_dv01 = v1[DIRTY] - v0[DIRTY]

        log_amount("Credit DV01", credit_dv01)

        return credit_dv01

    ###########################################################################

    def ir_dv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Calculation of the interest DV01 based on a simple bump of
        the discount factors and reconstruction of the CDS curve."""

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        with suppress_logging():

            v0 = self.value(
                value_dt,
                issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
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

        bumped_issuer_curve.libor_curve = original_discount_curve.bump_parallel(ONE_BP)
        bumped_issuer_curve.build_curve()

        with suppress_logging():

            v1 = self.value(
                value_dt,
                bumped_issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
            )

        interest_dv01 = v1[DIRTY] - v0[DIRTY]

        log_amount("Interest DV01", interest_dv01, 2)

        return interest_dv01

    ###########################################################################

    def recovery_dv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """PV change when contract and curve recovery both increase by 1%.

        The market CDS calibration quotes and discount curve are held fixed.

        Returns:
            PV(curve recovery + 0.01, contract recovery + 0.01)
            - PV(curve recovery, contract recovery)
        """

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        bump = ONE_PCT

        if not 0.0 <= recovery_rate < 1.0:
            raise FinError("Contract recovery rate must lie in [0.0, 1.0).")

        if recovery_rate + bump >= 1.0:
            raise FinError("Bumped contract recovery rate must be less than 1.0.")

        if not hasattr(issuer_curve, "recovery_rate"):
            raise FinError(
                f"Issuer curve type {type(issuer_curve).__name__} "
                "does not expose recovery_rate."
            )

        curve_recovery_rate = issuer_curve.recovery_rate

        if not 0.0 <= curve_recovery_rate < 1.0:
            raise FinError("Issuer curve recovery rate must lie in [0.0, 1.0).")

        if curve_recovery_rate + bump >= 1.0:
            raise FinError("Bumped issuer curve recovery rate must be less than 1.0.")

        # Base valuation:
        # - original survival curve
        # - original contract recovery
        with suppress_logging():

            value_base = self.value(
                value_dt,
                issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
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

        with suppress_logging():

            value_bumped = self.value(
                value_dt,
                bumped_issuer_curve,
                recovery_rate + bump,
                settle_dt,
                pay_accrued_on_default
            )

        recovery_dv01 = value_bumped[DIRTY] - value_base[DIRTY]

        log_amount("Recovery DV01", recovery_dv01, 2)

        return recovery_dv01

    ###########################################################################

    def upfront(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Return clean percentage PV expressed on T+3 settlement date."""

        discount_curve = issuer_curve.libor_curve

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        if settle_dt < discount_curve.anchor_dt:
            raise FinError(
                "Settlement date cannot precede the discount curve date: "
                f"settle_dt={settle_dt}, "
                f"curve.anchor_dt={discount_curve.anchor_dt}"
            )

        with suppress_logging():

            values = self.value(
                value_dt,
                issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
            )

        settlement_df = discount_curve.df(settle_dt)

        log_float("Settlement DF", settlement_df, 6, 3)

        if not np.isfinite(settlement_df) or settlement_df <= 0.0:
            raise FinError(
                "Settlement-date discount factor must be finite and positive: "
                f"df={settlement_df}"
            )

        upfront = values[CLEAN] / settlement_df / self.notional

        log_percentage("Upfront", upfront, 5, 2)

        return upfront

    ###########################################################################

    def cash_settlement_amount(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Return dirty amount paid on the T+3 settlement date."""

        discount_curve = issuer_curve.libor_curve

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        if settle_dt < discount_curve.anchor_dt:
            raise FinError(
                "Settlement date cannot precede the discount curve date: "
                f"settle_dt={settle_dt}, "
                f"curve.anchor_dt={discount_curve.anchor_dt}"
            )

        with suppress_logging():

            values = self.value(
                value_dt,
                issuer_curve,
                recovery_rate,
                settle_dt,
                pay_accrued_on_default
            )

        settlement_df = discount_curve.df(settle_dt)

        log_float("Settlement DF", settlement_df, 6, 3)

        if not np.isfinite(settlement_df) or settlement_df <= 0.0:
            raise FinError(
                "Settlement-date discount factor must be finite and positive: "
                f"df={settlement_df}"
            )

        csa = values[DIRTY] / settlement_df

        log_amount("Cash Settlement Amount", csa, 2)

        return csa

    ###########################################################################

    def clean_price(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):
        """Value of the CDS contract excluding accrued interest."""
        check_curve_dt(value_dt, issuer_curve)

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        rpv01 = self.rpv01(value_dt,
                           issuer_curve,
                           settle_dt,
                           pay_accrued_on_default)

        clean_rpv01 = rpv01[CLEAN]

        prot_pv = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            recovery_rate,
        )

        settlement_df = issuer_curve.libor_curve.df(settle_dt)

        clean_pv = prot_pv - self.running_cpn * clean_rpv01 * self.notional

        clean_price = 100.0 * (1.0 - clean_pv / self.notional / settlement_df)

        log_float("Clean Price", clean_price, level=2)

        return clean_price

    ###########################################################################

    def accrued_days(self, value_dt: Date):
        """Number of accrued premium days from the start of the
        current accrual period to the pricing step-in date."""

        step_in_dt = self._get_step_in_dt(value_dt)

        pcd, _ = self.get_pcd(step_in_dt)

        if self.accrue_from_step_in_dt:
            pcd = max(pcd, self.step_in_dt)

        accrued_days = step_in_dt - pcd

        log_int("Accrued Days", accrued_days, 2)

        return accrued_days

    ###########################################################################

    def accrued_interest(self, value_dt: Date):
        """Calculate the accrued premium rebate amount."""

        accrual_fraction = self._get_accrual_fraction(value_dt)

        accrued_interest = (
            accrual_fraction
            * self.notional
            * self.running_cpn
        )

        log_amount("Accrued Interest", accrued_interest, 2)

        return accrued_interest

    ###########################################################################

    def prot_leg_pv(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float = STANDARD_RECOVERY_RATE
    ):

        cds_model = self._get_cds_model(issuer_curve)

        """Calculates the protection leg PV of the CDS by calling into the
        fast NUMBA code that has been defined above."""
        check_curve_dt(value_dt, issuer_curve)

        eff_dt = self._get_eff_dt(value_dt, cds_model)
        t_eff = (eff_dt - value_dt) / G_DAYS_IN_YEAR
        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEAR

        libor_curve = issuer_curve.libor_curve

        if cds_model == CDSModelTypes.OKANE_TURNBULL:

            v = prot_leg_pv_ot_numba(
                t_eff,
                t_mat,
                libor_curve._times,
                libor_curve._dfs,
                issuer_curve._times,
                issuer_curve._qs,
                recovery_rate,
            )

        elif cds_model == CDSModelTypes.ISDA:

            v = prot_leg_pv_isda_numba(
                t_eff,
                t_mat,
                libor_curve._times,
                libor_curve._dfs,
                issuer_curve._times,
                issuer_curve._qs,
                recovery_rate,
            )

        else:

            raise FinError(f"Unknown CDS model: {cds_model}")

        return v * self.notional

    ###########################################################################

    def get_pcd(self, value_dt: Date):
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

        return pcd, start_index

    ###########################################################################

    def rpv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True,
    ):
        """Calculate the clean and dirty risky PV01 of a CDS."""

        check_curve_dt(value_dt, issuer_curve)

        cds_model = issuer_curve.cds_model

        libor_curve = issuer_curve.libor_curve

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        settlement_df = libor_curve.df(settle_dt)

        if not np.isfinite(settlement_df) or settlement_df <= 0.0:
            raise FinError(
                "Settlement discount factor must be finite and positive."
            )

        log_date("Settlement Date", settle_dt, level=3)
        log_float("Settlement DF", settlement_df, 6, 3)

        # Retrieve the cached numerical schedule.

        (
            step_in_dt,
            pcd,
            accrual_start_dt,
            accrual_factor_pcd_to_now,
            payment_dts,
            payment_times,
            accrual_start_times,
            accrual_end_times,
            year_fracs,
        ) = self._prepare_rpv01(value_dt)

        if len(payment_dts) == 0:
            return np.zeros(2)

        # Determine the effective time required by the kernel.

        eff_dt = self._get_eff_dt(value_dt, cds_model)

        t_eff = (eff_dt - value_dt) / G_DAYS_IN_YEAR

        log_value("CDS Model", cds_model.name, level=3)
        log_value("Previous Coupon Date", pcd, level=3)
        log_value("Accrual Start Date", accrual_start_dt, level=3)
        log_value("Pricing Step-in Date", step_in_dt, level=3)
        log_value("Kernel Effective Date", eff_dt, level=3)

        log_value(
            "Accrual Factor to Step-in Date",
            accrual_factor_pcd_to_now,
            level=3,
        )

        # Calculate the risky PV01 using the selected model.

        if cds_model == CDSModelTypes.OKANE_TURNBULL:

            value_rpv01, table = rpv01_ot_numba(
                t_eff,
                accrual_factor_pcd_to_now,
                payment_times,
                year_fracs,
                libor_curve._times,
                libor_curve._dfs,
                issuer_curve._times,
                issuer_curve._qs,
                accrual_start_times,
                accrual_end_times,
                settlement_df,
                pay_accrued_on_default,
            )

        elif cds_model == CDSModelTypes.ISDA:

            value_rpv01, table = rpv01_isda_numba(
                t_eff,
                accrual_factor_pcd_to_now,
                payment_times,
                year_fracs,
                libor_curve._times,
                libor_curve._dfs,
                issuer_curve._times,
                issuer_curve._qs,
                accrual_start_times,
                accrual_end_times,
                settlement_df,
                pay_accrued_on_default,
            )

        else:

            raise FinError(f"Unknown CDS model: {cds_model}")

        # Want to skip this loop entirely if logging is not on
        if logging_enabled(level=4):

            log_table(
                "Payment Date",
                "T",
                "Q",
                "Z",
                "Acc Factor",
                "Cumulative RPV01",
                level=4,
            )

            for i, row in enumerate(table):

                log_row(
                    payment_dts[i],
                    row[0],
                    row[1],
                    row[2],
                    row[3],
                    row[4],
                )

        return value_rpv01

    ###########################################################################

    def premium_leg_pv(self,
                       value_dt: Date,
                       issuer_curve: CDSCurve,
                       settle_dt: Date | None = None,
                       pay_accrued_on_default: bool = True):
        """Value of the premium leg of a CDS."""

        dirty_rpv01 = self.rpv01(value_dt,
                                 issuer_curve,
                                 settle_dt,
                                 pay_accrued_on_default)[DIRTY]

        v = dirty_rpv01 * self.notional * self.running_cpn
        return v

    ###########################################################################

    def par_spread(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float = STANDARD_RECOVERY_RATE,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True,
        all_values: bool = False
    ):
        """Breakeven CDS cpn that would make the value of the CDS contract
        equal to zero."""

        settle_dt = self.get_settle_dt(value_dt, settle_dt)

        clean_rpv01 = self.rpv01(value_dt,
                                 issuer_curve,
                                 settle_dt,
                                 pay_accrued_on_default)[CLEAN]

        prot_pv = self.prot_leg_pv(
            value_dt,
            issuer_curve,
            recovery_rate,
        )

        if clean_rpv01 <= 0.0:
            raise FinError("Clean RPV01 <= 0.0")

        # By convention this is calculated using the clean RPV01
        spd = prot_pv / clean_rpv01 / self.notional

        log_bps("Par CDS Spread", spd, level=2)

        if all_values:
            return spd, clean_rpv01, prot_pv
        else:
            return spd

    ###########################################################################

    def value_fast_approx(
        self,
        value_dt: Date,
        flat_cont_interest_rate: float,
        flat_cds_curve_spread,
        curve_recovery: float = STANDARD_RECOVERY_RATE,
        recovery_rate: float = STANDARD_RECOVERY_RATE,
        bump_size: float = ONE_BP,
        recovery_bump_size: float = ONE_PCT,
    ):
        """Fast approximate CDS valuation using flat hazard and discount curves.

        The hazard rate is implied from the flat CDS spread via the credit
        triangle ``h = (365/360) * spread / (1 - curve_recovery)`` and both survival
        and discounting are treated as flat exponentials. 

        The premium leg is approximated using continuous premium
        payments with an ACT/360 day-count adjustment. Accrued
        premium on default is not modelled separately.

        The accrued premium rebate is assumed to be exchanged
        on the valuation date rather than the cash settlement date.

        For contracts with accrue_from_step_in_dt=True, accrued
        premium is zero before the step-in date and begins
        accumulating from that date.

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
            A tuple containing:

                dirty_pv: Dirty CDS value in currency units.
                clean_pv: Clean CDS value in currency units.
                rpv01_dirty: Dirty risky PV01 per unit notional.
                rpv01_clean: Clean risky PV01 per unit notional.
                spread_dv01: Spread sensitivity per basis point.
                ir_dv01: Interest rate sensitivity per basis point.
                rec_dv01: Recovery sensitivity per percentage point.

            Elements are floats if all inputs are scalar, otherwise
            ndarrays with the common broadcast shape.

            Sensitivities are calculated using central finite
            differences and rescaled to the conventional risk units.

            Recovery sensitivity bumps both the curve and contractual
            recovery rates while holding the market CDS spread fixed.

        Raises:
            FinError: If any input is invalid, the inputs cannot be
                broadcast together, or the contract has matured.
        """
        if not isinstance(value_dt, Date):
            raise FinError("Valuation date must be a Date and not " + str(value_dt))

        rate = np.asarray(flat_cont_interest_rate, dtype=float)
        spread = np.asarray(flat_cds_curve_spread, dtype=float)
        curve_rec = np.asarray(curve_recovery, dtype=float)
        contract_rec = np.asarray(recovery_rate, dtype=float)
        bump = np.asarray(bump_size, dtype=float)
        rec_bump = np.asarray(recovery_bump_size, dtype=float)

        inputs = (rate, spread, curve_rec, contract_rec, bump, rec_bump)

        scalar_output = all(arr.ndim == 0 for arr in inputs)

        try:
            out_shape = np.broadcast_shapes(*(arr.shape for arr in inputs))
        except ValueError as exc:
            raise FinError("Inputs could not be broadcast together: " + str(exc))

        if np.any(curve_rec < 0.0) or np.any(curve_rec >= 1.0):
            raise FinError("Curve recovery must lie in [0.0, 1.0).")

        if np.any(contract_rec < 0.0) or np.any(contract_rec > 1.0):
            raise FinError("Contract recovery rate must lie in [0.0, 1.0].")

        if np.any(curve_rec - rec_bump < 0.0):
            raise FinError("Bumped curve recovery cannot be negative.")

        if np.any(contract_rec - rec_bump < 0.0):
            raise FinError("Bumped contract recovery cannot be negative.")

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

        for arr in inputs:
            if not np.all(np.isfinite(arr)):
                raise FinError("Inputs must be finite.")

        t_mat = (self.maturity_dt - value_dt) / G_DAYS_IN_YEAR

        if t_mat <= 0.0:
            raise FinError(
                "Contract has matured: maturity date is on or before the "
                "valuation date."
            )

        eff_dt = self._get_eff_dt(value_dt, CDSModelTypes.OKANE_TURNBULL)

        t_eff = (eff_dt - value_dt) / G_DAYS_IN_YEAR

        if t_mat <= t_eff:
            raise FinError(
                "Maturity date must be after the effective protection " "start date."
            )

        if self.long_protect:
            direction = 1.0
        else:
            direction = -1.0

        delta = self._get_accrual_fraction(value_dt)

        out = cds_fast_value(
            self.running_cpn,
            spread,
            t_eff,
            t_mat,
            delta,
            rate,
            curve_rec,
            contract_rec,
            self.notional,
            direction,
            bump,
            rec_bump,
        )

        if scalar_output:
            return tuple(float(x) for x in out)

        return tuple(
            np.broadcast_to(x, out_shape).copy()
            for x in out
        )

    ###########################################################################

    def _get_accrual_fraction(self, value_dt: Date):
        """Calculate the contractual accrued premium fraction."""

        step_in_dt = self._get_step_in_dt(value_dt)

        pcd, _ = self.get_pcd(step_in_dt)

        if self.accrue_from_step_in_dt:
            pcd = max(pcd, self.step_in_dt)

        day_count = DayCount(self.accrual_dc_type)

        accrual_fraction = day_count.year_frac(
            pcd,
            step_in_dt,
        )[0]

        return accrual_fraction

    ###########################################################################

    def get_settle_dt(
        self,
        value_dt: Date,
        settle_dt: Date | None = None,
    ):
        """Determine the cash settlement date."""

        if settle_dt is None:

            calendar = Calendar(self.cal_type)
            settle_dt = calendar.add_business_days(value_dt, 3)

        if settle_dt < value_dt:
            raise FinError("Settlement date cannot precede valuation date.")

        return settle_dt

    ###########################################################################

    def _get_eff_dt(
        self,
        value_dt: Date,
        cds_model: CDSModelTypes = None,
    ):
        """Determine the date argument required by the numerical kernel.
        Both models use the same pricing step-in date, but the ISDA
        model applies the one-day protection adjustment internally.
        """

        step_in_dt = self._get_step_in_dt(value_dt)

        if cds_model == CDSModelTypes.OKANE_TURNBULL:

            return max(value_dt, step_in_dt.add_days(-1))

        if cds_model == CDSModelTypes.ISDA:

            return step_in_dt

        raise FinError(f"Unknown CDS model: {cds_model}")

    ###########################################################################

    def _get_step_in_dt(self, value_dt: Date):
        """Determine the pricing step-in date."""

        if self.step_in_dt > value_dt:
            return self.step_in_dt

        return value_dt.add_days(1)

    ###########################################################################

    def _get_cds_model(
        self,
        issuer_curve: CDSCurve,
    ):
        """Resolve the pricing model and check calibration consistency."""

        cds_model = issuer_curve.cds_model

        return cds_model

    ###########################################################################

    def print_payments(self, value_dt: Date, issuer_curve: CDSCurve):
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
                obs_dt = self.accrual_end_dts[it]
                q = issuer_curve.survival_prob(obs_dt)
                print(
                    "%15s %10.6f %12.2f %12.6f %12.6f %12.2f"
                    % (dt, acc_factor, flow, z, q, flow * z * q)
                )

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self.step_in_dt)
        s += label_to_string("TERMINATION DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("LONG PROTECTION", self.long_protect)
        s += label_to_string("RUNNING SPREAD", self.running_cpn * 10000, "bp\n")
        s += label_to_string("DAY COUNT", self.accrual_dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY RULE", self.bd_type)
        s += label_to_string("DATE GEN RULE", self.dg_type)

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
