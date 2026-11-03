##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

from copy import deepcopy

import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCountTypes, DayCount
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.global_vars import ONE_MILLION
from ...utils.helpers import label_to_string, table_to_string
from ...utils.global_vars import DIRTY

from ...market.curves.cds_curve import CDSCurve
from .cds import CDS, GLOB_NUM_STEPS_PER_YEAR

USE_FLAT_HAZARD_RATE_INTEGRAL = True
STANDARD_RECOVERY_RATE = 0.40
ONE_BP = 0.0001
ONE_PCT = 0.01

# Premium accrues on an Act/360 basis while time is measured in
# calendar years, hence the 365/360 adjustment to the clean RPV01.
KAPPA = 365.0 / 360.0

########################################################################################
# TODO: Perform protection leg pv analytically using fact that hazard rate and
#       interest rates are flat between their combined node points. Right now I
#       do not find the protection leg PV calculations to be a bottleneck,
#       especially given the speedup benefits of using NUMBA.
########################################################################################


class CDSFixedRecovery(CDS):
    """Fixed-recovery credit default swap."""

    def __init__(
        self,
        step_in_dt: Date,
        maturity_dt_or_tenor: Date | str,
        running_cpn: float,
        fixed_recovery: float,
        notional: float = ONE_MILLION,
        long_protect: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes | list | tuple = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):

        if not 0.0 <= fixed_recovery <= 1.0:
            raise FinError("Fixed recovery must lie between 0 and 1.")

        super().__init__(
            step_in_dt,
            maturity_dt_or_tenor,
            running_cpn,
            notional,
            long_protect,
            freq_type,
            accrual_dc_type,
            cal_type,
            bd_type,
            dg_type,
        )

        self.fixed_recovery = fixed_recovery

    ###########################################################################

    def prot_leg_pv(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float = None,
        num_steps_per_year: int = GLOB_NUM_STEPS_PER_YEAR,
        prot_method=0,
    ):

        if recovery_rate is not None and recovery_rate != self.fixed_recovery:
            raise FinError("Recovery rate differs from contractual fixed recovery.")

        return super().prot_leg_pv(
            value_dt,
            issuer_curve,
            self.fixed_recovery,
            num_steps_per_year,
            prot_method,
        )

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float = None,
        pv01_method=0,
        prot_method=0,
        num_steps_per_year: int = GLOB_NUM_STEPS_PER_YEAR,
    ):

        if recovery_rate is not None and recovery_rate != self.fixed_recovery:
            raise FinError("Recovery rate differs from contractual fixed recovery.")

        return super().value(
            value_dt,
            issuer_curve,
            self.fixed_recovery,
            pv01_method,
            prot_method,
            num_steps_per_year,
        )

    ###########################################################################

    def recovery_dv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        pv01_method: int = 0,
        prot_method: int = 0,
        num_steps_per_year: int = GLOB_NUM_STEPS_PER_YEAR,
    ):
        """PV change for a 1% increase in the issuer curve recovery assumption.

        The contractual fixed recovery, market CDS quotes, and discount
        curve are held unchanged.

        Returns:
            PV(curve recovery + 0.01, fixed recovery)
            - PV(curve recovery, fixed recovery)
        """

        bump = ONE_PCT  # 0.01

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

        value_base = self.value(
            value_dt,
            issuer_curve,
            pv01_method=pv01_method,
            prot_method=prot_method,
            num_steps_per_year=num_steps_per_year,
        )

        bumped_issuer_curve = deepcopy(issuer_curve)
        bumped_issuer_curve.recovery_rate += ONE_PCT
        bumped_issuer_curve.build_curve()

        value_bumped = self.value(
            value_dt,
            bumped_issuer_curve,
            pv01_method=pv01_method,
            prot_method=prot_method,
            num_steps_per_year=num_steps_per_year,
        )

        return value_bumped[DIRTY] - value_base[DIRTY]

    ###########################################################################

    def value_fast_approx(
        self,
        value_dt: Date,
        flat_cont_interest_rate: float,
        flat_cds_curve_spread,
        curve_recovery: float = STANDARD_RECOVERY_RATE,
        bump_size: float = ONE_BP,
        recovery_bump_size: float = ONE_PCT,
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
            raise FinError("Valuation date must be a Date and not " + str(value_dt))

        rate = np.asarray(flat_cont_interest_rate, dtype=float)
        spread = np.asarray(flat_cds_curve_spread, dtype=float)
        curve_rec = np.asarray(curve_recovery, dtype=float)
        contract_rec = np.asarray(self.fixed_recovery, dtype=float)
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

        if np.any(bump <= 0.0):
            raise FinError("Bump size must be positive.")

        if np.any(rec_bump <= 0.0):
            raise FinError("Recovery bump size must be positive.")

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
                "Maturity date must be after the effective protection " "start date."
            )

        if self.long_protect:
            direction = 1.0
        else:
            direction = -1.0

        pcd, _ = self.get_pcd(value_dt)
        eff_dt = max(value_dt, self.step_in_dt)

        day_count = DayCount(self.accrual_dc_type)
        delta = day_count.year_frac(pcd, eff_dt)[0]

        horizon = t_mat - t_eff

        def _dirty_and_clean_pv(
            spread_, rate_: float, curve_recovery_: float, contract_recovery_: float
        ):
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
            annuity_per_df = np.where(is_zero, horizon, -np.expm1(-x) / safe_decay)

            risky_annuity = np.exp(-decay_rate * t_eff) * annuity_per_df

            rpv01_dirty = risky_annuity * KAPPA + delta
            rpv01_clean = risky_annuity * KAPPA

            protection_pv = (
                hazard_rate * (1.0 - contract_recovery_) * risky_annuity * self.notional
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

        dirty_pv_rate_up, _, _, _ = _dirty_and_clean_pv(
            spread, rate + bump, curve_rec, contract_rec
        )

        # Bump the recovery assumption used to imply the hazard rate.
        # The contractual fixed recovery remains unchanged.
        dirty_pv_recovery_up, _, _, _ = _dirty_and_clean_pv(
            spread,
            rate,
            curve_rec + rec_bump,
            contract_rec,
        )

        # One-sided differences, rescaled to the conventional units so the
        # reported figures stay "per 1 bp" / "per 1 pp" for any bump size.
        spread_dv01 = (dirty_pv_spread_up - dirty_pv) * (ONE_BP / bump)
        ir_dv01 = (dirty_pv_rate_up - dirty_pv) * (ONE_BP / bump)
        rec_dv01 = (dirty_pv_recovery_up - dirty_pv) * (ONE_PCT / rec_bump)

        results = (
            dirty_pv,
            clean_pv,
            rpv01_dirty,
            rpv01_clean,
            spread_dv01,
            ir_dv01,
            rec_dv01,
        )

        if scalar_output:
            return tuple(float(arr) for arr in results)

        return tuple(np.broadcast_to(arr, out_shape).copy() for arr in results)

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("STEP-IN DATE", self.step_in_dt)
        s += label_to_string("MATURITY", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("FIXED RECOVERY", self.fixed_recovery*100, "%\n")
        s += label_to_string("LONG PROTECTION", self.long_protect)
        s += label_to_string("RUNNING SPREAD", self.running_cpn * 10000, " bp\n")
        s += label_to_string("DAY COUNT", self.accrual_dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST", self.bd_type)
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
