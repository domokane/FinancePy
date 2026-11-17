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
from ...utils.global_types import CDSModelTypes
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.global_vars import ONE_MILLION
from ...utils.helpers import label_to_string, table_to_string
from ...utils.global_vars import DIRTY
from ...utils.valuation_log import log_value, log_table, log_row
from ...models.cds_fast_model import cds_fast_value

from ...market.curves.cds_curve import CDSCurve
from .cds import CDS, ONE_PCT, ONE_BP
from .cds import STANDARD_RECOVERY_RATE


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
    ):

        if recovery_rate is not None and recovery_rate != self.fixed_recovery:
            raise FinError("Recovery rate differs from contractual fixed recovery.")

        return super().prot_leg_pv(
            value_dt,
            issuer_curve,
            self.fixed_recovery,
        )

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
        recovery_rate: float = None,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True
    ):

        if recovery_rate is not None and recovery_rate != self.fixed_recovery:
            raise FinError("Recovery rate differs from contractual fixed recovery.")

        return super().value(
            value_dt,
            issuer_curve,
            self.fixed_recovery,
            settle_dt,
            pay_accrued_on_default
        )

    ###########################################################################

    def recovery_dv01(
        self,
        value_dt: Date,
        issuer_curve: CDSCurve,
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
        )

        bumped_issuer_curve = deepcopy(issuer_curve)
        bumped_issuer_curve.recovery_rate += ONE_PCT
        bumped_issuer_curve.build_curve()

        value_bumped = self.value(
            value_dt,
            bumped_issuer_curve,
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
            True
        )

        if scalar_output:
            return tuple(float(x) for x in out)

        return tuple(
            np.broadcast_to(x, out_shape).copy()
            for x in out
        )

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
