##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.error import FinError
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.interpolator import _uinterpolate, InterpTypes
from ...utils.helpers import input_time, table_to_string, times_from_dates
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.helpers import check_argument_types, _func_name
from ...utils.helpers import label_to_string
from ...utils.global_vars import CLEAN
from ...utils.global_types import CDSModelTypes

from ...utils.valuation_log import log_value, log_table, log_row, log_section
# from numba import njit, float64
from ...utils.valuation_log import suppress_logging

########################################################################################


def _f(q, *args):
    """Function that returns zero when the survival probability that gives a
    zero value of the CDS has been determined."""

    curve = args[0]
    value_dt = args[1]
    cds = args[2]
    recovery_rate = args[3]
    settle_dt = args[4]
    pay_accrued_on_default = args[5]

    curve.set_last_q(q)

    # This is important - we calibrate a curve that makes the clean PV of the
    # CDS equal to zero and so we select the second element of the value tuple

    with suppress_logging():
        v = cds.value(value_dt,
                      curve,
                      recovery_rate,
                      settle_dt,
                      pay_accrued_on_default)

    return v[CLEAN]


########################################################################################


class CDSCurve:
    """Generate a survival probability curve implied by the value of CDS
    contracts given a Ibor curve and an assumed recovery rate. The recovery
    rate corresponds to the seniority of the debt for these CDS. A scheme for
    the interpolation of the survival probabilities is also required."""

    def __init__(
        self,
        anchor_dt: Date,
        cds_contracts: list,
        libor_curve: DiscountCurve,
        recovery_rate: float,
        dc_type: DayCountTypes = DayCountTypes.ACT_365F,
        interp_method: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        cds_model: CDSModelTypes = CDSModelTypes.OKANE_TURNBULL,
        settle_dt: Date | None = None,
        pay_accrued_on_default: bool = True,
    ) -> None:
        """Construct a credit curve from a sequence of maturity-ordered CDS
        contracts and a Ibor curve using the same recovery rate and the
        same interpolation method."""

        check_argument_types(getattr(self, _func_name(), None), locals())

        if anchor_dt != libor_curve.anchor_dt:
            raise FinError(
                f"Anchor date mismatch: Issuer curve = {anchor_dt}, "
                f"Ibor curve = {libor_curve.anchor_dt}."
            )

        self.anchor_dt = anchor_dt
        self.cds_contracts = cds_contracts
        self.recovery_rate = recovery_rate
        self.libor_curve = libor_curve
        self.interp_method = interp_method
        self.dc_type = dc_type
        self.built_ok = False

        self.cds_model = cds_model
        self.settle_dt = settle_dt
        self.pay_accrued_on_default = pay_accrued_on_default

        self._times = []
        self._qs = []

        if len(self.cds_contracts) > 0:

            with suppress_logging():

                self.build_curve()

        else:
            pass  # In some cases we allow None to be passed

        return

    ###########################################################################

    def set_times(self, times: np.array):
        """Set the times vector"""
        self._times = times

    ###########################################################################

    def set_qs(self, qvector: np.array):
        """Set the survival probability curve"""
        self._qs = qvector

    ###########################################################################

    def set_q(self, index, q):
        """Set the survival probability at a specific index."""

        n_points = len(self._qs)

        if index < 0 or index >= n_points:
            raise IndexError("Index out of bounds")

        self._qs[index] = q

    ###########################################################################

    def set_last_q(self, q):
        """Set the survival probability factor at last index."""

        n_points = len(self._qs)
        self._qs[n_points - 1] = q

    ###########################################################################

    def _validate(self, cds_contracts):
        """Ensure that contracts are in increasing maturity."""

        if len(cds_contracts) == 0:
            raise FinError("No CDS contracts have been supplied.")

        maturity_dt = cds_contracts[0].maturity_dt

        for cds in cds_contracts[1:]:
            if cds.maturity_dt <= maturity_dt:
                raise FinError("CDS contracts not in increasing maturity.")

            maturity_dt = cds.maturity_dt

    ###########################################################################

    def survival_prob(self, dt):
        """Extract the survival probability to date dt. This function
        supports vectorisation."""

        if isinstance(dt, Date):
            t = times_from_dates(self.anchor_dt, dt, self.dc_type)
        elif isinstance(dt, list):
            if len(dt) > 0 and isinstance(dt[0], Date):
                t = times_from_dates(self.anchor_dt, dt, self.dc_type)
            else:
                t = np.array(dt, dtype=float)
        elif isinstance(dt, (float, int)):
            t = float(dt)
        elif isinstance(dt, np.ndarray):
            t = dt
        else:
            raise FinError("Unknown time type")

        if np.any(t < 0.0):
            raise FinError("Survival Date before curve anchor date")

        if isinstance(t, np.ndarray):
            n = len(t)
            qs = np.zeros(n)
            for i in range(0, n):
                qs[i] = _uinterpolate(t[i],
                                      self._times,
                                      self._qs,
                                      self.interp_method.value)
            return qs
        elif np.isscalar(t):
            q = _uinterpolate(float(t),
                              self._times,
                              self._qs,
                              self.interp_method.value)
            return q

        raise FinError("Unknown time type")

    ###########################################################################

    def df(self, dt):
        """Extract the discount factor from the underlying Ibor curve. This
        function supports vectorisation."""

        if isinstance(dt, Date):
            return self.libor_curve.df(dt)
        elif isinstance(dt, list):
            if len(dt) > 0 and isinstance(dt[0], Date):
                return self.libor_curve.df(dt)
            else:
                return self.libor_curve.df_t(np.array(dt, dtype=float))
        elif isinstance(dt, (float, int, np.ndarray)):
            return self.libor_curve.df_t(dt)
        else:
            raise FinError(f"Unknown date or time type: {type(dt)}")

    ###########################################################################

    def build_curve(self):
        """Construct the CDS survival curve from a set of CDS contracts"""

        self._validate(self.cds_contracts)
        num_times = len(self.cds_contracts)

        # we size the vectors to include time zero
        self._times = np.array([0.0])
        self._qs = np.array([1.0])

        for i in range(0, num_times):

            maturity_dt = self.cds_contracts[i].maturity_dt

            settle_dt = self.cds_contracts[i].get_settle_dt(self.anchor_dt)

            argtuple = (
                self,
                self.anchor_dt,
                self.cds_contracts[i],
                self.recovery_rate,
                settle_dt,
                self.pay_accrued_on_default
            )

            t_mat = times_from_dates(self.anchor_dt, maturity_dt, self.dc_type)

            # Estimate the initial hazard rate using the credit triangle.

            kappa = 365.0 / 360.0

            s_par = self.cds_contracts[i].running_cpn

            hazard_guess = kappa * s_par / (1.0 - self.recovery_rate)

            # Extrapolate the survival probability to the new maturity.

            dt = t_mat - self._times[-1]

            q_guess = self._qs[-1] * np.exp(-hazard_guess * dt)

            # Append the new curve node.

            self._times = np.append(self._times, t_mat)
            self._qs = np.append(self._qs, q_guess)

            q_star = optimize.newton(
                _f,
                x0=q_guess,
                fprime=None,
                args=argtuple,
                tol=1e-7,
                maxiter=50,
                fprime2=None,
            )

            if q_star < 0.0 or q_star > 1.0:
                raise FinError("Calibrated survival probability out of bounds")

            self.set_last_q(q_star)

    ###########################################################################

    def fwd(self, fwd_dt: Date):
        """Calculate the instantaneous forward rate at date fwd_dt
        using a numerical derivative."""

        t = input_time(fwd_dt, self)
        epsilon = 1e-4
        df1 = self.df(t) * self.survival_prob(t)
        df2 = self.df(t + epsilon) * self.survival_prob(t + epsilon)
        fwd = np.log(df1 / df2) / epsilon
        return fwd

    ###########################################################################

    def zero_rate(self, dt, freq_type=FrequencyTypes.CONTINUOUS):
        """Calculate the zero rate to date dt in the chosen compounding
        frequency where -1 is continuous is the default."""

        t = input_time(dt, self)
        f = annual_frequency(freq_type)
        df = self.df(t)
        q = self.survival_prob(t)
        dfq = df * q

        if f == 0:  # Simple interest
            zero_rate = (1.0 / dfq - 1.0) / t
        elif f == -1:  # Continuous
            zero_rate = -np.log(dfq) / t
        else:
            zero_rate = (dfq ** (-1.0 / (f * t)) - 1) * f

        return zero_rate

    ###########################################################################

    def log_curve(self):
        """Log CDS calibration inputs and the calibrated survival curve."""

        log_section("CDS Curve")

        log_value(
            "Anchor Date", self.anchor_dt,
            "Recovery Rate", f"{self.recovery_rate:.2%}",
            "Interpolation Type", self.interp_method,
            "Curve Day Count", self.dc_type,
            level=2,
        )

        # CDS contracts used to calibrate the curve
        log_table(
            "CDS Maturity",
            "CDS Spread (bp)",
            level=2,
        )

        for cds in self.cds_contracts:
            log_row(
                cds.maturity_dt,
                cds.running_cpn * 10000.0,
            )

        # Calibrated survival probability curve
        log_table(
            "Date",
            "Time",
            "Surv Prob",
            level=3,
        )

        # Anchor point
        log_row(
            self.anchor_dt,
            self._times[0],
            self._qs[0],
        )

        # Calibrated CDS maturity nodes
        for i, cds in enumerate(self.cds_contracts):
            log_row(
                cds.maturity_dt,
                self._times[i + 1],
                self._qs[i + 1],
            )

    ###########################################################################

    def __repr__(self):
        """Print out the details of the survival probability curve."""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        header = "TIME,SURVIVAL_PROBABILITY"

        value_table = [self._times, self._qs]
        precision = "10.7f"
        s += table_to_string(header, value_table, precision)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
