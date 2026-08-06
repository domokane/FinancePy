##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_vars import G_DAYS_IN_YEAR
from ...market.curves.interpolator import _uinterpolate, InterpTypes
from ...utils.helpers import input_time, table_to_string
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.helpers import check_argument_types, _func_name
from ...utils.helpers import label_to_string

DIRTY = 0
CLEAN = 1

# from numba import njit, float64

########################################################################################


def f(q, *args):
    """Function that returns zero when the survival probability that gives a
    zero value of the CDS has been determined."""

    curve = args[0]
    value_dt = args[1]
    cds = args[2]
    recovery_rate = args[3]

    curve.set_last_q(q)

    # This is important - we calibrate a curve that makes the clean PV of the
    # CDS equal to zero and so we select the second element of the value tuple
    obj_fn = cds.value(value_dt, curve, recovery_rate)[CLEAN]
    return obj_fn


########################################################################################


class CDSCurve:
    """Generate a survival probability curve implied by the value of CDS
    contracts given a Ibor curve and an assumed recovery rate. The recovery
    rate corresponds to the seniority of the debt for these CDS. A scheme for
    the interpolation of the survival probabilities is also required."""

    def __init__(
        self,
        value_dt: Date,
        cds_contracts: list,
        libor_curve,
        recovery_rate,
        interp_method: InterpTypes = InterpTypes.FLAT_FWD_RATES,
    ):
        """Construct a credit curve from a sequence of maturity-ordered CDS
        contracts and a Ibor curve using the same recovery rate and the
        same interpolation method."""

        check_argument_types(getattr(self, _func_name(), None), locals())

        if value_dt != libor_curve.value_dt:
            raise FinError(
                "Curve does not have same valuation date as Issuer curve."
            )

        self.value_dt = value_dt
        self.cds_contracts = cds_contracts
        self.recovery_rate = recovery_rate
        self.libor_curve = libor_curve
        self.interp_method = interp_method
        self.built_ok = False

        self._times = []
        self._qs = []

        if len(self.cds_contracts) > 0:
            self.build_curve()
        else:
            pass  # In some cases we allow None to be passed

        return

    ###########################################################################

    # @property
    # def times(self):
    #     return self._times.copy()

    # ###########################################################################

    # @property
    # def qs(self):
    #     return self._qs.copy()

    # @property
    # def dfs(self):
    #     return self.libor_curve.dfs

    def set_times(self, times: np.array):
        """Set the times vector"""
        self._times = times

    def set_qs(self, qvector: np.array):
        """Set the survival probability curve"""
        self._qs = qvector

    def set_q(self, index, q):
        """Set the survival probability at a specific index."""

        n_points = len(self._qs)

        if index < 0 or index >= n_points:
            raise IndexError("Index out of bounds")

        self._qs[index] = q

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
            t = (dt - self.value_dt) / G_DAYS_IN_YEAR
        elif isinstance(dt, list):
            t = np.array(dt)
        else:
            t = dt

        if np.any(t < 0.0):
            raise FinError("Survival Date before curve anchor date")

        if isinstance(t, np.ndarray):
            n = len(t)
            qs = np.zeros(n)
            for i in range(0, n):
                qs[i] = _uinterpolate(
                    t[i], self._times, self._qs, self.interp_method.value
                )
            return qs
        elif np.isscalar(t):
            q = _uinterpolate(
                t, self._times, self._qs, self.interp_method.value
            )
            return q

        raise FinError("Unknown time type")

    ###########################################################################

    def df(self, dt):
        """Extract the discount factor from the underlying Ibor curve. This
        function supports vectorisation."""

        if isinstance(dt, Date):
            t = (dt - self.value_dt) / G_DAYS_IN_YEAR
        elif isinstance(dt, list):
            t = np.array(dt)
        else:
            t = dt

        df = self.libor_curve.df_t(t)

        return df

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

            argtuple = (
                self,
                self.value_dt,
                self.cds_contracts[i],
                self.recovery_rate,
            )

            t_mat = (maturity_dt - self.value_dt) / G_DAYS_IN_YEAR
            q = self._qs[i]

            self._times = np.append(self._times, t_mat)
            self._qs = np.append(self._qs, q)

            q_star = optimize.newton(
                f,
                x0=q,
                fprime=None,
                args=argtuple,
                tol=1e-7,
                maxiter=50,
                fprime2=None,
            )

            if q_star < 0.0 or q_star > 1.0:
                raise FinError("Calibrated survival probability out of bounds")

            # TODO - DETERMINE WHY THIS FAILS
            #            if i > 0 and q_star > self._qs[i]:
            #                print(self._qs)
            #                raise FinError("Survival probabilities must be non-increasing")

            self.set_last_q(q_star)

    ###########################################################################

    def fwd(self, fwd_dt):
        """Calculate the instantaneous forward rate at date fwd_dt
        using a numerical derivative."""

        t = input_time(fwd_dt, self)
        epsilon = 1e-4
        df1 = self.df(t) * self.survival_prob(t)
        df2 = self.df(t + epsilon) * self.survival_prob(t + epsilon)
        fwd = np.log(df1 / df2) / epsilon
        return fwd

    ###########################################################################

    # def fwd_rate(self, date1, date2, dc_type):
    #     """Calculate the risky forward rate according between dates date1
    #     and date2 according to the specified day count convention."""

    #     print("WHY AM I USING THIS ???? fwd_rate cds_curve")

    #     if date1 < self.value_dt:
    #         raise FinError("Date1 before curve value date.")

    #     if date2 < date1:
    #         raise FinError("Date2 must not be before Date1")

    #     day_count = DayCount(dc_type)
    #     year_frac = day_count.year_frac(date1, date2)[0]
    #     df1 = self.df(date1)
    #     df2 = self.df(date2)
    #     fwd = (df1 / df2 - 1.0) / year_frac
    #     return fwd

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

    def __repr__(self):
        """Print out the details of the survival probability curve."""
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
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
