##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union
from numba import njit

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize


from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.math import scale
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates


from ...market.curves.discount_curve import DiscountCurve
from .curve_fits import CurveFitPolynomial

###############################################################################


def f_polynomial(params, *args):
    """
    Sum of squared dirty-price errors for a Polynomial zero curve.
    """

    times = args[0]
    flows = args[1]
    dirty_prices = args[2]
    curve_fit = args[3]
    t_max = args[4]

    # rescale to $1 of par
    dirty_prices = scale(dirty_prices, 0.01)

    # Fit zero rates with polynomial
    if isinstance(curve_fit, CurveFitPolynomial):
        fit = CurveFitPolynomial()
        fit.coeffs = np.array(params, dtype=float)

    pv_err = 0.0

    for i_bond in range(0, len(times)):

        flow_times = times[i_bond]
        flow_amounts = flows[i_bond]

        pv = 0.0
        n_flows = len(flow_times)

        for i_flow in range(0, n_flows):

            t = flow_times[i_flow]
            tau = t / t_max
            cf = flow_amounts[i_flow]
            z = fit.interp_rate(tau)

            df = np.exp(-z * t)
            pv += cf * df

        pv += df
        err = (pv - dirty_prices[i_bond]) / dirty_prices[i_bond]
        pv_err += err**2

    return pv_err


###############################################################################


class BondFittedZeroCurve(DiscountCurve):
    """Fit a parametric discount curve to dirty bond prices."""

    def __init__(
        self,
        settle_dt: Date,
        bonds: list,
        dirty_prices: Union[np.ndarray, list],
        curve_fit,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
    ):

        if len(bonds) != len(dirty_prices):
            raise FinError("Number of bonds and dirty prices must match")

        self.value_dt = settle_dt
        self.settle_dt = settle_dt
        self._curve_fit = curve_fit

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type
        self._bonds = bonds
        self._dirty_prices = np.array(dirty_prices, dtype=float)

        x0 = np.array(self._curve_fit.coeffs, dtype=float)

        #        t_max = (bonds[-1].maturity_dt - settle_dt) / G_DAYS_IN_YEAR

        t_max = times_from_dates(
            bonds[-1].maturity_dt, settle_dt, self.time_dc_type
        )

        # Precalculate bond times and flows
        time_grid = []
        flow_grid = []

        for i_bond in range(0, len(bonds)):
            bond = bonds[i_bond]
            time_grid.append(bond.times(settle_dt))
            flow_grid.append(bond.flows(settle_dt))

        self.t_max = t_max

        args = (time_grid, flow_grid, self._dirty_prices, curve_fit, t_max)

        result = minimize(
            f_polynomial,
            x0=x0,
            args=args,
            method="L-BFGS-B",
        )

        self._curve_fit.coeffs = np.array(result.x, dtype=float)

        self._times = np.array(
            [
                times_from_dates(
                    bond.maturity_dt,
                    self.settle_dt,
                    self.time_dc_type,
                )
                for bond in self._bonds
            ],
            dtype=float,
        )

        self._zero_rates = np.array(
            [self._curve_fit.interp_rate(t / t_max) for t in self._times],
            dtype=float,
        )

        self._dfs = np.exp(-self._times * self._zero_rates)

    ###########################################################################

    @property
    def curve_fit(self):
        """Accessor function for curve_fit."""
        return self._curve_fit

    ###########################################################################

    def zero_rate(self, t):
        """Return the fitted zero rate for scalar or array time in years."""

        if isinstance(t, Date):
            dt = t
            t = times_from_dates(dt, self.settle_dt, self.time_dc_type)
            t_arr = np.array([t], dtype=float)
            scalar_input = True
        else:
            scalar_input = np.isscalar(t)

            if scalar_input:
                t_arr = np.array([float(t)], dtype=float)
            else:
                t_arr = np.array(t, dtype=float)

        t_max = self.t_max

        z = np.array(
            [self._curve_fit.interp_rate(x / t_max) for x in t_arr],
            dtype=float,
        )

        if scalar_input:
            return float(z[0])

        return z

    ###########################################################################

    def par_rate(self, t, freq):
        """Return the fitted par rate for scalar or array time in years."""
        scalar_input = np.isscalar(t)

        if scalar_input:
            t_arr = np.array([float(t)], dtype=float)
        else:
            t_arr = np.array(t, dtype=float)

        par_rates = []

        for t in t_arr:

            n_flows = int(t * freq) + 1

            z = self.curve_fit.interp_rate(t / self.t_max)
            df = np.exp(-z * t)
            num = 1.0 - df

            # We assume first payment is proportional to stub length
            t_stub = t - (n_flows - 1) / freq
            z_stub = self.curve_fit.interp_rate(t_stub / self.t_max)
            df_stub = np.exp(-z_stub * t_stub)
            pv01 = df_stub * t_stub

            t_flow = t_stub

            for i in range(1, n_flows):
                t_flow += 1.0 / freq
                z_flow = self.curve_fit.interp_rate(t_flow / self.t_max)
                df = np.exp(-z_flow * t_flow)
                pv01 += df / freq

            par_rate = num / pv01
            par_rates.append(par_rate)

        if scalar_input:
            return float(par_rates[0])

        return np.array(par_rates, dtype=float)

    ###########################################################################

    def fwd_rate(self, t):
        """Return the fitted fwd rate for scalar or array time in years."""
        scalar_input = np.isscalar(t)

        if scalar_input:
            t_arr = np.array([float(t)], dtype=float)
        else:
            t_arr = np.array(t, dtype=float)

        fwds = []

        if 1 == 1:
            epsilon = 1e-6
            for i in range(len(t_arr)):
                t1 = t_arr[i]
                t2 = t1 + epsilon
                df1 = self.df_t(t1)
                df2 = self.df_t(t2)
                f = (df1 / df2 - 1.0) / (t2 - t1)
                fwds.append(f)
        else:
            for i in range(len(t_arr) - 1):
                t1 = t_arr[i]
                t2 = t_arr[i + 1]
                df1 = self.df_t(t1)
                df2 = self.df_t(t2)
                f = (df1 / df2 - 1.0) / (t2 - t1)
                fwds.append(f)

            fwds.append(f)

        if scalar_input:
            return float(fwds[0])

        return np.array(fwds, dtype=float)

    ###########################################################################

    def df_t(self, t):
        """Return discount factors implied by the fitted zero curve."""
        scalar_input = np.isscalar(t)

        if scalar_input:
            t_arr = np.array([float(t)], dtype=float)
        else:
            t_arr = np.array(t, dtype=float)

        z = self.zero_rate(t_arr)

        dfs = np.exp(-z * t_arr)

        if scalar_input:
            return float(dfs[0])

        return dfs

    ###########################################################################

    def plot_zero_rate(self, title, ylabel="Zero Rate (%)"):
        """Display fitted zero-rate curve."""
        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        t = self._times
        z = self.zero_rate(t)
        z = scale(z, 100.0)

        plt.plot(t, z, label=str(self._curve_fit))

        plt.legend(loc="lower right")
        plt.grid(True)
        return plt

    ###########################################################################

    def plot_fwd_rate(self, title, ylabel="Forward Rate (%)"):
        """Display fitted zero-rate curve."""
        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        t = self._times
        z = self.fwd_rate(t)
        z = scale(z, 100.0)

        plt.plot(t, z, label=str(self._curve_fit))

        plt.legend(loc="lower right")
        plt.ylim((min(z) - 0.3, max(z) * 1.1))
        plt.grid(True)
        return plt

    ###########################################################################

    def plot_par_rate(self, title, ylabel="Par Rate (%)"):
        """Display fitted par-rate curve."""
        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        freq = 2
        t = self._times
        z = self.par_rate(t, freq)
        z = scale(z, 100.0)

        plt.plot(t, z, label=str(self._curve_fit))

        plt.legend(loc="lower right")
        plt.ylim((min(z) - 0.3, max(z) * 1.1))
        plt.grid(True)
        return plt

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SETTLEMENT DATE", self.settle_dt)
        s += label_to_string("BONDS", self._bonds)
        s += label_to_string("DIRTY PRICES", self._dirty_prices)
        s += label_to_string("CURVE FIT", self._curve_fit)
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)

    ##############################################################################
