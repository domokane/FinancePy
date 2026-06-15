##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import least_squares

from ...utils.helpers import check_argument_types, _func_name

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.math import scale
from ...utils.helpers import label_to_string
from ...utils.helpers import times_from_dates

from ...market.curves.discount_curve import DiscountCurve

from .curve_fits import CurveFitMethod

########################################################################################

def get_fit_bounds(fit, n_params):

    bounds = fit.bounds

    lo = np.asarray(bounds[0], dtype=float)
    hi = np.asarray(bounds[1], dtype=float)

    if lo.ndim == 0:
        lo = np.full(n_params, lo)
    if hi.ndim == 0:
        hi = np.full(n_params, hi)

    if len(lo) != n_params or len(hi) != n_params:
        raise FinError(f"Bounds length {len(lo)},{len(hi)} does not match number of params {n_params}.")

    return lo, hi

########################################################################################

def f_fast(params, *args):

    flow_times, flow_amounts, dirty_prices, curve_fit, t_max = args

    if not np.all(np.isfinite(params)):
        return 1.0e25 * np.ones(len(flow_times))

    curve_fit.set_params(params)
    errors = np.zeros(len(flow_times))

    for i_bond in range(len(flow_times)):
        times_i = flow_times[i_bond]
        amounts_i = flow_amounts[i_bond]
        tau = times_i / t_max
        zero_rates = curve_fit.interp_rate(tau)

        if not np.all(np.isfinite(zero_rates)):
            errors[i_bond] = 1.0e25
            continue

        expo = np.clip(-zero_rates * times_i, -100.0, 100.0)
        dfs = np.exp(expo)

        pv = 100.0 * np.sum(amounts_i * dfs)

        if not np.isfinite(pv):
            errors[i_bond] = 1.0e25
            continue

        duration = np.sum(times_i * amounts_i * dfs) / np.sum(amounts_i * dfs)
        pv01 = pv * duration

        errors[i_bond] = (pv - dirty_prices[i_bond]) / max(pv01, 1e-8)

    return errors
########################################################################################


class BondParametricDiscountCurve(DiscountCurve):
    """Fit a parametric discount curve to dirty bond prices."""

    def __init__(
        self,
        value_dt: Date,
        bonds: list,
        clean_prices: list | np.ndarray,
        curve_fit: CurveFitMethod,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
        do_build: bool = True
        ):

        check_argument_types(getattr(self, _func_name(), None), locals())

        if len(bonds) != len(clean_prices):
            raise FinError("Num bonds does not equal number of prices.")

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type
        self.value_dt = value_dt
        self._curve_fit = curve_fit
        self._interp_type = None

        self.used_bonds = bonds
        self._validate_inputs()

        clean_prices = np.array(clean_prices, dtype=float)
        if np.any(clean_prices <= 0.0):
            raise FinError("Clean prices must be positive.")
        self.clean_prices = clean_prices

        self._precompute_bond_flows()

        self._t_mats = []
        for bond in bonds:
            t_mat = times_from_dates(self.value_dt,
                                     bond.maturity_dt,
                                     self.time_dc_type)
            self._t_mats.append(t_mat)

        self._t_max = max(np.max(self._t_mats), 1.0e-8)

        if do_build:
            self.build_curve()

    ###########################################################################

    @property
    def curve_fit(self):
        """Accessor function for curve_fit."""
        return self._curve_fit

    ####################################################################################

    def build_curve(self):

        if not self._bond_flow_times:
            raise FinError("No bond cash flows available for fitting.")

        dirty_prices = self.clean_prices + np.array(self._accrued)

        args = (
            self._bond_flow_times,
            self._bond_flow_amounts,
            dirty_prices,
            self._curve_fit,
            self._t_max,
        )

        x0 = self._curve_fit.get_params()

        bounds = get_fit_bounds(self._curve_fit, len(x0))

        lo, hi = bounds

        x0 = np.asarray(x0, dtype=float)

        lower_finite = np.isfinite(lo)
        upper_finite = np.isfinite(hi)

        x0[lower_finite] = np.maximum(x0[lower_finite], lo[lower_finite] + 1.0e-10)
        x0[upper_finite] = np.minimum(x0[upper_finite], hi[upper_finite] - 1.0e-10)

        result = least_squares(
            f_fast,
            x0=x0,
            args=args,
            bounds=bounds,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
            max_nfev=1000,
        )

        if result.status <= 0:
            raise FinError(f"Curve fitting failed: {result.message}")

        self._curve_fit.set_params(result.x)

        self._times = np.concatenate(([0.0], self._t_mats))
        taus = self._times / self._t_max
        self._zero_rates = self._curve_fit.interp_rate(taus)
        self._dfs = np.exp(-self._times * self._zero_rates)

    ####################################################################################

    def _validate_inputs(self):
        """Validate bond inputs are non-empty and in increasing maturity order."""

        num_bonds = len(self.used_bonds)
        bonds = self.used_bonds

        if num_bonds == 0:
            raise FinError("No calibration instruments.")

        if num_bonds > 1:
            # Bonds must be increasing in tenor/maturity
            prev_dt = bonds[0].maturity_dt
            for bond in bonds[1:]:
                next_dt = bond.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("Bonds must be in increasing maturity")
                prev_dt = next_dt

    ####################################################################################

    def _precompute_bond_flows(self):
        # We can speed up curve fitting by pre-computing and storing bond flows
        self._bond_flow_times = []
        self._bond_flow_amounts = []
        self._accrued = []

        for bond in self.used_bonds:
            bond.accrued_interest(self.value_dt, bond.par)
            self._accrued.append(bond.accrued_int)

            times = []
            amounts = []

            for cpn_dt, pmt_dt, flow in zip(
                bond.cpn_dts,
                bond.payment_dts,
                bond.flow_amounts,
            ):
                if cpn_dt > self.value_dt:
                    amt = flow
                    if pmt_dt == bond.payment_dts[-1]:
                        amt += bond.par / 100.0

                    t = times_from_dates(
                        self.value_dt,
                        pmt_dt,
                        self.time_dc_type,
                    )
                    times.append(t)
                    amounts.append(amt)

            self._bond_flow_times.append(np.array(times))
            self._bond_flow_amounts.append(np.array(amounts))

    ####################################################################################

    def df_t(self, t):
        """Discount factor from fitted parametric zero-rate curve."""

        times, scalar_input = self._to_time_array(t)
        taus = times / self._t_max
        zero_rates = self.curve_fit.interp_rate(taus)

        expo = np.clip(-zero_rates * times, -100.0, 100.0)
        dfs = np.exp(expo)
        dfs = np.maximum(dfs, 1.0e-300)

        if scalar_input:
            return float(dfs[0])
        return dfs

    ####################################################################################

    def bond_price_errors(self):

        fitted_clean_prices = []

        for i_bond in range(len(self.used_bonds)):

            times_i = self._bond_flow_times[i_bond]
            amounts_i = self._bond_flow_amounts[i_bond]

            tau = times_i / self._t_max
            zero_rates = self._curve_fit.interp_rate(tau)

            dfs = np.exp(-zero_rates * times_i)

            dirty_price_fit = 100.0 * np.sum(amounts_i * dfs)
            clean_price_fit = dirty_price_fit - self._accrued[i_bond]
            fitted_clean_prices.append(clean_price_fit)

        fitted_clean_prices = np.array(fitted_clean_prices)

        clean_errors = fitted_clean_prices - self.clean_prices

        return (
            fitted_clean_prices,
            clean_errors,
        )

    ####################################################################################

    def bond_yield_errors(self):

        n = len(self.used_bonds)

        maturities = np.array(self._t_mats, dtype=float)
        market_clean = self.clean_prices

        fitted_clean = np.zeros(n)
        market_ytm = np.zeros(n)
        fitted_ytm = np.zeros(n)

        for i, bond in enumerate(self.used_bonds):

            times_i = self._bond_flow_times[i]
            amounts_i = self._bond_flow_amounts[i]
            accrued_i = self._accrued[i]

            dfs_i = self.df_t(times_i)

            fitted_dirty = 100.0 * np.sum(amounts_i * dfs_i)
            fitted_clean[i] = fitted_dirty - accrued_i

            market_ytm[i] = bond.yield_to_maturity(
                self.value_dt,
                market_clean[i],
            )

            fitted_ytm[i] = bond.yield_to_maturity(
                self.value_dt,
                fitted_clean[i],
            )

        ytm_error = fitted_ytm - market_ytm

        return {
            "maturities": maturities,
            "market_ytm": market_ytm,
            "fitted_ytm": fitted_ytm,
            "ytm_error": ytm_error,
        }

    ####################################################################################

    def rms_yield_error(self):
        """RMS yield error in basis points."""
        out = self.bond_yield_errors()
        error_bp = 10000.0 * out["ytm_error"]
        return np.sqrt(np.mean(error_bp ** 2))

    ####################################################################################

    def rms_price_error(self):
        """RMS price error in cents (per $100 face)."""
        _, clean_errors = self.bond_price_errors()
        return np.sqrt(np.mean(clean_errors ** 2))

    ####################################################################################

    def plot_bond_yield_fit(self, title="Bond yield fit"):

        out = self.bond_yield_errors()

        t = out["maturities"]
        market_ytm = 100.0 * out["market_ytm"]
        fitted_ytm = 100.0 * out["fitted_ytm"]

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to maturity")
        plt.ylabel("Yield (%)")

        plt.plot(t, market_ytm, "o", label="Market YTM")
        plt.plot(t, fitted_ytm, "-", lw =2)

        plt.legend(loc="best")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    ####################################################################################

    def plot_bond_yield_errors(self, title="Bond yield fit errors"):

        out = self.bond_yield_errors()

        t = out["maturities"]
        error_bp = 10000.0 * out["ytm_error"]

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to maturity")
        plt.ylabel("Yield error (bp)")

        plt.axhline(0.0, linestyle="--", linewidth=1.0)
        plt.plot(t, error_bp, "o-", label="Fitted - market")

        rmse_bp = np.sqrt(np.mean(error_bp * error_bp))
        max_abs_bp = np.max(np.abs(error_bp))

        plt.legend(
            title=f"RMSE={rmse_bp:.3f} bp, MaxAbs={max_abs_bp:.3f} bp",
            loc="best",
        )

        plt.grid(True)
        plt.tight_layout()
        return plt

    ####################################################################################

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
        t = np.maximum(self._times, 1.0e-6)

        z = self.fwd_rate_inst_t(t)
        z = scale(z, 100.0)

        plt.plot(t, z, label=str(self._curve_fit))

        plt.legend(loc="lower right")
        plt.ylim((min(z) - 0.3, max(z) * 1.1))
        plt.grid(True)
        return plt

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUE DATE", self.value_dt)
#        s += label_to_string("BONDS", self.used_bonds)
        s += label_to_string("CLEAN PRICES", self.clean_prices)
        s += label_to_string("CURVE FIT", self._curve_fit)
        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)

    ##############################################################################
