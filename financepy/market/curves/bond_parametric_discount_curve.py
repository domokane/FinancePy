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

from .curve_fits import CurveFitTypes
from .curve_fits import CurveFitCubicPolynomial
from .curve_fits import CurveFitQuarticPolynomial
from .curve_fits import CurveFitQuinticPolynomial
from .curve_fits import CurveFitNelsonSiegel
from .curve_fits import CurveFitNelsonSiegelSvensson
from .curve_fits import CurveFitBSpline

from ...utils.format_graphs import set_plot_style
set_plot_style()

################################################################################


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


################################################################################


def f_fast(params, *args):

    flow_times, flow_amounts, dirty_prices, curve_fit = args

    if not np.all(np.isfinite(params)):
        return 1.0e25 * np.ones(len(flow_times))

    curve_fit.set_params(params)

    errors = np.zeros(len(flow_times))

    for i_bond in range(len(flow_times)):

        times_i = flow_times[i_bond]
        amounts_i = flow_amounts[i_bond]

        zero_rates = curve_fit.interp_rate(times_i)

        if not np.all(np.isfinite(zero_rates)):
            errors[i_bond] = 1.0e25
            continue

        expo = np.clip(-zero_rates * times_i, -100.0, 100.0)
        dfs = np.exp(expo)

        pv = 100.0 * np.sum(amounts_i * dfs)

        if not np.isfinite(pv):
            errors[i_bond] = 1.0e25
            continue

        denom = np.sum(amounts_i * dfs)

        if denom <= 0.0:
            errors[i_bond] = 1.0e25
            continue

        duration = np.sum(times_i * amounts_i * dfs) / denom

        pv01 = pv * duration

        errors[i_bond] = (pv - dirty_prices[i_bond]) / max(pv01, 1.0e-8)

    return errors


################################################################################


class BondParametricDiscountCurve(DiscountCurve):
    """Fit a parametric discount curve to bond prices."""

    def __init__(
        self,
        anchor_dt: Date,
        bonds: list,
        clean_prices: list | np.ndarray,
        curve_fit_type: CurveFitTypes,
        curve_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
        do_build: bool = True,
    ) -> None:

        check_argument_types(getattr(self, _func_name(), None), locals())

        if len(bonds) != len(clean_prices):
            raise FinError("Num bonds does not equal number of prices.")

        if not isinstance(curve_fit_type, CurveFitTypes):
            raise FinError("Invalid curve fit type.")

        if not isinstance(curve_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.anchor_dt = anchor_dt
        self.curve_dc_type = curve_dc_type
        self.curve_fit_type = curve_fit_type

        #######################################################################
        # Create fitter
        #######################################################################

        if curve_fit_type == CurveFitTypes.CUBIC_POLYNOMIAL:

            self._curve_fit = CurveFitCubicPolynomial()

        elif curve_fit_type == CurveFitTypes.QUARTIC_POLYNOMIAL:

            self._curve_fit = CurveFitQuarticPolynomial()

        elif curve_fit_type == CurveFitTypes.QUINTIC_POLYNOMIAL:

            self._curve_fit = CurveFitQuinticPolynomial()

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL:

            self._curve_fit = CurveFitNelsonSiegel()

        elif curve_fit_type == CurveFitTypes.NELSON_SIEGEL_SVENSSON:
            self._curve_fit = CurveFitNelsonSiegelSvensson()

        elif curve_fit_type == CurveFitTypes.BSPLINE:

            self._curve_fit = CurveFitBSpline()

        else:

            raise FinError("Unrecognised curve fit type.")

        self._interp_type = None

        self.used_bonds = bonds

        self._validate_inputs()

        clean_prices = np.asarray(clean_prices, dtype=float)

        if np.any(clean_prices <= 0.0):
            raise FinError("Clean prices must be positive.")

        self.clean_prices = clean_prices

        self._precompute_bond_flows()

        self._t_mats = []

        for bond in bonds:

            t_mat = times_from_dates(
                self.anchor_dt,
                bond.maturity_dt,
                self.curve_dc_type,
            )

            self._t_mats.append(t_mat)

        self._t_mats = np.asarray(self._t_mats, dtype=float)

        self.t_max = max(np.max(self._t_mats), 1.0e-8)

        #######################################################################
        # Set polynomial scale once
        #######################################################################

        if self.curve_fit_type in (
            CurveFitTypes.CUBIC_POLYNOMIAL,
            CurveFitTypes.QUARTIC_POLYNOMIAL,
            CurveFitTypes.QUINTIC_POLYNOMIAL,
        ):

            self._curve_fit.t_scale = self.t_max

        #######################################################################

        if do_build:
            self.build_curve()

    ###########################################################################

    @property
    def curve_fit(self):
        return self._curve_fit

    ###########################################################################

    def build_curve(self):

        if not self._bond_flow_times:
            raise FinError("No bond cash flows available for fitting.")

        dirty_prices = self.clean_prices + np.asarray(self._accrued)

        args = (
            self._bond_flow_times,
            self._bond_flow_amounts,
            dirty_prices,
            self._curve_fit,
        )

        x0 = self._curve_fit.get_params()

        bounds = get_fit_bounds(
            self._curve_fit,
            len(x0),
        )

        lo, hi = bounds

        x0 = np.asarray(x0, dtype=float)

        lower_finite = np.isfinite(lo)
        upper_finite = np.isfinite(hi)

        x0[lower_finite] = np.maximum(
            x0[lower_finite],
            lo[lower_finite] + 1.0e-10,
        )

        x0[upper_finite] = np.minimum(
            x0[upper_finite],
            hi[upper_finite] - 1.0e-10,
        )

        result = least_squares(
            f_fast,
            x0=x0,
            args=args,
            bounds=bounds,
            xtol=1.0e-10,
            ftol=1.0e-10,
            gtol=1.0e-10,
            max_nfev=1000,
        )

        if result.status <= 0:
            raise FinError(f"Curve fitting failed: {result.message}")

        self._curve_fit.set_params(result.x)

        self._times = np.concatenate(
            (
                [0.0],
                self._t_mats,
            )
        )

        self._zero_rates = self._curve_fit.interp_rate(self._times)

        self._dfs = np.exp(-self._times * self._zero_rates)

    ###########################################################################

    def _validate_inputs(self):

        num_bonds = len(self.used_bonds)

        if num_bonds == 0:
            raise FinError("No calibration instruments.")

        if num_bonds > 1:

            prev_dt = self.used_bonds[0].maturity_dt

            for bond in self.used_bonds[1:]:

                if bond.maturity_dt <= prev_dt:
                    raise FinError("Bonds must be in increasing maturity")

                prev_dt = bond.maturity_dt

    ###########################################################################

    def _precompute_bond_flows(self):

        self._bond_flow_times = []
        self._bond_flow_amounts = []
        self._accrued = []

        for bond in self.used_bonds:

            bond.accrued_interest(
                self.anchor_dt,
                bond.par,
            )

            self._accrued.append(bond.accrued_int)

            times = []
            amounts = []

            for cpn_dt, pmt_dt, flow in zip(
                bond.cpn_dts,
                bond.payment_dts,
                bond.flow_amounts,
            ):

                if cpn_dt > self.anchor_dt:

                    amt = flow

                    if pmt_dt == bond.payment_dts[-1]:
                        amt += bond.par / 100.0

                    t = times_from_dates(
                        self.anchor_dt,
                        pmt_dt,
                        self.curve_dc_type,
                    )

                    times.append(t)
                    amounts.append(amt)

            self._bond_flow_times.append(np.asarray(times))

            self._bond_flow_amounts.append(np.asarray(amounts))

    ###########################################################################

    def df_t(self, t):

        times, scalar_input = self._to_time_array(t)

        zero_rates = self._curve_fit.interp_rate(times)

        expo = np.clip(
            -zero_rates * times,
            -100.0,
            100.0,
        )

        dfs = np.exp(expo)

        dfs = np.maximum(
            dfs,
            1.0e-300,
        )

        if scalar_input:
            return float(dfs[0])

        return dfs

    ###########################################################################

    def bond_price_errors(self):

        fitted_clean_prices = []

        for i_bond in range(len(self.used_bonds)):

            times_i = self._bond_flow_times[i_bond]
            amounts_i = self._bond_flow_amounts[i_bond]

            zero_rates = self._curve_fit.interp_rate(times_i)

            dfs = np.exp(-zero_rates * times_i)

            dirty_price_fit = 100.0 * np.sum(amounts_i * dfs)

            clean_price_fit = dirty_price_fit - self._accrued[i_bond]

            fitted_clean_prices.append(clean_price_fit)

        fitted_clean_prices = np.asarray(fitted_clean_prices)

        clean_errors = fitted_clean_prices - self.clean_prices

        return fitted_clean_prices, clean_errors

    ###########################################################################

    def bond_yield_errors(self):

        n = len(self.used_bonds)

        maturities = np.asarray(
            self._t_mats,
            dtype=float,
        )

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
                self.anchor_dt,
                market_clean[i],
            )

            fitted_ytm[i] = bond.yield_to_maturity(
                self.anchor_dt,
                fitted_clean[i],
            )

        ytm_error = fitted_ytm - market_ytm

        return {
            "maturities": maturities,
            "market_ytm": market_ytm,
            "fitted_ytm": fitted_ytm,
            "ytm_error": ytm_error,
        }

    ###########################################################################

    def rms_yield_error(self):

        out = self.bond_yield_errors()

        error_bp = 10000.0 * out["ytm_error"]

        return np.sqrt(np.mean(error_bp**2))

    ###########################################################################

    def rms_price_error(self):

        _, clean_errors = self.bond_price_errors()

        return np.sqrt(np.mean(clean_errors**2))

    ###########################################################################

    def plot_bond_yield_fit(self, title="Bond yield fit"):

        out = self.bond_yield_errors()

        t = out["maturities"]

        market_ytm = 100.0 * out["market_ytm"]
        fitted_ytm = 100.0 * out["fitted_ytm"]

        plt.figure()
        plt.title(title)

        plt.xlabel("Time to maturity")
        plt.ylabel("Yield (%)")

        plt.plot(
            t,
            market_ytm,
            "o",
            label="Market YTM",
        )

        plt.plot(
            t,
            fitted_ytm,
            "-",
            lw=2,
            label="Fitted YTM",
        )

        plt.legend(loc="best")
        plt.grid(True)
        plt.show()

    ###########################################################################

    def plot_bond_yield_errors(self, title="Bond yield fit errors"):

        out = self.bond_yield_errors()

        t = out["maturities"]

        error_bp = 10000.0 * out["ytm_error"]

        plt.figure()

        plt.title(title)

        plt.xlabel("Time to maturity")
        plt.ylabel("Yield error (bp)")

        plt.axhline(
            0.0,
            linestyle="--",
            linewidth=1.0,
        )

        plt.plot(
            t,
            error_bp,
            "o-",
            label="Fitted - market",
        )

        rmse_bp = np.sqrt(np.mean(error_bp * error_bp))

        max_abs_bp = np.max(np.abs(error_bp))

        plt.legend(
            title=(f"RMSE={rmse_bp:.3f} bp, " f"MaxAbs={max_abs_bp:.3f} bp"),
            loc="best",
        )

        plt.grid(True)

        return plt

    ###########################################################################

    def plot_zero_rate(self, title, ylabel="Zero Rate (%)"):

        plt.figure()

        plt.title(title)

        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        t = self._times

        z = self.zero_rate(t)

        z = scale(z, 100.0)

        plt.plot(
            t,
            z,
            label=str(self._curve_fit),
        )

        plt.legend(loc="lower right")
        plt.grid(True)

        return plt

    ###########################################################################

    def plot_fwd_rate(self, title, ylabel="Forward Rate (%)"):

        plt.figure()

        plt.title(title)

        plt.xlabel("Time to Maturity (years)")
        plt.ylabel(ylabel)

        t = np.maximum(
            self._times,
            1.0e-6,
        )

        z = self.fwd_rate_inst_t(t)

        z = scale(z, 100.0)

        plt.plot(
            t,
            z,
            label=str(self._curve_fit),
        )

        plt.legend(loc="lower right")

        plt.ylim(
            (
                min(z) - 0.3,
                max(z) * 1.1,
            )
        )

        plt.grid(True)

        return plt

    ###########################################################################

    def __repr__(self):

        s = label_to_string(
            "OBJECT TYPE",
            type(self).__name__,
        )

        s += label_to_string(
            "ANCHOR DATE",
            self.anchor_dt,
        )

        s += label_to_string(
            "CLEAN PRICES",
            self.clean_prices,
        )

        s += label_to_string(
            "CURVE FIT TYPE",
            self.curve_fit_type,
        )

        s += label_to_string(
            "CURVE FIT",
            self._curve_fit,
        )

        return s

    ###########################################################################

    def _print(self):
        print(self)


##############################################################################
