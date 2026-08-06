##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Global fit due to V. Piterbarg 2024
##############################################################################

import matplotlib.pyplot as plt

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.helpers import check_argument_types, _func_name
from ...utils.helpers import times_from_dates
from ...utils.day_count import DayCountTypes
from ...utils.helpers import table_to_string
from .interpolator import InterpTypes, Interpolator
from ...utils.error import FinError
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string

FAST = True

########################################################################################


def _f_slow(df, *args):

    curve = args[0]
    value_dt = args[1]
    bond = args[2]
    mkt_clean_price = args[3]

    curve.set_last_df(df)
    curve.fit(curve._times, curve._dfs)

    bond_discount_price = bond.clean_price_from_discount_curve(value_dt, curve)
    obj_fn = bond_discount_price - mkt_clean_price

    return obj_fn

########################################################################################

def _f_fast(df, *args):

    curve, bond_index, mkt_clean_price = args

    curve.set_last_df(df)
    curve.fit(curve._times, curve._dfs)

    times = curve._bond_flow_times[bond_index]
    amounts = curve._bond_flow_amounts[bond_index]
    accrued = curve._accrued[bond_index]

    dfs = curve.df_t(times)
    dirty = 100.0 * np.sum(amounts * dfs)
    clean = dirty - accrued

    return clean - mkt_clean_price


########################################################################################


class BondBootstrapDiscountCurve(DiscountCurve):
    """Class to do bootstrap exact fitting of the bond discount curve."""

    def __init__(
        self,
        value_dt: Date,
        bonds: list,
        clean_prices: list | np.ndarray,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        time_dc_type: DayCountTypes = DayCountTypes.ACT_365F,
        check_refit_flag: bool = False,  # Set to True to test it works
        do_build: bool = True,
    ):

        check_argument_types(getattr(self, _func_name(), None), locals())

        if len(bonds) != len(clean_prices):
            raise FinError("Num bonds does not equal number of prices.")

        if not isinstance(time_dc_type, DayCountTypes):
            raise FinError("Invalid time day count type.")

        self.time_dc_type = time_dc_type

        self.check_refit_flag = check_refit_flag
        self._interp_type = interp_type
        self._interpolator = Interpolator(self._interp_type)
        self.is_built = False

        self.value_dt = value_dt

        self._validate_inputs(bonds)

        clean_prices = np.array(clean_prices, dtype=float)

        if np.any(clean_prices <= 0.0):
            raise FinError("Clean prices must be positive.")

        self.clean_prices = clean_prices

        if FAST:
            self._precompute_bond_flows()

        self._t_mats = []
        for bond in bonds:
            t_mat = times_from_dates(self.value_dt, bond.maturity_dt, self.time_dc_type)
            self._t_mats.append(t_mat)

        if do_build:
            self.build_curve()

    ###########################################################################

    def build_curve(self, **kwargs):
        """
        Build curve based on interpolation.

        Not all interpolators are suitable for the bootstrap/1d solver, only those
        that are local,
        where the value of df[i] does not affect discount factors for t<=t[i-1]
        """

        if Interpolator.suitable_for_bootstrap(self._interp_type):
            self._build_curve_using_1d_solver()
        else:
            self._build_curve_using_least_squares()

        self.is_built = True

    ####################################################################################

    def _validate_inputs(self, bonds):
        """Validate the inputs for each of the Ibor products."""

        num_bonds = len(bonds)

        if num_bonds == 0:
            raise FinError("No calibration instruments.")

        if num_bonds > 1:

            # Swaps must be increasing in tenor/maturity
            prev_dt = bonds[0].maturity_dt
            for bond in bonds[1:]:
                next_dt = bond.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("Bonds must be in increasing maturity")
                prev_dt = next_dt

        # Now determine which instruments are used
        self.used_bonds = bonds

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
                        amt += 1.0

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

    def _build_curve_using_1d_solver(self):

        self._interpolator = Interpolator(self._interp_type)
        self._times = np.array([])
        self._dfs = np.array([])
        self.is_built = True

        t_mat = 0.0
        df_mat = 1.0

        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, df_mat)
        self.fit(self._times, self._dfs)

        for i, bond in enumerate(self.used_bonds):

            clean_price = self.clean_prices[i]

            # The last payment for the bond might be adjusted so I use that rather
            # than using the maturity date itself. It's small but not trivial.
            maturity_dt = bond.payment_dts[-1]

            t_mat = times_from_dates(
                self.value_dt,
                maturity_dt,
                self.time_dc_type,
            )

            if FAST:

                argtuple = (self, i, clean_price)

                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)

                df_mat = optimize.newton(
                    _f_fast,
                    x0=df_mat,
                    fprime=None,
                    args=argtuple,
                    tol=1e-5,
                    maxiter=100,
                )

            else:
                argtuple = (self, self.value_dt, bond, clean_price)

                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)

                df_mat = optimize.newton(
                    _f_slow,
                    x0=df_mat,
                    fprime=None,
                    args=argtuple,
                    tol=1e-5,
                    maxiter=100,
                )

            # I set this explicitly to be sure
            self.set_last_df(df_mat)
            self.fit(self._times, self._dfs)

            if df_mat <= 0.0:
                raise FinError("Bootstrapped discount factor must be positive")

        if self.check_refit_flag is True:
            self.check_refit(1e-6)

    ###########################################################################

    def _build_curve_using_least_squares(self, **kwargs):
        """
        Construct the discount curve using a least-squares minimisation approach.
        This enables a more complex interpolation scheme.
        """

        def _obj_f_curve_build_ls_slow(dfs):
            """
            Objective function for fitting all knot dfs at once to the benchmark
            securities  -- suitable for non-local interpolators
            """

            bond_curve = self
            value_dt = bond_curve.value_dt
            bond_curve._dfs[1:] = dfs

            bond_curve.fit(bond_curve._times, bond_curve._dfs)

            n_bonds = len(self.used_bonds)
            out = np.zeros(n_bonds)
            idx = 0

            for i in range(0, n_bonds):
                p_fit = self.used_bonds[i].clean_price_from_discount_curve(
                    value_dt, bond_curve
                )
                p_mkt = self.clean_prices[i]
                out[idx] = p_fit - p_mkt
                idx = idx + 1

            return out

        def _obj_f_curve_build_ls_fast(dfs):
            """
            Objective function for fitting all knot dfs at once to the benchmark
            securities -- suitable for non-local interpolators.

            Uses precomputed bond flows to avoid calling
            bond.clean_price_from_discount_curve().
            """

            bond_curve = self
            bond_curve._dfs[1:] = dfs
            bond_curve.fit(bond_curve._times, bond_curve._dfs)

            n_bonds = len(self.used_bonds)
            out = np.zeros(n_bonds)

            for i in range(n_bonds):

                times = bond_curve._bond_flow_times[i]
                amounts = bond_curve._bond_flow_amounts[i]
                accrued = bond_curve._accrued[i]

                dfs_i = bond_curve.df_t(times)

                dirty = 100.0 * np.sum(amounts * dfs_i)
                clean = dirty - accrued

                out[i] = clean - self.clean_prices[i]

            return out

        self.is_built = True

        orig_check_refit_flag = self.check_refit_flag
        self.check_refit_flag = False

        # Basic bootstrap using non-local interpolation

        self._build_curve_using_1d_solver()
        self.check_refit_flag = orig_check_refit_flag

        if FAST:
            res = optimize.least_squares(
                _obj_f_curve_build_ls_fast,
                self._dfs[1:],
                bounds=(0, 1.0),
                ftol=1e-4,
                xtol=1e-5,
            )
        else:
            res = optimize.least_squares(
                _obj_f_curve_build_ls_slow,
                self._dfs[1:],
                bounds=(0, 1.0),
                ftol=1e-4,
                xtol=1e-5,
            )

        if not res.success:
            raise FinError(res.message)

        self._dfs[1:] = np.array(res.x)
        self.fit(self._times, self._dfs)

        if self.check_refit_flag is True:
            self.check_refit(1e-6)

    ####################################################################################

    def check_refit(self, bond_tol):
        """Ensure that the Bond curve refits the calibration instruments."""

        for i, bond in enumerate(self.used_bonds):
            # We value it as of the start date of the swap
            p_mod = bond.clean_price_from_discount_curve(self.value_dt, self)
            p_mkt = self.clean_prices[i]
            diff = p_mod - p_mkt

            if abs(diff) > bond_tol:
                print(self._interp_type)
                print(
                    f"Bond with maturity {bond.maturity_dt} "
                    f"not repriced. Model Price={p_mod:.8f} "
                    f"Market Price={p_mkt:.8f} "
                    f"Diff={diff:.8f}"
                )
                raise FinError("Bond not repriced.")

    ####################################################################################

    def plot_zero_rates(self, title: str):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Zero Rate (%)")

        tmax = np.max(self._t_mats)
        t = np.linspace(1.0e-6, float(tmax), 100)

        zero_rate = self.zero_rate_t(t)
        plt.plot(t, zero_rate * 100, label="Zero Rate Bootstrap", marker="o")
        plt.legend(loc="lower right")
        plt.grid(True)

    ###########################################################################

    def plot_fwd_rates(self, title: str):
        """Display yield curve."""

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel("Time to Maturity (years)")
        plt.ylabel("Forward Rate (%)")

        tmax = np.max(self._t_mats)
        t = np.linspace(1.0e-6, float(tmax), 100)

        fwd_rate = self.fwd_rate_inst_t(t)
        plt.plot(t, fwd_rate * 100, label="Fwd Rate Bootstrap", marker="o")
        plt.legend(loc="lower right")
        plt.grid(True)

    ###########################################################################

    def __repr__(self):
        # TODO
        header = "TIMES,DISCOUNT FACTORS"
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        value_table = [self._times, self._dfs]
        precision = "10.7f"
        s += table_to_string(header, value_table, precision)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
