##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types, _func_name
from ...utils.global_vars import g_days_in_year
from ...market.curves.interpolator import InterpTypes
from ...market.curves.discount_curve import DiscountCurve

SWAP_TOL = 1e-10

##############################################################################
# TODO: CHANGE times to df_times
##############################################################################


def _f(df, *args):
    """Root search objective function for swaps"""
    curve = args[0]
    value_dt = args[1]
    swap = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df
    v_swap = swap.value(value_dt, curve, curve, None, 1.0)
    v_swap /= swap.notional
    return v_swap


###############################################################################


def _g(df, *args):
    """Root search objective function for swaps"""
    curve = args[0]
    value_dt = args[1]
    fra = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df
    v_fra = fra.value(value_dt, curve)
    v_fra /= fra.notional
    return v_fra


###############################################################################


class InflationSwapCurve(DiscountCurve):
    """Constructs a discount curve as implied by the prices of Ibor
    deposits, FRAs and IRS. The curve date is the date on which we are
    performing the valuation based on the information available on the
    curve date. Typically it is the date on which an amount of 1 unit paid
    has a present value of 1. This class inherits from FinDiscountCurve
    and so it has all of the methods that that class has.

    There are two main curve-building approaches:

    1) The first uses a bootstrap that interpolates swap rates linearly for
    cpn dates that fall between the swap maturity dates. With this, we can
    solve for the discount factors iteratively without need of a solver. This
    will give us a set of discount factors on the grid dates that refit the
    market exactly. However, when extracting discount factors, we will then
    assume flat forward rates between these cpn dates. There is no
    contradiction as it is as though we had been quoted a swap curve with all
    of the market swap rates, and with an additional set as though the market
    quoted swap rates at a higher frequency than the market.

    2) The second uses a bootstrap that uses only the swap rates provided but
    which also assumes that forwards are flat between these swap maturity
    dates. This approach is non-linear and so requires a solver. Consequently
    it is slower. Its advantage is that we can switch interpolation schemes
    to provide a smoother or other functional curve shape which may have a more
    economically justifiable shape. However the root search makes it slower."""

    ###########################################################################

    def __init__(
        self,
        value_dt: Date,
        ibor_deposits: list,
        ibor_fras: list,
        ibor_swaps: list,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        check_refit: bool = False,
    ):  # Set to True to test it works
        """Create an instance of a Ibor curve given a valuation date and
        a set of ibor deposits, ibor FRAs and ibor_swaps. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on cpn dates and to then assume
        flat forwards between these cpn dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        """
        raise FinError("DO NOT USE")

        check_argument_types(getattr(self, _func_name(), None), locals())

        self.value_dt = value_dt
        self.validate_inputs(ibor_deposits, ibor_fras, ibor_swaps)
        self._interp_type = interp_type
        self.check_refit = check_refit
        self.build_curve()

    ###########################################################################

    def _build_curve(self):
        """Build curve based on interpolation."""
        self.build_curve_using_solver()

    ###########################################################################

    def _validate_inputs(self, ibor_deposits, ibor_fras, ibor_swaps):
        """Validate the inputs for each of the Ibor products."""

        num_depos = len(ibor_deposits)
        num_fras = len(ibor_fras)
        num_swaps = len(ibor_swaps)

        if num_depos + num_fras + num_swaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if num_depos > 0:
            for depo in ibor_deposits:
                start_dt = depo.start_dt
                if start_dt < self.value_dt:
                    raise FinError("First deposit starts before value date.")

            for depo in ibor_deposits:
                start_dt = depo.start_dt
                end_dt = depo.maturity_dt
                if start_dt >= end_dt:
                    raise FinError("First deposit ends on or before it begins")

        # Ensure order of depos
        if num_depos > 1:
            prev_dt = ibor_deposits[0].maturity_dt

            for depo in ibor_deposits[1:]:
                next_dt = depo.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("Deposits must be in increasing maturity")
                prev_dt = next_dt

        # Ensure that valuation date is on or after first deposit start date
        if num_depos > 1:
            if ibor_deposits[0].start_dt > self.value_dt:
                raise FinError(
                    "Valuation date must not be before first deposit settles."
                )

        if num_fras > 0:
            for fra in ibor_fras:
                start_dt = fra.start_dt
                if start_dt <= self.value_dt:
                    raise FinError("FRAs starts before valuation date")

        if num_fras > 1:
            prev_dt = ibor_fras[0].maturity_dt
            for fra in ibor_fras[1:]:
                next_dt = fra.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("FRAs must be in increasing maturity")
                prev_dt = next_dt

        if num_swaps > 0:
            for swap in ibor_swaps:
                start_dt = swap.start_dt
                if start_dt < self.value_dt:
                    raise FinError("Swaps starts before valuation date.")

        if num_swaps > 1:

            # Swaps must all start on the same date for the bootstrap
            start_dt = ibor_swaps[0].start_dt
            for swap in ibor_swaps[1:]:
                next_start_dt = swap.start_dt
                if next_start_dt != start_dt:
                    raise FinError("Swaps must all have same start date.")

            # Swaps must be increasing in tenor/maturity
            prev_dt = ibor_swaps[0].maturity_dt
            for swap in ibor_swaps[1:]:
                next_dt = swap.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("Swaps must be in increasing maturity")
                prev_dt = next_dt

            # Swaps must have same cash flows for bootstrap to work
            longest_swap = ibor_swaps[-1]
            longest_swap_cpn_dts = longest_swap.adjusted_fixed_dts
            for swap in ibor_swaps[0:-1]:
                swap_cpn_dts = swap.adjusted_fixed_dts
                num_flows = len(swap_cpn_dts)
                for i_flow in range(0, num_flows):
                    if swap_cpn_dts[i_flow] != longest_swap_cpn_dts[i_flow]:
                        raise FinError(
                            "Swap cpns are not on the same date grid."
                        )

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        last_deposit_maturity_dt = Date(1, 1, 1900)
        first_fra_maturity_dt = Date(1, 1, 1900)
        last_fra_maturity_dt = Date(1, 1, 1900)

        if num_depos > 0:
            last_deposit_maturity_dt = ibor_deposits[-1].maturity_dt

        if num_fras > 0:
            first_fra_maturity_dt = ibor_fras[0].maturity_dt
            last_fra_maturity_dt = ibor_fras[-1].maturity_dt

        if num_swaps > 0:
            first_swap_maturity_dt = ibor_swaps[0].maturity_dt

        if num_depos > 0 and num_fras > 0:
            if first_fra_maturity_dt <= last_deposit_maturity_dt:
                print("FRA Maturity Date:", first_fra_maturity_dt)
                print("Last Deposit Date:", last_deposit_maturity_dt)
                raise FinError("First FRA must end after last Deposit")

        if num_fras > 0 and num_swaps > 0:
            if first_swap_maturity_dt <= last_fra_maturity_dt:
                raise FinError("First Swap must mature after last FRA")

        # Now determine which instruments are used
        self.used_deposits = ibor_deposits
        self.used_fras = ibor_fras
        self.used_swaps = ibor_swaps
        self.dc_type = None

    ###########################################################################

    def _build_curve_using_solver(self):
        """Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver."""

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        t_mat = 0.0
        df_mat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, df_mat)

        for depo in self.used_deposits:
            df_settle = self.df(depo.start_dt)
            df_mat = depo.maturity_df() * df_settle
            t_mat = (depo.maturity_dt - self.value_dt) / g_days_in_year
            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self._dfs, df_mat)

        oldt_mat = t_mat

        for fra in self.used_fras:

            t_set = (fra.start_dt - self.value_dt) / g_days_in_year
            t_mat = (fra.maturity_dt - self.value_dt) / g_days_in_year

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if t_set < oldt_mat and t_mat > oldt_mat:
                df_mat = fra.maturity_df(self)
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)
            else:
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)

                argtuple = (self, self.value_dt, fra)
                df_mat = optimize.newton(
                    _g,
                    x0=df_mat,
                    fprime=None,
                    args=argtuple,
                    tol=SWAP_TOL,
                    maxiter=50,
                    fprime2=None,
                )

        for swap in self.used_swaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturity_dt = swap.lastPaymentDate
            t_mat = (maturity_dt - self.value_dt) / g_days_in_year

            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self._dfs, df_mat)

            argtuple = (self, self.value_dt, swap)

            df_mat = optimize.newton(
                _f,
                x0=df_mat,
                fprime=None,
                args=argtuple,
                tol=SWAP_TOL,
                maxiter=50,
                fprime2=None,
                full_output=False,
            )

        if self.check_refit is True:
            self.check_refits(1e-10, SWAP_TOL, 1e-5)

    ###########################################################################

    def _build_curve_linear_swap_rate_interpolation(self):
        """Construct the discount curve using a bootstrap approach. This is
        the linear swap rate method that is fast and exact as it does not
        require the use of a solver. It is also market standard."""

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        t_mat = 0.0
        df_mat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, df_mat)

        for depo in self.used_deposits:
            df_settle = self.df(depo.start_dt)
            df_mat = depo.maturity_df() * df_settle
            t_mat = (depo.maturity_dt - self.value_dt) / g_days_in_year
            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self._dfs, df_mat)

        oldt_mat = t_mat

        for fra in self.used_fras:

            t_set = (fra.start_dt - self.value_dt) / g_days_in_year
            t_mat = (fra.maturity_dt - self.value_dt) / g_days_in_year

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if t_set < oldt_mat and t_mat > oldt_mat:
                df_mat = fra.maturity_df(self)
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)
            else:
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self._dfs, df_mat)

                argtuple = (self, self.value_dt, fra)
                df_mat = optimize.newton(
                    _g,
                    x0=df_mat,
                    fprime=None,
                    args=argtuple,
                    tol=SWAP_TOL,
                    maxiter=50,
                    fprime2=None,
                )

        if len(self.used_swaps) == 0:
            if self.check_refit is True:
                self.check_refits(1e-10, SWAP_TOL, 1e-5)
            return

        #        print("CURVE SO FAR")
        #        print(self._times)
        #        print(self._dfs)

        #######################################################################
        # ADD SWAPS TO CURVE
        #######################################################################

        # Find where the FRAs and Depos go up to as this bit of curve is done
        found_start = False
        last_dt = self.value_dt
        if len(self.used_deposits) != 0:
            last_dt = self.used_deposits[-1].maturity_dt

        if len(self.used_fras) != 0:
            last_dt = self.used_fras[-1].maturity_dt

        # We use the longest swap assuming it has a superset of ALL of the
        # swap flow dates used in the curve construction
        longest_swap = self.used_swaps[-1]
        cpn_dts = longest_swap.adjusted_fixed_dts
        num_flows = len(cpn_dts)

        # Find where first cpn without discount factor starts
        start_index = 0
        for i in range(0, num_flows):
            if cpn_dts[i] > last_dt:
                start_index = i
                found_start = True
                break

        if found_start is False:
            raise FinError("Found start is false. Swaps payments inside FRAs")

        swap_rates = []
        swap_times = []

        # I use the last cpn date for the swap rate interpolation as this
        # may be different from the maturity date due to a holiday adjustment
        # and the swap rates need to align with the cpn payment dates
        for swap in self.used_swaps:
            swap_rate = swap.fixed_cpn
            maturity_dt = swap.adjusted_fixed_dts[-1]
            tswap = (maturity_dt - self.value_dt) / g_days_in_year
            swap_times.append(tswap)
            swap_rates.append(swap_rate)

        interpolated_swap_rates = [0.0]
        interpolated_swap_times = [0.0]

        for dt in cpn_dts[1:]:
            swap_years = (dt - self.value_dt) / g_days_in_year
            swap_rate = np.interp(swap_years, swap_times, swap_rates)
            interpolated_swap_rates.append(swap_rate)
            interpolated_swap_times.append(swap_years)

        # Do I need this line ?
        interpolated_swap_rates[0] = interpolated_swap_rates[1]

        #        print("Interpolated swap times:", interpolated_swap_times)
        #        print("Interpolated swap rates:", interpolated_swap_rates)

        accrual_factors = longest_swap.fixed_year_fracs

        acc = 0.0
        df = 1.0
        pv01 = 0.0
        df_settle = self.df(longest_swap.start_dt)
        #        print("SETTLE", df_settle)

        for i in range(1, start_index):
            dt = cpn_dts[i]
            df = self.df(dt)
            acc = accrual_factors[i - 1]
            pv01 += acc * df
        #            print("BEFORE", i, dt, df, acc, pv01)

        for i in range(start_index, num_flows):

            dt = cpn_dts[i]
            t_mat = (dt - self.value_dt) / g_days_in_year
            swap_rate = interpolated_swap_rates[i]
            acc = accrual_factors[i - 1]
            pv01_end = acc * swap_rate + 1.0

            df_mat = (df_settle - swap_rate * pv01) / pv01_end

            #  print("IN: %12s %12.10f %12.10f %12.10f %12.10f OUT: %14.12f" %
            #                  (dt, swap_rate, acc, pv01, pv01_end, df_mat))

            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self._dfs, df_mat)

            pv01 += acc * df_mat

        #        print(self._times)
        #        print(self._dfs)

        if self.check_refit is True:
            self.check_refits(1e-10, SWAP_TOL, 1e-5)

    ###########################################################################

    def _check_refits(self, depo_tol, fra_tol, swap_tol):
        """Ensure that the Ibor curve refits the calibration instruments."""
        for depo in self.used_deposits:
            v = depo.value(self.value_dt, self) / depo.notional
            if abs(v - 1.0) > depo_tol:
                print("Value", v)
                raise FinError("Deposit not repriced.")

        for fra in self.used_fras:
            v = fra.value(self.value_dt, self) / fra.notional
            if abs(v) > fra_tol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self.used_swaps:
            # We value it as of the start date of the swap
            v = swap.value(swap.start_dt, self, self, None, principal=0.0)
            v = v / swap.notional
            if abs(v) > swap_tol:
                print(
                    "Swap with maturity "
                    + str(swap.maturity_dt)
                    + " Not Repriced. Has Value",
                    v,
                )
                swap.print_fixed_leg_pv()
                swap.print_float_leg_pv()
                raise FinError("Swap not repriced.")

    ###########################################################################

    def __repr__(self):
        """Print out the details of the Ibor curve."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUATION DATE", self.value_dt)

        for depo in self.used_deposits:
            s += label_to_string("DEPOSIT", "")
            s += depo.__repr__()

        for fra in self.used_fras:
            s += label_to_string("FRA", "")
            s += fra.__repr__()

        for swap in self.used_swaps:
            s += label_to_string("SWAP", "")
            s += swap.__repr__()

        num_points = len(self._times)

        s += label_to_string("INTERP TYPE", self._interp_type)

        s += label_to_string("GRID TIMES", "GRID DFS")
        for i in range(0, num_points):
            s += label_to_string(
                "% 10.6f" % self._times[i], "%12.10f" % self._dfs[i]
            )

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
