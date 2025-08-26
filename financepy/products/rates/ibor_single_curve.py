##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import copy
import numpy as np
from scipy import optimize

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.date import datediff
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types, _func_name
from ...utils.global_vars import G_DAYS_IN_YEARS
from ...market.curves.interpolator import InterpTypes, Interpolator
from ...market.curves.discount_curve import DiscountCurve
from ...products.rates.ibor_deposit import IborDeposit
from ...products.rates.ibor_fra import IborFRA
from ...products.rates.ibor_swap import IborSwap

SWAP_TOL = 1e-10

##############################################################################
# TODO: CHANGE times to df_times
##############################################################################


def _f(
    df,
    *args,
):
    """Root search objective function for IRS"""

    curve = args[0]
    value_dt = args[1]
    swap = args[2]
    curve.set_last_df(df)

    # For discount that need a fit function, we fit it now
    curve.fit(curve.times, curve.dfs)
    v_swap = swap.value(value_dt, curve, curve, None)
    notional = swap.fixed_leg.notional
    v_swap /= notional
    return v_swap


########################################################################################


def _g(df, *args):
    """Root search objective function for FRAs"""
    curve = args[0]
    value_dt = args[1]
    fra = args[2]
    curve.set_last_df(df)
    # For discount that need a fit function, we fit it now
    curve.fit(curve.times, curve.dfs)
    v_fra = fra.value(value_dt, curve)
    v_fra /= fra.notional
    return v_fra


########################################################################################


def _cost_function(dfs, *args):
    """Objective function for fitting all knot dfs at once to the benchmark
    securities  -- suitable for non-local interpolators"""

    #    print("Discount factors:", dfs)

    libor_curve = args[0]
    value_dt = libor_curve.value_dt
    libor_curve._dfs = dfs

    # For discount that need a fit function, we fit it now
    libor_curve.fit(libor_curve.times, libor_curve.dfs)

    cost = 0.0
    for depo in libor_curve.used_deposits:
        v = depo.value(value_dt, libor_curve) / depo.notional
        #        print("DEPO:", depo.maturity_dt, v)
        cost += (v - 1.0) ** 2

    for fra in libor_curve.used_fras:
        v = fra.value(value_dt, libor_curve) / fra.notional
        #        print("FRA:", fra.maturity_dt, v)
        cost += v * v

    for swap in libor_curve.used_swaps:
        v = (
            swap.value(value_dt, libor_curve)
            / swap.fixed_leg.notional
            / swap.pv01(value_dt, libor_curve)
        )
        #        print("SWAP:", swap.maturity_dt, v)
        cost += v * v

    print("Cost:", cost)
    return cost


########################################################################################


class IborSingleCurve(DiscountCurve):
    """Constructs one discount and index curve as implied by prices of Ibor
    deposits, FRAs and IRS. Discounting is assumed to be at Libor and the value
    of the floating leg (including a notional) is assumed to be par. This
    approach has been overtaken since 2008 as OIS discounting has become the
    agreed discounting approach for ISDA derivatives. This curve method is
    therefore intended for those happy to assume simple Libor discounting.

    The curve date is the date on which we are performing the valuation based
    on the information available on the curve date. Typically it is the date on
    which an amount of 1 unit paid has a present value of 1. This class
    inherits from DiscountCurve and so it has all of the methods that that
    class has.

    There are two main curve-building approaches:

    1) The first uses a bootstrap that interpolates swap rates linearly for
    coupon dates that fall between the swap maturity dates. With this, we can
    solve for the discount factors iteratively without need of a solver. This
    will give us a set of discount factors on the grid dates that refit the
    market exactly. However, when extracting discount factors, we will then
    assume flat forward rates between these coupon dates. There is no
    contradiction as it is as though we had been quoted a swap curve with all
    of the market swap rates, and with an additional set as though the market
    quoted swap rates at a higher frequency than the market.

    2) The second uses a bootstrap that uses only the swap rates provided but
    which also assumes that forwards are flat between these swap maturity
    dates. This approach is non-linear and so requires a solver. Consequently
    it is slower. Its advantage is that we can switch interpolation schemes
    to provide a smoother or other functional curve shape which may have a more
    economically justifiable shape. However the root search makes it slower."""

    ####################################################################################

    def __init__(
        self,
        value_dt: Date,  # This is the trade date (not T+2)
        ibor_deposits: list,
        ibor_fras: list,
        ibor_swaps: list,
        interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
        check_refit_flag: bool = False,  # Set to True to test it works
        do_build: bool = True,
        **kwargs,
    ):
        """Create an instance of a FinIbor curve given a valuation date and
        a set of ibor deposits, ibor FRAs and ibor_swaps. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on coupon dates and to then assume
        flat forwards between these coupon dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        If no instrument is starting on the valuation date, the curve is then
        assumed to be flat out to the first instrument using its zero rate.
        """

        check_argument_types(getattr(self, _func_name(), None), locals())

        self.value_dt = value_dt
        self._interp_type = interp_type
        self.check_refit_flag = check_refit_flag
        self._interpolator = Interpolator(self._interp_type, **kwargs)
        self.is_built = False
        self._optional_interp_params = kwargs

        self._validate_inputs(
            ibor_deposits=ibor_deposits,
            ibor_fras=ibor_fras,
            ibor_swaps=ibor_swaps,
        )

        if do_build:
            self.build_curve(**kwargs)

    ####################################################################################

    def build_curve(self, **kwargs):
        """
        Build curve based on interpolation.

        Not all interpolators are suitable for the boostrap/1d solver, only those
        that are local,
        where the value of df[i] does not affect discount factors for t<=t[i-1]
        """

        if Interpolator.suitable_for_bootstrap(self._interp_type):
            self._build_curve_using_1d_solver(**kwargs)
        else:
            self._build_curve_using_least_squares(**kwargs)

        self.is_built = True

    ####################################################################################

    def _validate_inputs(self, ibor_deposits, ibor_fras, ibor_swaps):
        """Validate the inputs for each of the Ibor products."""

        num_depos = len(ibor_deposits)
        num_fras = len(ibor_fras)
        num_swaps = len(ibor_swaps)

        depo_start_dt = self.value_dt
        swap_start_dt = self.value_dt

        if num_depos + num_fras + num_swaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if num_depos > 0:

            depo_start_dt = ibor_deposits[0].start_dt

            for depo in ibor_deposits:

                if isinstance(depo, IborDeposit) is False:
                    raise FinError("Deposit is not of type IborDeposit")

                start_dt = depo.start_dt

                if start_dt < self.value_dt:
                    raise FinError("First deposit starts before valuation date.")

                if start_dt < depo_start_dt:
                    depo_start_dt = start_dt

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

        # REMOVED THIS AS WE WANT TO ANCHOR CURVE AT VALUATION DATE
        # USE A SYNTHETIC DEPOSIT TO BRIDGE GAP FROM VALUE DATE TO SETTLEMENT DATE
        # Ensure that valuation date is on or after first deposit start date
        # if num_depos > 1:
        #    if ibor_deposits[0].effective_dt > self.value_dt:
        #     raise FinError("Valuation date must not be before first deposit settles.")

        if num_fras > 0:
            for fra in ibor_fras:
                if isinstance(fra, IborFRA) is False:
                    raise FinError("FRA is not of type IborFRA")

                start_dt = fra.start_dt
                if start_dt < self.value_dt:
                    raise FinError("FRAs starts before valuation date")

        if num_fras > 1:
            prev_dt = ibor_fras[0].maturity_dt
            for fra in ibor_fras[1:]:
                next_dt = fra.maturity_dt
                if next_dt <= prev_dt:
                    raise FinError("FRAs must be in increasing maturity")
                prev_dt = next_dt

        if num_swaps > 0:

            swap_start_dt = ibor_swaps[0].effective_dt

            for swap in ibor_swaps:

                # is False and isinstance(swap, IborSwap) is False:
                if isinstance(swap, IborSwap) is False:
                    raise FinError("Swap is not of type IborSwap")

                start_dt = swap.effective_dt
                if start_dt < self.value_dt:
                    raise FinError("Swaps starts before valuation date.")

                if swap.effective_dt < swap_start_dt:
                    swap.start_dt = swap.effective_dt

        if num_swaps > 1:

            # Swaps must all start on the same date for the bootstrap
            start_dt = ibor_swaps[0].effective_dt
            for swap in ibor_swaps[1:]:
                next_start_dt = swap.effective_dt
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

            longest_swap_cpn_dates = longest_swap.fixed_leg.payment_dts

            for swap in ibor_swaps[0:-1]:

                swap_cpn_dts = swap.fixed_leg.payment_dts

                num_flows = len(swap_cpn_dts)
                for i_flow in range(0, num_flows):
                    if swap_cpn_dts[i_flow] != longest_swap_cpn_dates[i_flow]:
                        raise FinError("Swap coupons are not on the same date grid.")

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

        first_swap_maturity_dt = None

        if num_swaps > 0:
            first_swap_maturity_dt = ibor_swaps[0].maturity_dt

        if num_depos > 0 and num_fras > 0:
            if first_fra_maturity_dt <= last_deposit_maturity_dt:
                print("FRA Maturity Date:", first_fra_maturity_dt)
                print("Last Deposit Date:", last_deposit_maturity_dt)
                raise FinError("First FRA must end after last Deposit")

        if num_fras > 0 and num_swaps > 0:
            if first_swap_maturity_dt <= last_fra_maturity_dt:
                raise FinError("First Swap must mature after last FRA ends")

        # If both depos and swaps start after T, we need a rate to get them to
        # the first deposit. So we create a synthetic deposit rate contract.

        if swap_start_dt > self.value_dt:

            if num_depos == 0:
                raise FinError("Need a deposit rate to pin down short end.")

            if depo_start_dt > self.value_dt:
                first_depo = ibor_deposits[0]
                if first_depo.start_dt > self.value_dt:
                    synthetic_deposit = copy.deepcopy(first_depo)
                    synthetic_deposit.start_dt = self.value_dt
                    synthetic_deposit.maturity_dt = first_depo.start_dt
                    ibor_deposits.insert(0, synthetic_deposit)
                    num_depos += 1

        # Now determine which instruments are used
        self.used_deposits = ibor_deposits
        self.used_fras = ibor_fras
        self.used_swaps = ibor_swaps

        #  Need the floating leg basis for the curve
        if len(self.used_swaps) > 0:
            self.dc_type = ibor_swaps[0].float_leg.dc_type
        else:
            self.dc_type = None

    ####################################################################################

    def _build_curve_using_1d_solver(self, **kwargs):
        """Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver."""

        self._interpolator = Interpolator(self._interp_type, **kwargs)
        self._times = np.array([])
        self._dfs = np.array([])
        self.is_built = True

        # time zero is now.
        t_mat = 0.0
        df_mat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self.dfs, df_mat)
        self.fit(self.times, self._dfs)

        for depo in self.used_deposits:

            df_settle_dt = self.df(depo.start_dt)
            df_mat = depo.maturity_df() * df_settle_dt
            t_mat = (depo.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS

            self._times = np.append(self._times, t_mat)
            self._dfs = np.append(self.dfs, df_mat)
            self.fit(self.times, self.dfs)

        oldt_mat = t_mat

        for fra in self.used_fras:

            t_set = (fra.start_dt - self.value_dt) / G_DAYS_IN_YEARS
            t_mat = (fra.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if t_set < oldt_mat and t_mat > oldt_mat:
                df_mat = fra.maturity_df(self)
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self.dfs, df_mat)
            else:
                self._times = np.append(self._times, t_mat)
                self._dfs = np.append(self.dfs, df_mat)
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
            # I use the last_payment_dt in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturity_dt = swap.fixed_leg.payment_dts[-1]
            t_mat = (maturity_dt - self.value_dt) / G_DAYS_IN_YEARS

            self._times = np.append(self.times, t_mat)
            self._dfs = np.append(self.dfs, df_mat)

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

        if self.check_refit_flag is True:
            # self.check_refit(1e-10, swaptol, 1e-5)
            self.check_refit(1e-5, 1e-5, 1e-5)

    ####################################################################################

    def _build_curve_using_least_squares(self, **kwargs):
        """
        Construct the discount curve using a least-squares minimisation approach.
        This enables a more complex interpolation scheme.
        """

        def _obj_f_curve_build_ls(dfs):
            """
            Objective function for fitting all knot dfs at once to the benchmark
            securities  -- suitable for non-local interpolators
            """

            libor_curve = self
            value_dt = libor_curve.value_dt
            libor_curve._dfs[1:] = dfs

            libor_curve.fit(libor_curve.times, libor_curve.dfs)

            out = np.zeros(
                len(libor_curve.used_deposits)
                + len(libor_curve.used_fras)
                + len(libor_curve.used_swaps)
            )

            idx = 0
            for depo in libor_curve.used_deposits:
                # do not need to be too exact here
                acc_factor = datediff(depo.start_dt, depo.maturity_dt)
                # as rate
                r = (
                    -np.log(depo.value(value_dt, libor_curve) / depo.notional)
                    / acc_factor
                )
                out[idx] = r
                idx = idx + 1

            for fra in libor_curve.used_fras:
                # do not need to be too exact here
                acc_factor = datediff(fra.start_dt, fra.maturity_dt)
                v = fra.value(value_dt, libor_curve) / fra.notional / acc_factor
                out[idx] = v
                idx = idx + 1

            for swap in libor_curve.used_swaps:
                v = (
                    swap.value(value_dt, libor_curve)
                    / swap.fixed_leg.notional
                    / swap.pv01(value_dt, libor_curve)
                )
                out[idx] = v
                idx = idx + 1

            return out

        bootstrap_first = True
        self.is_built = True

        if bootstrap_first:
            orig_check_refit_flag = self.check_refit_flag
            self.check_refit_flag = False
            self._build_curve_using_1d_solver(**kwargs)
            self.check_refit_flag = orig_check_refit_flag
        else:
            tmat = 0.0
            df_mat = 1.0

            grid_times = [tmat]
            grid_dfs = [df_mat]

            self._interpolator = Interpolator(self._interp_type, **kwargs)

            for depo in self.used_deposits:
                tmat = (depo.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
                grid_times.append(tmat)

            for depo in self.used_deposits:
                tmat = (depo.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
                grid_times.append(tmat)

            for fra in self.used_fras:
                tmat = (fra.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
                grid_times.append(tmat)
                grid_dfs.append(df_mat)

            for swap in self.used_swaps:
                tmat = (swap.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
                grid_times.append(tmat)

            self._times = np.array(grid_times)
            self._dfs = np.exp(-self._times * 0.05)

        res = optimize.least_squares(
            _obj_f_curve_build_ls,
            self._dfs[1:],
            bounds=(0, np.inf),
            ftol=1e-4,
            xtol=1e-6,
        )

        self._dfs[1:] = np.array(res.x)
        self.fit(self._times, self._dfs)

        if self.check_refit_flag is True:
            self.check_refit(1e-5, 1e-5, 1e-5)

    ####################################################################################

    def _build_curve_linear_swap_rate_interpolation(self):
        """Construct the discount curve using a bootstrap approach. This is
        the linear swap rate method that is fast and exact as it does not
        require the use of a solver. It is also market standard."""

        self._interpolator = Interpolator(self._interp_type)
        self.is_built = True

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        t_mat = 0.0
        df_mat = 1.0
        self._times = np.append(self.times, 0.0)
        self._dfs = np.append(self.dfs, df_mat)
        self.fit(self.times, self.dfs)

        for depo in self.used_deposits:
            df_settle_dt = self.df(depo.start_dt)
            df_mat = depo.maturity_df() * df_settle_dt
            t_mat = (depo.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
            self._times = np.append(self.times, t_mat)
            self._dfs = np.append(self.dfs, df_mat)
            self.fit(self.times, self.dfs)

        oldt_mat = t_mat

        for fra in self.used_fras:

            t_set = (fra.start_dt - self.value_dt) / G_DAYS_IN_YEARS
            t_mat = (fra.maturity_dt - self.value_dt) / G_DAYS_IN_YEARS

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if t_set < oldt_mat and t_mat > oldt_mat:

                df_mat = fra.maturity_df(self)

                self._times = np.append(self.times, t_mat)
                self._dfs = np.append(self.dfs, df_mat)
                self.fit(self.times, self._dfs)

            else:

                self._times = np.append(self.times, t_mat)
                self._dfs = np.append(self.dfs, df_mat)
                self.fit(self.times, self.dfs)

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
            if self.check_refit_flag is True:
                self.check_refit(1e-10, SWAP_TOL, 1e-5)
            return

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

        # Find where first coupon without discount factor starts
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

        # I use the last coupon date for the swap rate interpolation as this
        # may be different from the maturity date due to a holiday adjustment
        # and the swap rates need to align with the coupon payment dates
        for swap in self.used_swaps:
            swap_rate = swap.fixed_leg.coupon
            maturity_dt = swap.adjusted_fixed_dts[-1]
            tswap = (maturity_dt - self.value_dt) / G_DAYS_IN_YEARS
            swap_times.append(tswap)
            swap_rates.append(swap_rate)

        interpolated_swap_rates = [0.0]
        interpolated_swap_times = [0.0]

        for dt in cpn_dts[1:]:
            swap_years = (dt - self.value_dt) / G_DAYS_IN_YEARS
            swap_rate = np.interp(swap_years, swap_times, swap_rates)
            interpolated_swap_rates.append(swap_rate)
            interpolated_swap_times.append(swap_years)

        # Do I need this line ?
        interpolated_swap_rates[0] = interpolated_swap_rates[1]

        accrual_factors = longest_swap.fixed_year_fracs

        acc = 0.0
        df = 1.0
        pv01 = 0.0
        df_settle_dt = self.df(longest_swap.effective_dt)

        for i in range(1, start_index):
            dt = cpn_dts[i]
            df = self.df(dt)
            acc = accrual_factors[i - 1]
            pv01 += acc * df

        for i in range(start_index, num_flows):

            dt = cpn_dts[i]
            t_mat = (dt - self.value_dt) / G_DAYS_IN_YEARS
            swap_rate = interpolated_swap_rates[i]
            acc = accrual_factors[i - 1]
            pv01_end = acc * swap_rate + 1.0

            df_mat = (df_settle_dt - swap_rate * pv01) / pv01_end

            self._times = np.append(self.times, t_mat)
            self._dfs = np.append(self.dfs, df_mat)
            self._interpolator.fit(self.times, self.dfs)

            pv01 += acc * df_mat

        if self.check_refit_flag is True:
            self.check_refit(1e-10, SWAP_TOL, 1e-5)

    ####################################################################################

    def check_refit(self, depo_tol, fra_tol, swap_tol):
        """Ensure that the Ibor curve refits the calibration instruments."""
        for depo in self.used_deposits:
            v = depo.value(self.value_dt, self) / depo.notional
            if abs(v - 1.0) > depo_tol:
                raise FinError(
                    f"Deposit not repriced, error = {abs(v - 1.0)} vs tol={depo_tol}"
                )

        for fra in self.used_fras:
            v = fra.value(self.value_dt, self, self) / fra.notional
            if abs(v) > fra_tol:
                raise FinError(f"FRA not repriced, error = {abs(v) } vs tol={fra_tol}")

        for swap in self.used_swaps:
            # We value it as of the start date of the swap
            v = swap.value(swap.effective_dt, self, self, None)
            v = (
                v / swap.fixed_leg.notional / swap.pv01(self.value_dt, self)
            )  # express in terms of the rate
            #            print("REFIT SWAP VALUATION:", swap._adjusted_maturity_dt, v)
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

    ####################################################################################

    def _df(self, t: Union[float, np.ndarray]):
        """
        Override from DiscountCurve so we can check if the curve
        has actually been built.
        """
        assert (
            self.is_built
        ), "The curve has not yet been built, call build_curve() first"
        return super().df_t(t)

    ####################################################################################

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

        s += label_to_string("INTERP TYPE", self._interp_type)
        s += label_to_string("IS BUILT", self.is_built)

        if self.is_built:
            num_points = len(self._times)
            s += label_to_string("GRID TIMES", "GRID DFS")
            for i in range(num_points):
                s += label_to_string(
                    f"{self._times[i]:10.6f}", f"{self._dfs[i]:12.10f}"
                )

        return s

    ####################################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
