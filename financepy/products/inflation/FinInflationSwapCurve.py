##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types, _func_name
from ...utils.global_vars import gDaysInYear
from ...market.curves.interpolator import InterpTypes
from ...market.curves.discount_curve import DiscountCurve

swaptol = 1e-8

##############################################################################
# TODO: CHANGE times to df_times
##############################################################################


def _f(df, *args):
    """ Root search objective function for swaps """
    curve = args[0]
    valuation_date = args[1]
    swap = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df
    v_swap = swap.value(valuation_date, curve, curve, None, 1.0)
    v_swap /= swap._notional
    return v_swap

###############################################################################


def _g(df, *args):
    """ Root search objective function for swaps """
    curve = args[0]
    valuation_date = args[1]
    fra = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df
    v_fra = fra.value(valuation_date, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class InflationSwapCurve(DiscountCurve):
    """ Constructs a discount curve as implied by the prices of Ibor
    deposits, FRAs and IRS. The curve date is the date on which we are
    performing the valuation based on the information available on the
    curve date. Typically it is the date on which an amount of 1 unit paid
    has a present value of 1. This class inherits from FinDiscountCurve
    and so it has all of the methods that that class has.

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

###############################################################################

    def __init__(self,
                 valuation_date: Date,
                 ibor_deposits: list,
                 ibor_fras: list,
                 ibor_swaps: list,
                 interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
                 check_refit: bool = False):  # Set to True to test it works
        """ Create an instance of a Ibor curve given a valuation date and
        a set of ibor deposits, ibor FRAs and ibor_swaps. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on coupon dates and to then assume
        flat forwards between these coupon dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        """
        raise FinError("DO NOT USE")

        check_argument_types(getattr(self, _func_name(), None), locals())

        self._valuation_date = valuation_date
        self._validate_inputs(ibor_deposits, ibor_fras, ibor_swaps)
        self._interp_type = interp_type
        self._check_refit = check_refit
        self._build_curve()

###############################################################################

    def _build_curve(self):
        """ Build curve based on interpolation. """
        self._build_curve_using_solver()

###############################################################################

    def _validate_inputs(self,
                         ibor_deposits,
                         ibor_fras,
                         ibor_swaps):
        """ Validate the inputs for each of the Ibor products. """

        num_depos = len(ibor_deposits)
        num_fras = len(ibor_fras)
        num_swaps = len(ibor_swaps)

        if num_depos + num_fras + num_swaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if num_depos > 0:
            for depo in ibor_deposits:
                startDt = depo._start_date
                if startDt < self._valuation_date:
                    raise FinError("First deposit starts before value date.")

            for depo in ibor_deposits:
                startDt = depo._start_date
                endDt = depo._maturity_date
                if startDt >= endDt:
                    raise FinError("First deposit ends on or before it begins")

        # Ensure order of depos
        if num_depos > 1:
            prev_dt = ibor_deposits[0]._maturity_date

            for depo in ibor_deposits[1:]:
                next_dt = depo._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("Deposits must be in increasing maturity")
                prev_dt = next_dt

        # Ensure that valuation date is on or after first deposit start date
        if num_depos > 1:
            if ibor_deposits[0]._start_date > self._valuation_date:
                raise FinError(
                    "Valuation date must not be before first deposit settles.")

        if num_fras > 0:
            for fra in ibor_fras:
                startDt = fra._start_date
                if startDt <= self._valuation_date:
                    raise FinError("FRAs starts before valuation date")

        if num_fras > 1:
            prev_dt = ibor_fras[0]._maturity_date
            for fra in ibor_fras[1:]:
                next_dt = fra._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("FRAs must be in increasing maturity")
                prev_dt = next_dt

        if num_swaps > 0:
            for swap in ibor_swaps:
                startDt = swap._start_date
                if startDt < self._valuation_date:
                    raise FinError("Swaps starts before valuation date.")

        if num_swaps > 1:

            # Swaps must all start on the same date for the bootstrap
            startDt = ibor_swaps[0]._start_date
            for swap in ibor_swaps[1:]:
                nextStartDt = swap._start_date
                if nextStartDt != startDt:
                    raise FinError("Swaps must all have same start date.")

            # Swaps must be increasing in tenor/maturity
            prev_dt = ibor_swaps[0]._maturity_date
            for swap in ibor_swaps[1:]:
                next_dt = swap._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("Swaps must be in increasing maturity")
                prev_dt = next_dt

            # Swaps must have same cash flows for bootstrap to work
            longestSwap = ibor_swaps[-1]
            longestSwapCpnDates = longestSwap._adjustedFixedDates
            for swap in ibor_swaps[0:-1]:
                swapCpnDates = swap._adjustedFixedDates
                num_flows = len(swapCpnDates)
                for iFlow in range(0, num_flows):
                    if swapCpnDates[iFlow] != longestSwapCpnDates[iFlow]:
                        raise FinError(
                            "Swap coupons are not on the same date grid.")

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        lastDepositMaturityDate = Date(1, 1, 1900)
        firstFRAMaturityDate = Date(1, 1, 1900)
        lastFRAMaturityDate = Date(1, 1, 1900)

        if num_depos > 0:
            lastDepositMaturityDate = ibor_deposits[-1]._maturity_date

        if num_fras > 0:
            firstFRAMaturityDate = ibor_fras[0]._maturity_date
            lastFRAMaturityDate = ibor_fras[-1]._maturity_date

        if num_swaps > 0:
            firstSwapMaturityDate = ibor_swaps[0]._maturity_date

        if num_depos > 0 and num_fras > 0:
            if firstFRAMaturityDate <= lastDepositMaturityDate:
                print("FRA Maturity Date:", firstFRAMaturityDate)
                print("Last Deposit Date:", lastDepositMaturityDate)
                raise FinError("First FRA must end after last Deposit")

        if num_fras > 0 and num_swaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        # Now determine which instruments are used
        self._usedDeposits = ibor_deposits
        self._usedFRAs = ibor_fras
        self._usedSwaps = ibor_swaps
        self._day_count_type = None

###############################################################################

    def _build_curve_using_solver(self):
        """ Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver. """

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, dfMat)

        for depo in self._usedDeposits:
            dfSettle = self.df(depo._start_date)
            dfMat = depo._maturity_df() * dfSettle
            tmat = (depo._maturity_date - self._valuation_date) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._start_date - self._valuation_date) / gDaysInYear
            tmat = (fra._maturity_date - self._valuation_date) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturity_df(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)

                argtuple = (self, self._valuation_date, fra)
                dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        for swap in self._usedSwaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturity_date = swap._lastPaymentDate
            tmat = (maturity_date - self._valuation_date) / gDaysInYear

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

            argtuple = (self, self._valuation_date, swap)

            dfMat = optimize.newton(_f, x0=dfMat, fprime=None, args=argtuple,
                                    tol=swaptol, maxiter=50, fprime2=None,
                                    full_output=False)

        if self._check_refit is True:
            self._check_refits(1e-10, swaptol, 1e-5)

###############################################################################

    def _build_curve_linear_swap_rate_interpolation(self):
        """ Construct the discount curve using a bootstrap approach. This is
        the linear swap rate method that is fast and exact as it does not
        require the use of a solver. It is also market standard. """

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, dfMat)

        for depo in self._usedDeposits:
            dfSettle = self.df(depo._start_date)
            dfMat = depo._maturity_df() * dfSettle
            tmat = (depo._maturity_date - self._valuation_date) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._start_date - self._valuation_date) / gDaysInYear
            tmat = (fra._maturity_date - self._valuation_date) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturity_df(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)

                argtuple = (self, self._valuation_date, fra)
                dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        if len(self._usedSwaps) == 0:
            if self._check_refit is True:
                self._check_refits(1e-10, swaptol, 1e-5)
            return

#        print("CURVE SO FAR")
#        print(self._times)
#        print(self._dfs)

        #######################################################################
        # ADD SWAPS TO CURVE
        #######################################################################

        # Find where the FRAs and Depos go up to as this bit of curve is done
        foundStart = False
        lastDate = self._valuation_date
        if len(self._usedDeposits) != 0:
            lastDate = self._usedDeposits[-1]._maturity_date

        if len(self._usedFRAs) != 0:
            lastDate = self._usedFRAs[-1]._maturity_date

        # We use the longest swap assuming it has a superset of ALL of the
        # swap flow dates used in the curve construction
        longestSwap = self._usedSwaps[-1]
        couponDates = longestSwap._adjustedFixedDates
        num_flows = len(couponDates)

        # Find where first coupon without discount factor starts
        start_index = 0
        for i in range(0, num_flows):
            if couponDates[i] > lastDate:
                start_index = i
                foundStart = True
                break

        if foundStart is False:
            raise FinError("Found start is false. Swaps payments inside FRAs")

        swap_rates = []
        swapTimes = []

        # I use the last coupon date for the swap rate interpolation as this
        # may be different from the maturity date due to a holiday adjustment
        # and the swap rates need to align with the coupon payment dates
        for swap in self._usedSwaps:
            swap_rate = swap._fixed_coupon
            maturity_date = swap._adjustedFixedDates[-1]
            tswap = (maturity_date - self._valuation_date) / gDaysInYear
            swapTimes.append(tswap)
            swap_rates.append(swap_rate)

        interpolatedSwapRates = [0.0]
        interpolatedSwapTimes = [0.0]

        for dt in couponDates[1:]:
            swapTime = (dt - self._valuation_date) / gDaysInYear
            swap_rate = np.interp(swapTime, swapTimes, swap_rates)
            interpolatedSwapRates.append(swap_rate)
            interpolatedSwapTimes.append(swapTime)

        # Do I need this line ?
        interpolatedSwapRates[0] = interpolatedSwapRates[1]

#        print("Interpolated swap times:", interpolatedSwapTimes)
#        print("Interpolated swap rates:", interpolatedSwapRates)

        accrual_factors = longestSwap._fixed_year_fracs

        acc = 0.0
        df = 1.0
        pv01 = 0.0
        dfSettle = self.df(longestSwap._start_date)
#        print("SETTLE", dfSettle)

        for i in range(1, start_index):
            dt = couponDates[i]
            df = self.df(dt)
            acc = accrual_factors[i-1]
            pv01 += acc * df
#            print("BEFORE", i, dt, df, acc, pv01)

        for i in range(start_index, num_flows):

            dt = couponDates[i]
            tmat = (dt - self._valuation_date) / gDaysInYear
            swap_rate = interpolatedSwapRates[i]
            acc = accrual_factors[i-1]
            pv01End = (acc * swap_rate + 1.0)

            dfMat = (dfSettle - swap_rate * pv01) / pv01End

#            print("IN: %12s %12.10f %12.10f %12.10f %12.10f OUT: %14.12f" %
#                  (dt, swap_rate, acc, pv01, pv01End, dfMat))

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

            pv01 += acc * dfMat

#        print(self._times)
#        print(self._dfs)

        if self._check_refit is True:
            self._check_refits(1e-10, swaptol, 1e-5)

###############################################################################

    def _check_refits(self, depoTol, fraTol, swapTol):
        """ Ensure that the Ibor curve refits the calibration instruments. """
        for depo in self._usedDeposits:
            v = depo.value(self._valuation_date, self) / depo._notional
            if abs(v - 1.0) > depoTol:
                print("Value", v)
                raise FinError("Deposit not repriced.")

        for fra in self._usedFRAs:
            v = fra.value(self._valuation_date, self) / fra._notional
            if abs(v) > fraTol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self._usedSwaps:
            # We value it as of the start date of the swap
            v = swap.value(swap._start_date, self, self, None, principal=0.0)
            v = v / swap._notional
#            print("REFIT SWAP VALUATION:", swap._adjustedMaturityDate, v)
            if abs(v) > swapTol:
                print("Swap with maturity " + str(swap._maturity_date)
                      + " Not Repriced. Has Value", v)
                swap.print_fixed_leg_pv()
                swap.print_float_leg_pv()
                raise FinError("Swap not repriced.")

###############################################################################

    def __repr__(self):
        """ Print out the details of the Ibor curve. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUATION DATE", self._valuation_date)

        for depo in self._usedDeposits:
            s += label_to_string("DEPOSIT", "")
            s += depo.__repr__()

        for fra in self._usedFRAs:
            s += label_to_string("FRA", "")
            s += fra.__repr__()

        for swap in self._usedSwaps:
            s += label_to_string("SWAP", "")
            s += swap.__repr__()

        num_points = len(self._times)

        s += label_to_string("INTERP TYPE", self._interp_type)

        s += label_to_string("GRID TIMES", "GRID DFS")
        for i in range(0, num_points):
            s += label_to_string("% 10.6f" % self._times[i],
                                 "%12.10f" % self._dfs[i])

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
