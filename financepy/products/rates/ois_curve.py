##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
import copy

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.helpers import label_to_string
from ...utils.helpers import check_argument_types, _func_name
from ...utils.global_vars import gDaysInYear
from ...market.curves.interpolator import InterpTypes, Interpolator
from ...market.curves.discount_curve import DiscountCurve

from ...products.rates.ibor_deposit import IborDeposit
from ...products.rates.ois import OIS

swaptol = 1e-10

##############################################################################
# TODO: CHANGE times to df_times
##############################################################################


def _fois(oir, *args):
    """ Extract the implied overnight index rate assuming it is flat over
    period in question. """

    target_ois_rate = args[0]
    day_counter = args[1]
    date_schedule = args[2]

    start_date = date_schedule[0]
    end_date = date_schedule[-1]

    df = 1.0
    prev_dt = date_schedule[0]
    for dt in date_schedule[1:]:
        year_frac = day_counter.year_frac(prev_dt, dt)
        df = df * (1.0 + oir * year_frac)

    period = day_counter.year_frac(start_date, end_date)

    ois_rate = (df - 1.0) / period
    diff = ois_rate - target_ois_rate
    return diff

###############################################################################


def _f(df, *args):
    """ Root search objective function for OIS """

    curve = args[0]
    value_date = args[1]
    swap = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df

    # For discount that need a fit function, we fit it now
    curve._interpolator.fit(curve._times, curve._dfs)
    v_swap = swap.value(value_date, curve, None)
    notional = swap._fixed_leg._notional
    v_swap /= notional
    return v_swap

###############################################################################


def _g(df, *args):
    """ Root search objective function for swaps """
    curve = args[0]
    value_date = args[1]
    fra = args[2]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df

    # For discount that need a fit function, we fit it now
    curve._interpolator.fit(curve._times, curve._dfs)
    v_fra = fra.value(value_date, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class OISCurve(DiscountCurve):
    """ Constructs a discount curve as implied by the prices of Overnight
    Index Rate swaps. The curve date is the date on which we are
    performing the valuation based on the information available on the
    curve date. Typically it is the date on which an amount of 1 unit paid
    has a present value of 1. This class inherits from FinDiscountCurve
    and so it has all of the methods that that class has.

    The construction of the curve is assumed to depend on just the OIS curve,
    i.e. it does not include information from Ibor-OIS basis swaps. For this
    reason I call it a one-curve.
    """

###############################################################################

    def __init__(self,
                 value_date: Date,
                 ois_deposits: list,
                 ois_fras: list,
                 ois_swaps: list,
                 interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES,
                 check_refit: bool = False):  # Set to True to test it works
        """ Create an instance of an overnight index rate swap curve given a
        valuation date and a set of OIS rates. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on coupon dates and to then assume
        flat forwards between these coupon dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        """

        check_argument_types(getattr(self, _func_name(), None), locals())

        self._value_date = value_date
        self._validate_inputs(ois_deposits, ois_fras, ois_swaps)
        self._interp_type = interp_type
        self._check_refit = check_refit
        self._interpolator = None
        self._build_curve()

###############################################################################

    def _build_curve(self):
        """ Build curve based on interpolation. """

        self._build_curve_using_1d_solver()

###############################################################################

    def _validate_inputs(self,
                         oisDeposits,
                         oisFRAs,
                         oisSwaps):
        """ Validate the inputs for each of the Libor products. """

        num_depos = len(oisDeposits)
        num_fras = len(oisFRAs)
        num_swaps = len(oisSwaps)

        depo_start_date = self._value_date
        swap_start_date = self._value_date

        if num_depos + num_fras + num_swaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if num_depos > 0:

            depo_start_date = oisDeposits[0]._start_date

            for depo in oisDeposits:

                if isinstance(depo, IborDeposit) is False:
                    raise FinError("Deposit is not of type IborDeposit")

                start_date = depo._start_date

                if start_date < self._value_date:
                    raise FinError("First deposit starts before value date.")

                if start_date < depo_start_date:
                    depo_start_date = start_date

            for depo in oisDeposits:
                startDt = depo._start_date
                endDt = depo._maturity_date
                if startDt >= endDt:
                    raise FinError("First deposit ends on or before it begins")

        # Ensure order of depos
        if num_depos > 1:

            prev_dt = oisDeposits[0]._maturity_date
            for depo in oisDeposits[1:]:
                next_dt = depo._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("Deposits must be in increasing maturity")
                prev_dt = next_dt

        # REMOVED THIS AS WE WANT TO ANCHOR CURVE AT VALUATION DATE
        # USE A SYNTHETIC DEPOSIT TO BRIDGE GAP FROM VALUE DATE TO SETTLEMENT DATE
        # Ensure that valuation date is on or after first deposit start date
        # Ensure that valuation date is on or after first deposit start date
        # if num_depos > 1:
        #    if oisDeposits[0]._start_date > self._value_date:
        #        raise FinError("Valuation date must not be before first deposit settles.")

        # Validation of the inputs.
        if num_fras > 0:
            for fra in oisFRAs:
                startDt = fra._start_date
                if startDt <= self._value_date:
                    raise FinError("FRAs starts before valuation date")

                startDt = fra._start_date
                if startDt < self._value_date:
                    raise FinError("FRAs starts before valuation date")

        if num_fras > 1:
            prev_dt = oisFRAs[0]._maturity_date
            for fra in oisFRAs[1:]:
                next_dt = fra._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("FRAs must be in increasing maturity")
                prev_dt = next_dt

        if num_swaps > 0:

            swap_start_date = oisSwaps[0]._effective_date

            for swap in oisSwaps:

                if isinstance(swap, OIS) is False:
                    raise FinError("Swap is not of type FinOIS")

                startDt = swap._effective_date
                if startDt < self._value_date:
                    raise FinError("Swaps starts before valuation date.")

                if swap._effective_date < swap_start_date:
                    swap_start_date = swap._effective_date

        if num_swaps > 1:

            # Swaps must be increasing in tenor/maturity
            prev_dt = oisSwaps[0]._maturity_date
            for swap in oisSwaps[1:]:
                next_dt = swap._maturity_date
                if next_dt <= prev_dt:
                    raise FinError("Swaps must be in increasing maturity")
                prev_dt = next_dt

        # TODO: REINSTATE THESE CHECKS ?
            # Swaps must have same cash flows for linear swap bootstrap to work
#            longestSwap = oisSwaps[-1]
#            longestSwapCpnDates = longestSwap._adjusted_fixed_dates
#            for swap in oisSwaps[0:-1]:
#                swapCpnDates = swap._adjusted_fixed_dates
#                num_flows = len(swapCpnDates)
#                for iFlow in range(0, num_flows):
#                    if swapCpnDates[iFlow] != longestSwapCpnDates[iFlow]:
#                        raise FinError("Swap coupons are not on the same date grid.")

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        lastDepositMaturityDate = Date(1, 1, 1900)
        firstFRAMaturityDate = Date(1, 1, 1900)
        lastFRAMaturityDate = Date(1, 1, 1900)

        if num_depos > 0:
            lastDepositMaturityDate = oisDeposits[-1]._maturity_date

        if num_fras > 0:
            firstFRAMaturityDate = oisFRAs[0]._maturity_date
            lastFRAMaturityDate = oisFRAs[-1]._maturity_date

        if num_swaps > 0:
            first_swap_maturity_date = oisSwaps[0]._maturity_date

        if num_fras > 0 and num_swaps > 0:
            if first_swap_maturity_date <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA ends")

        if num_depos > 0 and num_fras > 0:
            if firstFRAMaturityDate <= lastDepositMaturityDate:
                print("FRA Maturity Date:", firstFRAMaturityDate)
                print("Last Deposit Date:", lastDepositMaturityDate)
                raise FinError("First FRA must end after last Deposit")

        if num_fras > 0 and num_swaps > 0:
            if first_swap_maturity_date <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        if swap_start_date > self._value_date:

            if num_depos == 0:
                raise FinError("Need a deposit rate to pin down short end.")

            if depo_start_date > self._value_date:
                firstDepo = oisDeposits[0]
                if firstDepo._start_date > self._value_date:
                    print("Inserting synthetic deposit")
                    syntheticDeposit = copy.deepcopy(firstDepo)
                    syntheticDeposit._start_date = self._value_date
                    syntheticDeposit._maturity_date = firstDepo._start_date
                    oisDeposits.insert(0, syntheticDeposit)
                    num_depos += 1

        # Now determine which instruments are used
        self._usedDeposits = oisDeposits
        self._usedFRAs = oisFRAs
        self._usedSwaps = oisSwaps

       # Need the floating leg basis for the curve
        if len(self._usedSwaps) > 0:
            self._dc_type = oisSwaps[0]._float_leg._dc_type
        else:
            self._dc_type = None

###############################################################################

    def _build_curve_using_1d_solver(self):
        """ Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver. """

        self._interpolator = Interpolator(self._interp_type)
        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        df_mat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, df_mat)
        self._interpolator.fit(self._times, self._dfs)

        for depo in self._usedDeposits:
            df_settle = self.df(depo._start_date)
            df_mat = depo._maturity_df() * df_settle
            tmat = (depo._maturity_date - self._value_date) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, df_mat)
            self._interpolator.fit(self._times, self._dfs)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._start_date - self._value_date) / gDaysInYear
            tmat = (fra._maturity_date - self._value_date) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                df_mat = fra.maturity_df(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, df_mat)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, df_mat)
                argtuple = (self, self._value_date, fra)
                df_mat = optimize.newton(_g, x0=df_mat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        for swap in self._usedSwaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturity_date = swap._fixed_leg._payment_dates[-1]
            tmat = (maturity_date - self._value_date) / gDaysInYear

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, df_mat)

            argtuple = (self, self._value_date, swap)

            df_mat = optimize.newton(_f, x0=df_mat, fprime=None, args=argtuple,
                                    tol=swaptol, maxiter=50, fprime2=None,
                                    full_output=False)

        if self._check_refit is True:
            self._check_refits(1e-10, swaptol, 1e-5)

###############################################################################

    def _build_curve_linear_swap_rate_interpolation(self):
        """ Construct the discount curve using a bootstrap approach. This is
        the linear swap rate method that is fast and exact as it does not
        require the use of a solver. It is also market standard. """

        self._interpolator = Interpolator(self._interp_type)

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        df_mat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, df_mat)

        for depo in self._usedDeposits:
            df_settle = self.df(depo._start_date)
            df_mat = depo._maturity_df() * df_settle
            tmat = (depo._maturity_date - self._value_date) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, df_mat)
            self._interpolator.fit(self._times, self._dfs)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._start_date - self._value_date) / gDaysInYear
            tmat = (fra._maturity_date - self._value_date) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                df_mat = fra.maturity_df(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, df_mat)
                self._interpolator.fit(self._times, self._dfs)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, df_mat)
                self._interpolator.fit(self._times, self._dfs)

                argtuple = (self, self._value_date, fra)
                df_mat = optimize.newton(_g, x0=df_mat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        if len(self._usedSwaps) == 0:
            if self._check_refit is True:
                self._check_refits(1e-10, swaptol, 1e-5)
            return

        #######################################################################
        # ADD SWAPS TO CURVE
        #######################################################################

        # Find where the FRAs and Depos go up to as this bit of curve is done
        foundStart = False
        lastDate = self._value_date
        if len(self._usedDeposits) != 0:
            lastDate = self._usedDeposits[-1]._maturity_date

        if len(self._usedFRAs) != 0:
            lastDate = self._usedFRAs[-1]._maturity_date

        # We use the longest swap assuming it has a superset of ALL of the
        # swap flow dates used in the curve construction
        longestSwap = self._usedSwaps[-1]
        cpn_dates = longestSwap._adjusted_fixed_dates
        num_flows = len(cpn_dates)

        # Find where first coupon without discount factor starts
        start_index = 0
        for i in range(0, num_flows):
            if cpn_dates[i] > lastDate:
                start_index = i
                foundStart = True
                break

        if foundStart is False:
            raise FinError("Found start is false. Swaps payments inside FRAs")

        swap_rates = []
        swap_times = []

        # I use the last coupon date for the swap rate interpolation as this
        # may be different from the maturity date due to a holiday adjustment
        # and the swap rates need to align with the coupon payment dates
        for swap in self._usedSwaps:
            swap_rate = swap._fixed_coupon
            maturity_date = swap._adjusted_fixed_dates[-1]
            tswap = (maturity_date - self._value_date) / gDaysInYear
            swap_times.append(tswap)
            swap_rates.append(swap_rate)

        interpolatedSwapRates = [0.0]
        interpolatedswap_times = [0.0]

        for dt in cpn_dates[1:]:
            swapTime = (dt - self._value_date) / gDaysInYear
            swap_rate = np.interp(swapTime, swap_times, swap_rates)
            interpolatedSwapRates.append(swap_rate)
            interpolatedswap_times.append(swapTime)

        # Do I need this line ?
        interpolatedSwapRates[0] = interpolatedSwapRates[1]

        accrual_factors = longestSwap._fixed_year_fracs

        acc = 0.0
        df = 1.0
        pv01 = 0.0
        df_settle = self.df(longestSwap._start_date)

        for i in range(1, start_index):
            dt = cpn_dates[i]
            df = self.df(dt)
            acc = accrual_factors[i-1]
            pv01 += acc * df

        for i in range(start_index, num_flows):

            dt = cpn_dates[i]
            tmat = (dt - self._value_date) / gDaysInYear
            swap_rate = interpolatedSwapRates[i]
            acc = accrual_factors[i-1]
            pv01End = (acc * swap_rate + 1.0)

            df_mat = (df_settle - swap_rate * pv01) / pv01End

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, df_mat)
            self._interpolator.fit(self._times, self._dfs)

            pv01 += acc * df_mat

        if self._check_refit is True:
            self._check_refits(1e-10, swaptol, 1e-5)

###############################################################################

    def _check_refits(self, depoTol, fraTol, swapTol):
        """ Ensure that the Libor curve refits the calibration instruments. """

        for fra in self._usedFRAs:
            v = fra.value(self._value_date, self) / fra._notional
            if abs(v) > fraTol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self._usedSwaps:
            # We value it as of the start date of the swap
            v = swap.value(swap._effective_date, self,
                           self, None, principal=0.0)
            v = v / swap._notional
            if abs(v) > swapTol:
                print("Swap with maturity " + str(swap._maturity_date)
                      + " Not Repriced. Has Value", v)
                swap.print_fixed_leg_pv()
                swap.print_float_leg_pv()
                raise FinError("Swap not repriced.")

###############################################################################

    # def overnight_rate(self,
    #                   settle_date: Date,
    #                   start_date: Date,
    #                   maturity_date: (Date, list),
    #                   dc_type: DayCountTypes=DayCountTypes.THIRTY_E_360):
    #     """ get a vector of dates and values for the overnight rate implied by
    #     the OIS rate term structure. """

    #     # Note that this function does not call the IborSwap class to
    #     # calculate the swap rate since that will create a circular dependency.
    #     # I therefore recreate the actual calculation of the swap rate here.

    #     if isinstance(maturity_date, Date):
    #         maturity_dates = [maturity_date]
    #     else:
    #         maturity_dates = maturity_date

    #     overnightRates = []

    #     dfValuationDate = self.df(settle_date)

    #     for maturityDt in maturity_dates:

    #         schedule = FinSchedule(start_date,
    #                                maturityDt,
    #                                frequencyType)

    #         flow_dates = schedule._generate()
    #         flow_dates[0] = start_date

    #         day_counter = FinDayCount(dc_type)
    #         prev_dt = flow_dates[0]
    #         pv01 = 0.0
    #         df = 1.0

    #         for next_dt in flow_dates[1:]:
    #             if next_dt > settle_date:
    #                 df = self.df(next_dt) / dfValuationDate
    #                 alpha = day_counter.year_frac(prev_dt, next_dt)[0]
    #                 pv01 += alpha * df

    #             prev_dt = next_dt

    #         if abs(pv01) < gSmall:
    #             par_rate = None
    #         else:
    #             df_start = self.df(start_date)
    #             par_rate = (df_start - df) / pv01

    #         par_rates.append(par_rate)

    #     if isinstance(maturity_date, Date):
    #         return par_rates[0]
    #     else:
    #         return par_rates

###############################################################################

    def __repr__(self):
        """ Print out the details of the Libor curve. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUATION DATE", self._value_date)

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
