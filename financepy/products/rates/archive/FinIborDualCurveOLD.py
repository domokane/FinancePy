##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
import copy

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import check_argument_types, _funcName
from ...finutils.FinGlobalVariables import gDaysInYear
from ...market.curves.FinInterpolator import FinInterpTypes, FinInterpolator
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...products.rates.FinIborDeposit import FinIborDeposit
from ...products.rates.FinIborFRA import FinIborFRA
from ...products.rates.FinIborSwapOLD import FinIborSwapOLD

swaptol = 1e-10

###############################################################################
# TODO: CHANGE times to dfTimes
###############################################################################


def _f(df, *args):
    """ Root search objective function for swaps """
    discount_curve = args[0]
    curve = args[1]
    valuation_date = args[2]
    swap = args[3]

    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df

    # For curves that need a fit function, we fit it now 
    curve._interpolator.fit(curve._times, curve._dfs)     
    v_swap = swap.value(valuation_date, discount_curve, curve, None)

    notional = swap._notional
    v_swap /= notional
    return v_swap

###############################################################################


def _g(df, *args):
    """ Root search objective function for swaps """

    discount_curve = args[0]
    curve = args[1]
    valuation_date = args[2]
    fra = args[3]
    num_points = len(curve._times)
    curve._dfs[num_points - 1] = df

    # For curves that need a fit function, we fit it now 
    curve._interpolator.fit(curve._times, curve._dfs)     
    v_fra = fra.value(valuation_date, discount_curve, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class FinIborDualCurveOLD(FinDiscountCurve):
    """ Constructs an index curve as implied by the prices of Ibor
    deposits, FRAs and IRS. Discounting is assumed to be at a discount rate
    that is an input and usually derived from OIS rates. """

###############################################################################

    def __init__(self,
                 valuation_date: FinDate,
                 discount_curve: FinDiscountCurve,
                 iborDeposits: list,
                 iborFRAs: list,
                 iborSwaps: list,
                 interp_type: FinInterpTypes = FinInterpTypes.FLAT_FWD_RATES,
                 checkRefit: bool = False):  # Set to True to test it works
        """ Create an instance of a FinIbor curve given a valuation date and
        a set of ibor deposits, ibor FRAs and iborSwaps. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on coupon dates and to then assume
        flat forwards between these coupon dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        """

        check_argument_types(getattr(self, _funcName(), None), locals())

        self._valuation_date = valuation_date
        self._discount_curve = discount_curve
        self._validateInputs(iborDeposits, iborFRAs, iborSwaps)
        self._interp_type = interp_type
        self._checkRefit = checkRefit
        self._buildCurve()

###############################################################################

    def _buildCurve(self):
        """ Build curve based on interpolation. """

        self._buildCurveUsing1DSolver()

###############################################################################

    def _validateInputs(self,
                        iborDeposits,
                        iborFRAs,
                        iborSwaps):
        """ Validate the inputs for each of the Ibor products. """

        numDepos = len(iborDeposits)
        numFRAs = len(iborFRAs)
        numSwaps = len(iborSwaps)

        depoStartDate = self._valuation_date
        swapStartDate = self._valuation_date

        if numDepos + numFRAs + numSwaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if numDepos > 0:

            depoStartDate = iborDeposits[0]._start_date

            for depo in iborDeposits:

                if isinstance(depo, FinIborDeposit) is False:
                    raise FinError("Deposit is not of type FinIborDeposit")

                startDt = depo._start_date

                if startDt < self._valuation_date:
                    raise FinError("First deposit starts before value date.")

                if startDt < depoStartDate:
                    depoStartDate = start_date

            for depo in iborDeposits:
                startDt = depo._start_date
                endDt = depo._maturity_date
                if startDt >= endDt:
                    raise FinError("First deposit ends on or before it begins")

        # Ensure order of depos
        if numDepos > 1:

            prevDt = iborDeposits[0]._maturity_date
            for depo in iborDeposits[1:]:
                nextDt = depo._maturity_date
                if nextDt <= prevDt:
                    raise FinError("Deposits must be in increasing maturity")
                prevDt = nextDt

        # REMOVED THIS AS WE WANT TO ANCHOR CURVE AT VALUATION DATE 
        # USE A SYNTHETIC DEPOSIT TO BRIDGE GAP FROM VALUE DATE TO SETTLEMENT DATE
        # Ensure that valuation date is on or after first deposit start date
        # if numDepos > 1:
        #    if iborDeposits[0]._effective_date > self._valuation_date:
        #        raise FinError("Valuation date must not be before first deposit settles.")

        if numFRAs > 0:
            for fra in iborFRAs:
                if isinstance(fra, FinIborFRA) is False:
                    raise FinError("FRA is not of type FinIborFRA")

                startDt = fra._start_date
                if startDt <= self._valuation_date:
                    raise FinError("FRAs starts before valuation date")

        if numFRAs > 1:
            prevDt = iborFRAs[0]._maturity_date
            for fra in iborFRAs[1:]:
                nextDt = fra._maturity_date
                if nextDt <= prevDt:
                    raise FinError("FRAs must be in increasing maturity")
                prevDt = nextDt

        if numSwaps > 0:

            swapStartDate = iborSwaps[0]._effective_date

            for swap in iborSwaps:

                if isinstance(swap, FinIborSwapOLD) is False:
                    raise FinError("Swap is not of type FinIborSwap")

                startDt = swap._effective_date
                if startDt < self._valuation_date:
                    raise FinError("Swaps starts before valuation date.")

                if swap._effective_date < swapStartDate:
                    swapStartDate = swap._effective_date

        if numSwaps > 1:

            # Swaps must all start on the same date for the bootstrap
            startDt = iborSwaps[0]._effective_date
            for swap in iborSwaps[1:]:
                nextStartDt = swap._effective_date
                if nextStartDt != startDt:
                    raise FinError("Swaps must all have same start date.")

            # Swaps must be increasing in tenor/maturity
            prevDt = iborSwaps[0]._maturity_date
            for swap in iborSwaps[1:]:
                nextDt = swap._maturity_date
                if nextDt <= prevDt:
                    raise FinError("Swaps must be in increasing maturity")
                prevDt = nextDt

            # Swaps must have same cash flows for bootstrap to work
            longestSwap = iborSwaps[-1]
            longestSwapCpnDates = longestSwap._adjustedFixedDates
            for swap in iborSwaps[0:-1]:
                swapCpnDates = swap._adjustedFixedDates
                num_flows = len(swapCpnDates)
                for iFlow in range(0, num_flows):
                    if swapCpnDates[iFlow] != longestSwapCpnDates[iFlow]:
                        raise FinError("Swap coupons are not on the same date grid.")

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        lastDepositMaturityDate = FinDate(1, 1, 1900)
        firstFRAMaturityDate = FinDate(1, 1, 1900)
        lastFRAMaturityDate = FinDate(1, 1, 1900)

        if numDepos > 0:
            lastDepositMaturityDate = iborDeposits[-1]._maturity_date

        if numFRAs > 0:
            firstFRAMaturityDate = iborFRAs[0]._maturity_date
            lastFRAMaturityDate = iborFRAs[-1]._maturity_date

        if numSwaps > 0:
            firstSwapMaturityDate = iborSwaps[0]._maturity_date

        if numDepos > 0 and numFRAs > 0:
            if firstFRAMaturityDate <= lastDepositMaturityDate:
                print("FRA Maturity Date:", firstFRAMaturityDate)
                print("Last Deposit Date:", lastDepositMaturityDate)
                raise FinError("First FRA must end after last Deposit")

        if numFRAs > 0 and numSwaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        # If both depos and swaps start after T, we need a rate to get them to
        # the first deposit. So we create a synthetic deposit rate contract.
        
        if swapStartDate > self._valuation_date:

            if numDepos == 0:
                raise FinError("Need a deposit rate to pin down short end.")

            if depoStartDate > self._valuation_date:
                firstDepo = iborDeposits[0]
                if firstDepo._start_date > self._valuation_date:
                    print("Inserting synthetic deposit")
                    syntheticDeposit = copy.deepcopy(firstDepo)
                    syntheticDeposit._effective_date = self._valuation_date
                    syntheticDeposit._maturity_date = firstDepo._start_date
                    iborDeposits.insert(0, syntheticDeposit)
                    numDepos += 1

        # Now determine which instruments are used
        self._usedDeposits = iborDeposits
        self._usedFRAs = iborFRAs
        self._usedSwaps = iborSwaps
        self._day_count_type = None

###############################################################################

    def _buildCurveUsing1DSolver(self):
        """ Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver. """

        self._interpolator = FinInterpolator(self._interp_type)

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, dfMat)
        self._interpolator.fit(self._times, self._dfs)

        # A deposit is not margined and not indexed to Libor so should
        # probably not be used to build an indexed Libor curve from
        for depo in self._usedDeposits:
            dfSettle = self.df(depo._start_date)
            dfMat = depo._maturityDf() * dfSettle
            tmat = (depo._maturity_date - self._valuation_date) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)
            self._interpolator.fit(self._times, self._dfs)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._start_date - self._valuation_date) / gDaysInYear
            tmat = (fra._maturity_date - self._valuation_date) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturityDf(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
                argtuple = (self._discount_curve, self, self._valuation_date, fra)
                dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        for swap in self._usedSwaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturity_date = swap._adjustedFixedDates[-1]
            tmat = (maturity_date - self._valuation_date) / gDaysInYear

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

            argtuple = (self._discount_curve, self, self._valuation_date, swap)

            dfMat = optimize.newton(_f, x0=dfMat, fprime=None, args=argtuple,
                                    tol=swaptol, maxiter=50, fprime2=None,
                                    full_output=False)

        if self._checkRefit is True:
            self._checkRefits(1e-10, swaptol, 1e-5)

###############################################################################

    # def _buildCurveLinearSwapRateInterpolation(self):
    #     """ Construct the discount curve using a bootstrap approach. This is
    #     the linear swap rate method that is fast and exact as it does not
    #     require the use of a solver. It is also market standard. """

    #     self._interpolator = FinInterpolator(self._interp_type)

    #     self._times = np.array([])
    #     self._dfs = np.array([])

    #     # time zero is now.
    #     tmat = 0.0
    #     dfMat = 1.0
    #     self._times = np.append(self._times, 0.0)
    #     self._dfs = np.append(self._dfs, dfMat)
    #     self._interpolator.fit(self._times, self._dfs)

    #     for depo in self._usedDeposits:
    #         dfSettle = self.df(depo._effective_date)
    #         dfMat = depo._maturityDf() * dfSettle
    #         tmat = (depo._maturity_date - self._valuation_date) / gDaysInYear
    #         self._times = np.append(self._times, tmat)
    #         self._dfs = np.append(self._dfs, dfMat)
    #         self._interpolator.fit(self._times, self._dfs)

    #     oldtmat = tmat

    #     for fra in self._usedFRAs:

    #         tset = (fra._start_date - self._valuation_date) / gDaysInYear
    #         tmat = (fra._maturity_date - self._valuation_date) / gDaysInYear

    #         # if both dates are after the previous FRA/FUT then need to
    #         # solve for 2 discount factors simultaneously using root search

    #         if tset < oldtmat and tmat > oldtmat:
    #             dfMat = fra.maturityDf(self)
    #             self._times = np.append(self._times, tmat)
    #             self._dfs = np.append(self._dfs, dfMat)
    #             self._interpolator.fit(self._times, self._dfs)
    #         else:
    #             self._times = np.append(self._times, tmat)
    #             self._dfs = np.append(self._dfs, dfMat)
    #             self._interpolator.fit(self._times, self._dfs)

    #             argtuple = (self._discount_curve, self, self._valuation_date, fra)
    #             dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
    #                                     args=argtuple, tol=swaptol,
    #                                     maxiter=50, fprime2=None)

    #     if len(self._usedSwaps) == 0:
    #         if self._checkRefit is True:
    #             self._checkRefits(1e-10, swaptol, 1e-5)
    #         return

    #     #######################################################################
    #     # ADD SWAPS TO CURVE
    #     #######################################################################

    #     # Find where the FRAs and Depos go up to as this bit of curve is done
    #     foundStart = False
    #     lastDate = self._valuation_date
    #     if len(self._usedDeposits) != 0:
    #         lastDate = self._usedDeposits[-1]._maturity_date

    #     if len(self._usedFRAs) != 0:
    #         lastDate = self._usedFRAs[-1]._maturity_date

    #     # We use the longest swap assuming it has a superset of ALL of the
    #     # swap flow dates used in the curve construction
    #     longestSwap = self._usedSwaps[-1]
    #     couponDates = longestSwap._adjustedFixedDates
    #     num_flows = len(couponDates)

    #     # Find where first coupon without discount factor starts
    #     start_index = 0
    #     for i in range(0, num_flows):
    #         if couponDates[i] > lastDate:
    #             start_index = i
    #             foundStart = True
    #             break

    #     if foundStart is False:
    #         raise FinError("Found start is false. Swaps payments inside FRAs")

    #     swap_rates = []
    #     swapTimes = []

    #     # I use the last coupon date for the swap rate interpolation as this
    #     # may be different from the maturity date due to a holiday adjustment
    #     # and the swap rates need to align with the coupon payment dates
    #     for swap in self._usedSwaps:
    #         swap_rate = swap._fixedCoupon
    #         maturity_date = swap._adjustedFixedDates[-1]
    #         tswap = (maturity_date - self._valuation_date) / gDaysInYear
    #         swapTimes.append(tswap)
    #         swap_rates.append(swap_rate)

    #     interpolatedSwapRates = [0.0]
    #     interpolatedSwapTimes = [0.0]

    #     for dt in couponDates[1:]:
    #         swapTime = (dt - self._valuation_date) / gDaysInYear
    #         swap_rate = np.interp(swapTime, swapTimes, swap_rates)
    #         interpolatedSwapRates.append(swap_rate)
    #         interpolatedSwapTimes.append(swapTime)

    #     # Do I need this line ?
    #     interpolatedSwapRates[0] = interpolatedSwapRates[1]

    #     accrualFactors = longestSwap._fixedYearFracs

    #     acc = 0.0
    #     df = 1.0
    #     pv01 = 0.0
    #     dfSettle = self.df(longestSwap._effective_date)

    #     for i in range(1, start_index):
    #         dt = couponDates[i]
    #         df = self.df(dt)
    #         acc = accrualFactors[i-1]
    #         pv01 += acc * df

    #     for i in range(start_index, num_flows):

    #         dt = couponDates[i]
    #         tmat = (dt - self._valuation_date) / gDaysInYear
    #         swap_rate = interpolatedSwapRates[i]
    #         acc = accrualFactors[i-1]
    #         pv01End = (acc * swap_rate + 1.0)
    #         dfMat = (dfSettle - swap_rate * pv01) / pv01End

    #         self._times = np.append(self._times, tmat)
    #         self._dfs = np.append(self._dfs, dfMat)
    #         self._interpolator.fit(self._times, self._dfs)

    #         pv01 += acc * dfMat

    #     if self._checkRefit is True:
    #         self._checkRefits(1e-10, swaptol, 1e-5)

###############################################################################

    def _checkRefits(self, depoTol, fraTol, swapTol):
        """ Ensure that the Ibor curve refits the calibration instruments. """
        for depo in self._usedDeposits:
            v = depo.value(self._valuation_date, self) / depo._notional
            if abs(v - 1.0) > depoTol:
                print("Value", v)
                raise FinError("Deposit not repriced.")

        for fra in self._usedFRAs:
            v = fra.value(self._valuation_date, 
                          self._discount_curve, self) / fra._notional
            if abs(v) > fraTol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self._usedSwaps:
            # We value it as of the start date of the swap
            v = swap.value(swap._effective_date, self._discount_curve,
                           self, None)
            v = v / swap._notional
            if abs(v) > swapTol:
                print("Swap with maturity " + str(swap._maturity_date)
                      + " Not Repriced. Has Value", v)
                swap.printFixedLegPV()
                swap.printFloatLegPV()
                raise FinError("Swap not repriced.")

###############################################################################

    def __repr__(self):
        """ Print out the details of the Ibor curve. """

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUATION DATE", self._valuation_date)

        for depo in self._usedDeposits:
            s += labelToString("DEPOSIT", "")
            s += depo.__repr__()

        for fra in self._usedFRAs:
            s += labelToString("FRA", "")
            s += fra.__repr__()

        for swap in self._usedSwaps:
            s += labelToString("SWAP", "")
            s += swap.__repr__()

        num_points = len(self._times)

        s += labelToString("INTERP TYPE", self._interp_type)
        s += labelToString("GRID TIMES", "GRID DFS")

        for i in range(0, num_points):
            s += labelToString("% 10.6f" % self._times[i],
                               "%12.10f" % self._dfs[i])
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
