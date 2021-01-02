##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize
import copy

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinHelperFunctions import checkArgumentTypes, _funcName
from ...finutils.FinGlobalVariables import gDaysInYear
from ...market.curves.FinInterpolator import FinInterpTypes, FinInterpolator
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...products.rates.FinIborDeposit import FinIborDeposit
from ...products.rates.FinOIS import FinOIS

swaptol = 1e-10

##############################################################################
# TODO: CHANGE times to dfTimes
##############################################################################


def _fois(oir, *args):
    ''' Extract the implied overnight index rate assuming it is flat over 
    period in question. '''

    targetOISRate = args[0]
    dayCounter = args[1]
    dateSchedule = args[2]

    startDate = dateSchedule[0]
    endDate = dateSchedule[-1]

    df = 1.0
    prevDt = dateSchedule[0]
    for dt in dateSchedule[1:]:
        yearFrac = dayCounter.yearFrac(prevDt, dt)
        df = df * (1.0 + oir * yearFrac)

    period = dayCounter.yearFrac(startDate, endDate)
    
    OISRate = (df - 1.0) / period
    diff = OISRate - targetOISRate
    return diff

###############################################################################

def _f(df, *args):
    ''' Root search objective function for OIS '''

    curve = args[0]
    valueDate = args[1]
    swap = args[2]
    numPoints = len(curve._times)
    curve._dfs[numPoints - 1] = df

    # For curves that need a fit function, we fit it now 
    curve._interpolator.fit(curve._times, curve._dfs)     
    v_swap = swap.value(valueDate, curve, None)
    notional = swap._fixedLeg._notional
    v_swap /= notional
    return v_swap

###############################################################################


def _g(df, *args):
    ''' Root search objective function for swaps '''
    curve = args[0]
    valueDate = args[1]
    fra = args[2]
    numPoints = len(curve._times)
    curve._dfs[numPoints - 1] = df

    # For curves that need a fit function, we fit it now 
    curve._interpolator.fit(curve._times, curve._dfs)     
    v_fra = fra.value(valueDate, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class FinOISCurve(FinDiscountCurve):
    ''' Constructs a discount curve as implied by the prices of Overnight
    Index Rate swaps. The curve date is the date on which we are
    performing the valuation based on the information available on the
    curve date. Typically it is the date on which an amount of 1 unit paid
    has a present value of 1. This class inherits from FinDiscountCurve
    and so it has all of the methods that that class has.
    
    The construction of the curve is assumed to depend on just the OIS curve, 
    i.e. it does not include information from Ibor-OIS basis swaps. For this
    reason I call it a one-curve.
    '''

###############################################################################

    def __init__(self,
                 valuationDate: FinDate,
                 oisDeposits: list,
                 oisFRAs: list,
                 oisSwaps: list,
                 interpType: FinInterpTypes = FinInterpTypes.FLAT_FWD_RATES,
                 checkRefit: bool = False):  # Set to True to test it works
        ''' Create an instance of an overnight index rate swap curve given a
        valuation date and a set of OIS rates. Some of these may
        be left None and the algorithm will just use what is provided. An
        interpolation method has also to be provided. The default is to use a
        linear interpolation for swap rates on coupon dates and to then assume
        flat forwards between these coupon dates.

        The curve will assign a discount factor of 1.0 to the valuation date.
        '''

        checkArgumentTypes(getattr(self, _funcName(), None), locals())

        self._valuationDate = valuationDate
        self._validateInputs(oisDeposits, oisFRAs, oisSwaps)
        self._interpType = interpType
        self._checkRefit = checkRefit
        self._interpolator = None
        self._buildCurve()

###############################################################################

    def _buildCurve(self):
        ''' Build curve based on interpolation. '''
            
        self._buildCurveUsing1DSolver()

###############################################################################

    def _validateInputs(self,
                        oisDeposits,
                        oisFRAs,
                        oisSwaps):
        ''' Validate the inputs for each of the Libor products. '''

        numDepos = len(oisDeposits)
        numFRAs = len(oisFRAs)
        numSwaps = len(oisSwaps)

        depoStartDate = self._valuationDate
        swapStartDate = self._valuationDate

        if numDepos + numFRAs + numSwaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if numDepos > 0:

            depoStartDate = oisDeposits[0]._startDate

            for depo in oisDeposits:

                if isinstance(depo, FinIborDeposit) is False:
                    raise FinError("Deposit is not of type FinIborDeposit")

                startDate = depo._startDate

                if startDate < self._valuationDate:
                    raise FinError("First deposit starts before value date.")

                if startDate < depoStartDate:
                    depoStartDate = startDate

            for depo in oisDeposits:
                startDt = depo._startDate
                endDt = depo._maturityDate
                if startDt >= endDt:
                    raise FinError("First deposit ends on or before it begins")

        # Ensure order of depos
        if numDepos > 1:

            prevDt = oisDeposits[0]._maturityDate
            for depo in oisDeposits[1:]:
                nextDt = depo._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Deposits must be in increasing maturity")
                prevDt = nextDt

        # REMOVED THIS AS WE WANT TO ANCHOR CURVE AT VALUATION DATE 
        # USE A SYNTHETIC DEPOSIT TO BRIDGE GAP FROM VALUE DATE TO SETTLEMENT DATE
        # Ensure that valuation date is on or after first deposit start date
        # Ensure that valuation date is on or after first deposit start date
        # if numDepos > 1:
        #    if oisDeposits[0]._startDate > self._valuationDate:
        #        raise FinError("Valuation date must not be before first deposit settles.")

        # Validation of the inputs.
        if numFRAs > 0:
            for fra in oisFRAs:
                startDt = fra._startDate
                if startDt <= self._valuationDate:
                    raise FinError("FRAs starts before valuation date")

                startDt = fra._startDate
                if startDt < self._valuationDate:
                    raise FinError("FRAs starts before valuation date")

        if numFRAs > 1:
            prevDt = oisFRAs[0]._maturityDate
            for fra in oisFRAs[1:]:
                nextDt = fra._maturityDate
                if nextDt <= prevDt:
                    raise FinError("FRAs must be in increasing maturity")
                prevDt = nextDt

        if numSwaps > 0:

            swapStartDate = oisSwaps[0]._effectiveDate

            for swap in oisSwaps:

                if isinstance(swap, FinOIS) is False:
                    raise FinError("Swap is not of type FinOIS")

                startDt = swap._effectiveDate
                if startDt < self._valuationDate:
                    raise FinError("Swaps starts before valuation date.")

                if swap._effectiveDate < swapStartDate:
                    swapStartDate = swap._effectiveDate

        if numSwaps > 1:

            # Swaps must be increasing in tenor/maturity
            prevDt = oisSwaps[0]._maturityDate
            for swap in oisSwaps[1:]:
                nextDt = swap._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Swaps must be in increasing maturity")
                prevDt = nextDt

        # TODO: REINSTATE THESE CHECKS ?
            # Swaps must have same cashflows for linear swap bootstrap to work
#            longestSwap = oisSwaps[-1]
#            longestSwapCpnDates = longestSwap._adjustedFixedDates
#            for swap in oisSwaps[0:-1]:
#                swapCpnDates = swap._adjustedFixedDates
#                numFlows = len(swapCpnDates)
#                for iFlow in range(0, numFlows):
#                    if swapCpnDates[iFlow] != longestSwapCpnDates[iFlow]:
#                        raise FinError("Swap coupons are not on the same date grid.")

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        lastDepositMaturityDate = FinDate(1, 1, 1900)
        firstFRAMaturityDate = FinDate(1, 1, 1900)
        lastFRAMaturityDate = FinDate(1, 1, 1900)

        if numDepos > 0:
            lastDepositMaturityDate = oisDeposits[-1]._maturityDate

        if numFRAs > 0:
            firstFRAMaturityDate = oisFRAs[0]._maturityDate
            lastFRAMaturityDate = oisFRAs[-1]._maturityDate

        if numSwaps > 0:
            firstSwapMaturityDate = oisSwaps[0]._maturityDate

        if numFRAs > 0 and numSwaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA ends")

        if numDepos > 0 and numFRAs > 0:
            if firstFRAMaturityDate <= lastDepositMaturityDate:
                print("FRA Maturity Date:", firstFRAMaturityDate)
                print("Last Deposit Date:", lastDepositMaturityDate)
                raise FinError("First FRA must end after last Deposit")

        if numFRAs > 0 and numSwaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        if swapStartDate > self._valuationDate:

            if numDepos == 0:
                raise FinError("Need a deposit rate to pin down short end.")

            if depoStartDate > self._valuationDate:
                firstDepo = oisDeposits[0]
                if firstDepo._startDate > self._valuationDate:
                    print("Inserting synthetic deposit")
                    syntheticDeposit = copy.deepcopy(firstDepo)
                    syntheticDeposit._startDate = self._valuationDate
                    syntheticDeposit._maturityDate = firstDepo._startDate
                    oisDeposits.insert(0, syntheticDeposit)
                    numDepos += 1

        # Now determine which instruments are used
        self._usedDeposits = oisDeposits
        self._usedFRAs = oisFRAs
        self._usedSwaps = oisSwaps
        self._dayCountType = None

###############################################################################

    def _buildCurveUsing1DSolver(self):
        ''' Construct the discount curve using a bootstrap approach. This is
        the non-linear slower method that allows the user to choose a number
        of interpolation approaches between the swap rates and other rates. It
        involves the use of a solver. '''

        self._interpolator = FinInterpolator(self._interpType)
        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, dfMat)
        self._interpolator.fit(self._times, self._dfs)

        for depo in self._usedDeposits:
            dfSettle = self.df(depo._startDate)
            dfMat = depo._maturityDf() * dfSettle
            tmat = (depo._maturityDate - self._valuationDate) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)
            self._interpolator.fit(self._times, self._dfs)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._startDate - self._valuationDate) / gDaysInYear
            tmat = (fra._maturityDate - self._valuationDate) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturityDf(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
                argtuple = (self, self._valuationDate, fra)
                dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        for swap in self._usedSwaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturityDate = swap._fixedLeg._paymentDates[-1]
            tmat = (maturityDate - self._valuationDate) / gDaysInYear

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)

            argtuple = (self, self._valuationDate, swap)

            dfMat = optimize.newton(_f, x0=dfMat, fprime=None, args=argtuple,
                                    tol=swaptol, maxiter=50, fprime2=None,
                                    full_output=False)

        if self._checkRefit is True:
            self._checkRefits(1e-10, swaptol, 1e-5)

###############################################################################

    def _buildCurveLinearSwapRateInterpolation(self):
        ''' Construct the discount curve using a bootstrap approach. This is
        the linear swap rate method that is fast and exact as it does not
        require the use of a solver. It is also market standard. '''

        self._interpolator = FinInterpolator(self._interpType)

        self._times = np.array([])
        self._dfs = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._dfs = np.append(self._dfs, dfMat)

        for depo in self._usedDeposits:
            dfSettle = self.df(depo._startDate)
            dfMat = depo._maturityDf() * dfSettle
            tmat = (depo._maturityDate - self._valuationDate) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)
            self._interpolator.fit(self._times, self._dfs)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._startDate - self._valuationDate) / gDaysInYear
            tmat = (fra._maturityDate - self._valuationDate) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturityDf(self)
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
                self._interpolator.fit(self._times, self._dfs)
            else:
                self._times = np.append(self._times, tmat)
                self._dfs = np.append(self._dfs, dfMat)
                self._interpolator.fit(self._times, self._dfs)

                argtuple = (self, self._valuationDate, fra)
                dfMat = optimize.newton(_g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        if len(self._usedSwaps) == 0:
            if self._checkRefit is True:
                self._checkRefits(1e-10, swaptol, 1e-5)
            return

        #######################################################################
        # ADD SWAPS TO CURVE
        #######################################################################

        # Find where the FRAs and Depos go up to as this bit of curve is done
        foundStart = False
        lastDate = self._valuationDate
        if len(self._usedDeposits) != 0:
            lastDate = self._usedDeposits[-1]._maturityDate

        if len(self._usedFRAs) != 0:
            lastDate = self._usedFRAs[-1]._maturityDate

        # We use the longest swap assuming it has a superset of ALL of the
        # swap flow dates used in the curve construction
        longestSwap = self._usedSwaps[-1]
        couponDates = longestSwap._adjustedFixedDates
        numFlows = len(couponDates)

        # Find where first coupon without discount factor starts
        startIndex = 0
        for i in range(0, numFlows):
            if couponDates[i] > lastDate:
                startIndex = i
                foundStart = True
                break

        if foundStart is False:
            raise FinError("Found start is false. Swaps payments inside FRAs")

        swapRates = []
        swapTimes = []

        # I use the last coupon date for the swap rate interpolation as this
        # may be different from the maturity date due to a holiday adjustment
        # and the swap rates need to align with the coupon payment dates
        for swap in self._usedSwaps:
            swapRate = swap._fixedCoupon
            maturityDate = swap._adjustedFixedDates[-1]
            tswap = (maturityDate - self._valuationDate) / gDaysInYear
            swapTimes.append(tswap)
            swapRates.append(swapRate)

        interpolatedSwapRates = [0.0]
        interpolatedSwapTimes = [0.0]

        for dt in couponDates[1:]:
            swapTime = (dt - self._valuationDate) / gDaysInYear
            swapRate = np.interp(swapTime, swapTimes, swapRates)
            interpolatedSwapRates.append(swapRate)
            interpolatedSwapTimes.append(swapTime)

        # Do I need this line ?
        interpolatedSwapRates[0] = interpolatedSwapRates[1]

        accrualFactors = longestSwap._fixedYearFracs

        acc = 0.0
        df = 1.0
        pv01 = 0.0
        dfSettle = self.df(longestSwap._startDate)

        for i in range(1, startIndex):
            dt = couponDates[i]
            df = self.df(dt)
            acc = accrualFactors[i-1]
            pv01 += acc * df

        for i in range(startIndex, numFlows):

            dt = couponDates[i]
            tmat = (dt - self._valuationDate) / gDaysInYear
            swapRate = interpolatedSwapRates[i]
            acc = accrualFactors[i-1]
            pv01End = (acc * swapRate + 1.0)

            dfMat = (dfSettle - swapRate * pv01) / pv01End

            self._times = np.append(self._times, tmat)
            self._dfs = np.append(self._dfs, dfMat)
            self._interpolator.fit(self._times, self._dfs)

            pv01 += acc * dfMat

        if self._checkRefit is True:
            self._checkRefits(1e-10, swaptol, 1e-5)

###############################################################################

    def _checkRefits(self, depoTol, fraTol, swapTol):
        ''' Ensure that the Libor curve refits the calibration instruments. '''

        for fra in self._usedFRAs:
            v = fra.value(self._valuationDate, self) / fra._notional
            if abs(v) > fraTol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self._usedSwaps:
            # We value it as of the start date of the swap
            v = swap.value(swap._effectiveDate, self, self, None, principal=0.0)
            v = v / swap._notional
            if abs(v) > swapTol:
                print("Swap with maturity " + str(swap._maturityDate)
                      + " Not Repriced. Has Value", v)
                swap.printFixedLegPV()
                swap.printFloatLegPV()
                raise FinError("Swap not repriced.")

###############################################################################

    # def overnightRate(self,
    #                   settlementDate: FinDate,
    #                   startDate: FinDate,
    #                   maturityDate: (FinDate, list),
    #                   dayCountType: FinDayCountTypes=FinDayCountTypes.THIRTY_E_360):
    #     ''' get a vector of dates and values for the overnight rate implied by
    #     the OIS rate term structure. '''

    #     # Note that this function does not call the FinIborSwap class to
    #     # calculate the swap rate since that will create a circular dependency.
    #     # I therefore recreate the actual calculation of the swap rate here.

    #     if isinstance(maturityDate, FinDate):
    #         maturityDates = [maturityDate]
    #     else:
    #         maturityDates = maturityDate

    #     overnightRates = []

    #     dfValuationDate = self.df(settlementDate)

    #     for maturityDt in maturityDates:

    #         schedule = FinSchedule(startDate,
    #                                maturityDt,
    #                                frequencyType)

    #         flowDates = schedule._generate()
    #         flowDates[0] = startDate

    #         dayCounter = FinDayCount(dayCountType)
    #         prevDt = flowDates[0]
    #         pv01 = 0.0
    #         df = 1.0

    #         for nextDt in flowDates[1:]:                
    #             if nextDt > settlementDate:
    #                 df = self.df(nextDt) / dfValuationDate
    #                 alpha = dayCounter.yearFrac(prevDt, nextDt)[0]
    #                 pv01 += alpha * df

    #             prevDt = nextDt

    #         if abs(pv01) < gSmall:
    #             parRate = None
    #         else:
    #             dfStart = self.df(startDate)
    #             parRate = (dfStart - df) / pv01

    #         parRates.append(parRate)

    #     if isinstance(maturityDate, FinDate):
    #         return parRates[0]
    #     else:
    #         return parRates

###############################################################################

    def __repr__(self):
        ''' Print out the details of the Libor curve. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUATION DATE", self._valuationDate)

        for depo in self._usedDeposits:
            s += labelToString("DEPOSIT", "")
            s += depo.__repr__()

        for fra in self._usedFRAs:
            s += labelToString("FRA", "")
            s += fra.__repr__()

        for swap in self._usedSwaps:
            s += labelToString("SWAP", "")
            s += swap.__repr__()

        numPoints = len(self._times)

        s += labelToString("INTERP TYPE", self._interpType)

        s += labelToString("GRID TIMES", "GRID DFS")
        for i in range(0, numPoints):
            s += labelToString("% 10.6f" % self._times[i],
                               "%12.10f" % self._dfs[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
