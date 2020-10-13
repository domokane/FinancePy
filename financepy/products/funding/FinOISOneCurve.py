##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString, gridIndex
from ...finutils.FinHelperFunctions import checkArgumentTypes, _funcName
from ...finutils.FinGlobalVariables import gDaysInYear
from ...market.curves.FinInterpolate import FinInterpTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

swaptol = 1e-8

##############################################################################
# TODO: CHANGE times to dfTimes
##############################################################################


def _f(df, *args):
    ''' Root search objective function for swaps '''
    curve = args[0]
    valueDate = args[1]
    swap = args[2]
    numPoints = len(curve._times)
    curve._dfValues[numPoints - 1] = df
    v_swap = swap.value(valueDate, curve, curve, None, 1.0)
    v_swap /= swap._notional
    return v_swap

###############################################################################


def _g(df, *args):
    ''' Root search objective function for swaps '''
    curve = args[0]
    valueDate = args[1]
    fra = args[2]
    numPoints = len(curve._times)
    curve._dfValues[numPoints - 1] = df
    v_fra = fra.value(valueDate, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class FinOISSingleCurve(FinDiscountCurve):
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
                 oisFRAs: list,
                 oisSwaps: list,
                 interpType: FinInterpTypes = FinInterpTypes.LINEAR_SWAP_RATES,
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
        self._validateInputs(oisFRAs, oisSwaps)
        self._interpType = interpType
        self._checkRefit = checkRefit
        self._buildCurve()

###############################################################################

    def _buildCurve(self):
        ''' Build curve based on interpolation. '''
        if self._interpType == FinInterpTypes.LINEAR_SWAP_RATES:
            self._buildCurveLinearSwapRateInterpolation()
        else:
            self._buildCurveUsingSolver()

###############################################################################

    def _validateInputs(self,
                        oisFRAs,
                        oisSwaps):
        ''' Validate the inputs for each of the Libor products. '''

        numFRAs = len(oisFRAs)
        numSwaps = len(oisSwaps)

        if numFRAs + numSwaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if numFRAs > 0:
            for fra in oisFRAs:
                startDt = fra._startDate
                if startDt <= self._valuationDate:
                    raise FinError("FRAs starts before valuation date")

        if numFRAs > 1:
            prevDt = oisFRAs[0]._maturityDate
            for fra in oisFRAs[1:]:
                nextDt = fra._maturityDate
                if nextDt <= prevDt:
                    raise FinError("FRAs must be in increasing maturity")
                prevDt = nextDt

        if numSwaps > 0:
            for swap in oisSwaps:
                startDt = swap._startDate
                if startDt < self._valuationDate:
                    raise FinError("Swaps starts before valuation date.")

        if numSwaps > 1:

            # Swaps must all start on the same date for the bootstrap
            startDt = oisSwaps[0]._startDate
            for swap in oisSwaps[1:]:
                nextStartDt = swap._startDate
                if nextStartDt != startDt:
                    raise FinError("Swaps must all have same start date.")

            # Swaps must be increasing in tenor/maturity
            prevDt = oisSwaps[0]._maturityDate
            for swap in oisSwaps[1:]:
                nextDt = swap._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Swaps must be in increasing maturity")
                prevDt = nextDt

            # Swaps must have same cashflows for bootstrap to work
            longestSwap = oisSwaps[-1]
            longestSwapCpnDates = longestSwap._adjustedFixedDates
            for swap in oisSwaps[0:-1]:
                swapCpnDates = swap._adjustedFixedDates
                numFlows = len(swapCpnDates)
                for iFlow in range(0, numFlows):
                    if swapCpnDates[iFlow] != longestSwapCpnDates[iFlow]:
                        raise FinError("Swap coupons are not on the same date grid.")

        #######################################################################
        # Now we have ensure they are in order check for overlaps and the like
        #######################################################################

        lastFRAMaturityDate = FinDate(1, 1, 1900)

        if numFRAs > 0:
            lastFRAMaturityDate = oisFRAs[-1]._maturityDate

        if numSwaps > 0:
            firstSwapMaturityDate = oisSwaps[0]._maturityDate

        if numFRAs > 0 and numSwaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        # Now determine which instruments are used
        self._usedFRAs = oisFRAs
        self._usedSwaps = oisSwaps
        self._dayCountType = None

###############################################################################

    def _buildCurveUsingLinAlg(self):
        ''' Construct the discount curve using a linear algebra approach. This
        is exact and allows spot and forward starting OIS to be fitted. It also
        handles FRA contracts. '''

        numInstruments = len(self._usedFRAs) + len(self._usedSwaps)

        # generate the time grid - start with today
        gridTimes = [0.0]

        for fra in self._usedFRAs:
            tset = (fra._startDate - self._valuationDate) / gDaysInYear
            tmat = (fra._maturityDate - self._valuationDate) / gDaysInYear
            gridTimes.append(tset)
            gridTimes.append(tmat)

        for swap in self._usedSwaps:
            for fixedFlowDate in swap._adjustedFixedDates:
                tflow = (fixedFlowDate - self._valuationDate) / gDaysInYear
                gridTimes.append(tflow)

        gridTimes = list(set(sorted(gridTimes)))
        numTimes = len(gridTimes)
        
        flows = np.zeros(numInstruments, numTimes)
        rhs = np.zeros(numInstruments)

        instrumentCounter = 0        
        for fra in self._usedFRAs:
            tset = (fra._startDate - self._valuationDate) / gDaysInYear
            tmat = (fra._maturityDate - self._valuationDate) / gDaysInYear
            iset = gridIndex(tset, gridTimes)
            imat = gridIndex(tset, gridTimes)
            flows[instrumentCounter, iset] += -1.0
            flows[instrumentCounter, imat] += 1.0 + fra._fraRate
            rhs[instrumentCounter] = 0
            instrumentCounter += 1

        for swap in self._usedSwaps:
            for fixedFlowDate in swap._adjustedFixedDates:
                tflow = (fixedFlowDate - self._valuationDate) / gDaysInYear
                iFlow = gridIndex(tflow, gridTimes)
                accFactor = 1.0
                flows[instrumentCounter][iFlow] = swap._fixedCoupon * accFactor

            tmat = (swap._maturityDt - self._valuationDate) / gDaysInYear
            imat = gridIndex(tmat, gridTimes)
            flows[imat] += 1.0
            rhs[instrumentCounter] = 1.0
            instrumentCounter += 1
            
        print(flows)
        
        dfs = np.linalg.solve(flows, rhs)
        self._times = gridTimes
        self._dfValues = dfs

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
            v = swap.value(swap._startDate, self, self, None, principal=0.0)
            v = v / swap._notional
#            print("REFIT SWAP VALUATION:", swap._adjustedMaturityDate, v)
            if abs(v) > swapTol:
                print("Swap with maturity " + str(swap._maturityDate)
                      + " Not Repriced. Has Value", v)
                swap.printFixedLegPV()
                swap.printFloatLegPV()
                raise FinError("Swap not repriced.")

###############################################################################

    def __repr__(self):
        ''' Print out the details of the Libor curve. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUATION DATE", self._valuationDate)

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
                               "%12.10f" % self._dfValues[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
