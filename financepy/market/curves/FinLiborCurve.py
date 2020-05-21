##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy import optimize

from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...market.curves.FinInterpolate import FinInterpMethods
from ...finutils.FinHelperFunctions import labelToString

swaptol = 1e-8

##############################################################################


def f(df, *args):
    ''' Root search objective function for swaps '''
    curve = args[0]
    valueDate = args[1]
    swap = args[2]
    numPoints = len(curve._times)
    curve._values[numPoints - 1] = df
    v_swap = swap.value(valueDate, curve, curve, None, 1.0)
    v_swap /= swap._notional
    return v_swap

###############################################################################


def g(df, *args):
    ''' Root search objective function for swaps '''
    curve = args[0]
    valueDate = args[1]
    fra = args[2]
    numPoints = len(curve._times)
    curve._values[numPoints - 1] = df
    v_fra = fra.value(valueDate, curve)
    v_fra /= fra._notional
    return v_fra

###############################################################################


class FinLiborCurve(FinDiscountCurve):
    ''' Constructs a discount curve as implied by the prices of Libor
    deposits, FRAs and IRS. The curve date is the date on which we
    are performing the valuation based on the information available on the
    curve date. Typically it is the date on which an amount of $1 paid
    has a present value of $1.

    This class inherits from FinDiscountCurve so has all of the methods
    that class has. '''

    def __init__(self,
                 name,
                 curveDate,
                 liborDeposits,
                 liborFRAs,
                 liborSwaps,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):

        self._name = name
        self._times = []
        self._values = []
        self._curveDate = curveDate
        self._interpMethod = interpMethod
        self.validateInputs(liborDeposits, liborFRAs, liborSwaps)
        self.buildCurve()

    ######################################################################

    def validateInputs(self,
                       liborDeposits,
                       liborFRAs,
                       liborSwaps):
        ''' Construct the discount curve using a bootstrap approach. '''

        numDepos = len(liborDeposits)
        numFRAs = len(liborFRAs)
        numSwaps = len(liborSwaps)

        if numDepos + numFRAs + numSwaps == 0:
            raise FinError("No calibration instruments.")

        # Validation of the inputs.
        if numDepos != 0:
            prevDt = liborDeposits[0]._maturityDate
            for depo in liborDeposits[1:]:
                nextDt = depo._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Deposits must be in increasing maturity")
                prevDt = nextDt
 
        if numFRAs != 0:
            prevDt = liborFRAs[0]._maturityDate
            for fra in liborFRAs[1:]:
                nextDt = fra._maturityDate
                if nextDt <= prevDt:
                    raise FinError("FRAs must be in increasing maturity")
                prevDt = nextDt
 
        if numSwaps != 0:
            prevDt = liborSwaps[0]._maturityDate
            for swap in liborSwaps[1:]:
                nextDt = swap._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Swaps must be in increasing maturity")
                prevDt = nextDt
 
        # Now we have ensure they are in order check for overlaps and the like

        lastDepositMaturityDate = FinDate(1, 1, 1900)
        firstFRAMaturityDate = FinDate(1, 1, 1900)
        lastFRAMaturityDate = FinDate(1, 1, 1900)

        if numDepos > 0:
            lastDepositMaturityDate = liborDeposits[-1]._maturityDate

        if numFRAs > 0:
            firstFRAMaturityDate = liborFRAs[0]._maturityDate
            lastFRAMaturityDate = liborFRAs[-1]._maturityDate

        if numSwaps > 0:
            firstSwapMaturityDate = liborSwaps[0]._maturityDate

        if numDepos > 0 and numFRAs > 0:
            if firstFRAMaturityDate <= lastDepositMaturityDate:
                raise FinError("First FRA must end after last Deposit")

        if numFRAs > 0 and numSwaps > 0:
            if firstSwapMaturityDate <= lastFRAMaturityDate:
                raise FinError("First Swap must mature after last FRA")

        # Now determine which instruments are used
        self._usedDeposits = liborDeposits
        self._usedFRAs = liborFRAs
        self._usedSwaps = liborSwaps

    #######################################################################

    def buildCurve(self):
        ''' Construct the discount curve using a bootstrap approach. '''

        self._times = np.array([])
        self._values = np.array([])

        # time zero is now.
        tmat = 0.0
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._values = np.append(self._values, dfMat)

        for depo in self._usedDeposits:
            tmat = (depo._maturityDate - self._curveDate) / gDaysInYear
            dfMat = depo.maturityDf()
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, dfMat)

        oldtmat = tmat

        for fra in self._usedFRAs:

            tset = (fra._startDate - self._curveDate) / gDaysInYear
            tmat = (fra._maturityDate - self._curveDate) / gDaysInYear

            # if both dates are after the previous FRA/FUT then need to
            # solve for 2 discount factors simultaneously using root search

            if tset < oldtmat and tmat > oldtmat:
                dfMat = fra.maturityDf(self)
                self._times = np.append(self._times, tmat)
                self._values = np.append(self._values, dfMat)
            else:
                self._times = np.append(self._times, tmat)
                self._values = np.append(self._values, dfMat)

                argtuple = (self, self._curveDate, fra)
                dfMat = optimize.newton(g, x0=dfMat, fprime=None,
                                        args=argtuple, tol=swaptol,
                                        maxiter=50, fprime2=None)

        for swap in self._usedSwaps:
            # I use the lastPaymentDate in case a date has been adjusted fwd
            # over a holiday as the maturity date is usually not adjusted CHECK
            maturityDate = swap._lastPaymentDate
            tmat = (maturityDate - self._curveDate) / gDaysInYear

            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, dfMat)

            argtuple = (self, self._curveDate, swap)

            dfMat = optimize.newton(f, x0=dfMat, fprime=None, args=argtuple,
                                    tol=swaptol, maxiter=50, fprime2=None,
                                    full_output=False)

        self.checkRefits()

###############################################################################

    def checkRefits(self):

        for depo in self._usedDeposits:
            v = depo.value(self._curveDate, self) / depo._notional
            if abs(v - 1.0) > swaptol:
                print("Value", v)
                raise FinError("Deposit not repriced.")

        for fra in self._usedFRAs:
            v = fra.value(self._curveDate, self) / fra._notional
            if abs(v) > swaptol:
                print("Value", v)
                raise FinError("FRA not repriced.")

        for swap in self._usedSwaps:
            v = swap.value(self._curveDate, self, self) / swap._notional
            if abs(v) > swaptol*100:
                print("Value", v)
                swap.printFixedLeg(self._curveDate)
                swap.printFloatLeg(self._curveDate)
                raise FinError("Swap not repriced.")

###############################################################################

    def __repr__(self):
        ''' Print out the details of the Libor curve. '''
        numPoints = len(self._times)

        s = labelToString("TIME", "DISCOUNT FACTOR")
        for i in range(0, numPoints):
            s += labelToString(self._times[i], self._values[i])
        return s

###############################################################################
