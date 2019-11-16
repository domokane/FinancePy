# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from scipy import optimize

from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinInterpolate import FinInterpMethods
from ...finutils.FinError import FinError
from ...market.curves.FinDiscountCurve import FinDiscountCurve

#######################################################################


def f(df, *args):
    curve = args[0]
    valueDate = args[1]
    liborSwap = args[2]
    numPoints = len(curve._times)
    curve._values[numPoints - 1] = df
    objFn = liborSwap.fixedLegValue(valueDate, curve, principal=1.0) - 1.0
    return objFn

###############################################################################


class FinLiborOneCurve(FinDiscountCurve):
    ''' Constructs a discount curve as implied by the prices of Libor deposits,
    FRAs and interest rate swaps. The curve date is the date on which we
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

        self.buildCurve(liborDeposits, liborFRAs, liborSwaps)

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

        if numFRAs != 0:
            prevDt = liborFRAs[0]._maturityDate
            for fra in liborFRAs[1:]:
                nextDt = fra._maturityDate
                if nextDt <= prevDt:
                    raise FinError("FRAs must be in increasing maturity")

        if numSwaps != 0:
            prevDt = liborSwaps[0]._maturityDate
            for swap in liborSwaps[1:]:
                nextDt = swap._maturityDate
                if nextDt <= prevDt:
                    raise FinError("Swaps must be in increasing maturity")

        # Now we have ensure they are in order check for overlaps and the like

        if numDepos > 0:
            firstDepositMaturityDate = liborDeposits[0]._maturityDate
            lastDepositMaturityDate = liborDeposits[-1]._maturityDate
        else:
            firstDepositMaturityDate = self._curveDate
            lastDepositMaturityDate = self._curveDate

        if numFRAs > 0:
            firstFRAMaturityDate = liborFRAs[0]._maturityDate
            lastFRAMaturityDate = liborFRAs[-1]._maturityDate
        else:
            firstFRAMaturityDate = self._curveDate
            lastFRAMaturityDate = self._curveDate

        if numSwaps > 0:
            firstSwapMaturityDate = liborSwaps[0]._maturityDate
            lastSwapMaturityDate = liborSwaps[-1]._maturityDate
        else:
            firstSwapMaturityDate = self._curveDate
            lastSwapMaturityDate = self._curveDate

        # Now determine which instruments are used
        self._usedDeposits = liborDeposits
        self._usedFRAs = liborFRAs
        self._usedSwaps = liborSwaps

    #######################################################################

    def buildCurve(self,
                   liborDeposits,
                   liborFRAs,
                   liborSwaps):
        ''' Construct the discount curve using a bootstrap approach. '''

        self.validateInputs(liborDeposits, liborFRAs, liborSwaps)

        self._times = np.array([])
        self._values = np.array([])

        # time zero is now.
        dfMat = 1.0
        self._times = np.append(self._times, 0.0)
        self._values = np.append(self._values, dfMat)

        for depo in self._usedDeposits:
            tmat = (depo._maturityDate - self._curveDate) / gDaysInYear
            dfMat = depo.maturityDf()
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, dfMat)

        for fra in self._usedFRAs:
            tmat = (fra._maturityDate - self._curveDate) / gDaysInYear
            dfMat = fra.maturityDf(self)
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, dfMat)

        for swap in self._usedSwaps:
            maturityDate = swap._maturityDate

            argtuple = (self, self._curveDate, swap)
            tmat = (maturityDate - self._curveDate) / gDaysInYear
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, dfMat)

            optimize.newton(f, x0=dfMat, fprime=None, args=argtuple,
                            tol=1e-8, maxiter=100, fprime2=None)

##########################################################################
