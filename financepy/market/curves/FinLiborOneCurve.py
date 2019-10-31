# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from scipy import optimize

from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDate import FinDate
from ...finutils.FinInterpolate import interpolate
from ...finutils.FinInterpolate import FinInterpMethods
from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...products.libor.FinLiborSwap import FinLiborSwap

#######################################################################

def f(df, *args):
    self = args[0]
    valueDate = args[1]
    liborSwap = args[2]
    numPoints = len(self._times)    
    self._values[numPoints-1] = df
    objFn = liborSwap.fixedLegValue(valueDate, self, principal = 1.0) - 1.0
    return objFn

###############################################################################

class FinLiborOneCurve(FinDiscountCurve):
    ''' Constructs a discount curve as implied by the prices of Libor deposits, 
    futures and interest rate swaps. 
    TODO: Complete the addition of interest rate futures and add OIS '''

    def __init__(self, 
                 name, 
                 curveDate, 
                 liborDeposits,
                 liborFRAs,
                 liborSwaps,
                 interpolationMethod = FinInterpMethods.FLAT_FORWARDS ):

        self._name = name
        self._times = []
        self._values = []
        self._curveDate = curveDate
        self._interpolationMethod = interpolationMethod

        self.buildCurve(
                  liborDeposits,
                  liborFRAs,
                  liborSwaps)

######################################################################

    def buildCurve(self,
                  liborDeposits,
                  liborFRAs,
                  liborSwaps):
        ''' Construct the discount curve using a bootstrap approach. '''

        self._times = np.array([])
        self._values = np.array([])

        self._times = np.append(self._times,0.0)
        self._values = np.append(self._values,1.0)

        df = 1.0

        for depo in liborDeposits:

            tmat = (depo._maturityDate - depo._settlementDate) / gDaysInYear
            dcc = FinDayCount(depo._dayCountType)
            acc = dcc.yearFrac(depo._settlementDate,depo._maturityDate)
            df = 1.0/(1.0 + acc * depo._depositRate)  # this is the interest x accrual

            self._times = np.append(self._times,tmat)
            self._values = np.append(self._values,df)
            
        for fra in liborFRAs:
            pass

        for swap in liborSwaps:

            swapRate = swap._fixedCoupon
            maturityDate = swap._maturityDate
            swapFrequencyType = swap._fixedFrequencyType
            swapBasisType = swap._fixedDayCountType

            liborSwapFixedLeg = FinLiborSwap(self._curveDate,
                                             maturityDate,
                                             swapRate,
                                             swapFrequencyType,
                                             swapBasisType)

            argtuple = (self, self._curveDate, liborSwapFixedLeg)

            tmat = (maturityDate - self._curveDate) / gDaysInYear

            self._times = np.append(self._times,tmat)
            self._values = np.append(self._values,df)

            optimize.newton(f, x0=df, fprime=None, args=argtuple, 
                            tol=1e-8, maxiter=100, fprime2=None)

            df = self._values[-1]

#######################################################################

    def buildCurve2(self,
                  liborDeposits,
                  liborFRAs,
                  liborSwaps):
        ''' Construct the discount curve using a bootstrap approach. '''

        self._times = np.array([])
        self._values = np.array([])

        self._times = np.append(self._times,0.0)
        self._values = np.append(self._values,1.0)

        df = 1.0

        for depo in liborDeposits:

            tmat = (depo._maturityDate - depo._settlementDate) / gDaysInYear
            dcc = FinDayCount(depo._dayCountType)
            acc = dcc.yearFrac(depo._settlementDate,depo._maturityDate)
            df = 1.0/(1.0 + acc * depo._depositRate)  # this is the interest x accrual

            self._times = np.append(self._times,tmat)
            self._values = np.append(self._values,df)
            
        for fra in liborFRAs:
            pass

        # calculate a grid of cashflow dates for all swaps
        # they should lie on a regular grid with overlapping flow dates
        # use linear interpolation of the swap rates to populate all maturity points
        # build curve by iterating over maturity dates
        for swap in liborSwaps:

            swapRate = swap._fixedCoupon
            maturityDate = swap._maturityDate
            swapFrequencyType = swap._fixedFrequencyType
            swapBasisType = swap._fixedDayCountType

            liborSwapFixedLeg = FinLiborSwap(self._curveDate,
                                             maturityDate,
                                             swapRate,
                                             swapFrequencyType,
                                             swapBasisType)

            argtuple = (self, self._curveDate, liborSwapFixedLeg)

            tmat = (maturityDate - self._curveDate) / gDaysInYear

            self._times = np.append(self._times,tmat)
            self._values = np.append(self._values,df)

            optimize.newton(f, x0=df, fprime=None, args=argtuple, 
                            tol=1e-8, maxiter=100, fprime2=None)

            df = self._values[-1]

#######################################################################

    def fwd(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according to the corresponding day count 
        convention. '''

        if date1 < self._curveDate:
            raise FinError("Date " + str(date1) + " before curve value date " + str(self._curveDate))

        if date2 < date1:
            raise FinError("Date2 before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1,date2)

        df1 = self.df(date1)
        df2 = self.df(date2)
        
        fwd = (df1/df2-1.0)/yearFrac

        return fwd

################################################################################

    def df(self, t):
        ''' Determine the discount factor by interpolation. '''

        if type(t) is FinDate:
            t = (t - self._curveDate)/gDaysInYear

        df = interpolate(t, 
                         self._times,
                         self._values,
                         self._interpolationMethod.value)
        return df

################################################################################

    def print(self):
        ''' Print the details '''
        numPoints = len(self._times)
        
        for i in range(0,numPoints):
            print(self._times[i], self._values[i])

################################################################################