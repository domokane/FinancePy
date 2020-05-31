##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################



import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from ...finutils.FinMath import scale, testMonotonicity
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinHelperFunctions import inputTime, inputFrequency
from ...market.curves.FinInterpolate import FinInterpMethods, interpolate
from ...finutils.FinHelperFunctions import labelToString

##############################################################################


def f(df, *args):
    curve = args[0]
    valueDate = args[1]
    bond = args[2]
    marketCleanPrice = args[3]
    numPoints = len(curve._times)
    curve._values[numPoints - 1] = df
    bondDiscountPrice = bond.cleanValueFromDiscountCurve(valueDate, curve)
    objFn = bondDiscountPrice - marketCleanPrice
    return objFn

##############################################################################


class FinBondZeroCurve():
    ''' Class to do bootstrap exact fitting of the bond zero rate curve. '''

    def __init__(self, settlementDate, bonds, cleanPrices,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Fit a discount curve to a set of bond yields using the type of 
        curve specified. '''

        if len(bonds) != len(cleanPrices):
            raise ValueError("Num bonds does not equal number of prices.")

        self._settlementDate = settlementDate
        self._curveDate = settlementDate
        self._bonds = bonds
        self._cleanPrices = np.array(cleanPrices)
        self._discountCurve = None
        self._interpMethod = interpMethod

        times = []
        for bond in self._bonds:
            tmat = (bond._maturityDate-self._settlementDate)/gDaysInYear
            times.append(tmat)

        times = np.array(times)
        if testMonotonicity(times) is False:
            raise ValueError("Times are not sorted in increasing order")

        self._yearsToMaturity = np.array(times)

        self.bootstrapZeroRates()

##############################################################################

    def bootstrapZeroRates(self):

        self._times = np.array([0.0])
        self._values = np.array([1.0])
        df = 1.0

        for i in range(0, len(self._bonds)):
            bond = self._bonds[i]
            maturityDate = bond._maturityDate
            cleanPrice = self._cleanPrices[i]
            tmat = (maturityDate - self._settlementDate) / gDaysInYear
            argtuple = (self, self._settlementDate, bond, cleanPrice)
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, df)

            optimize.newton(f, x0=df, fprime=None, args=argtuple,
                            tol=1e-8, maxiter=100, fprime2=None)

##########################################################################

    def zeroRate(self, dt, compoundingFreq=-1):
        ''' Calculate the zero rate to maturity date. '''
        t = inputTime(dt, self)
        f = inputFrequency(compoundingFreq)
        df = self.df(t)

        if f == 0:  # Simple interest
            zeroRate = (1.0/df-1.0)/t
        if f == -1:  # Continuous
            zeroRate = -np.log(df) / t
        else:
            zeroRate = (df**(-1.0/t) - 1) * f
        return zeroRate

##########################################################################

    def df(self, dt):
        t = inputTime(dt, self)
        z = interpolate(t, self._times, self._values,
                        self._interpMethod.value)
        return z

##########################################################################

    def survProb(self, dt):
        t = inputTime(dt, self)
        q = interpolate(t, self._times, self._values,
                        self._interpMethod.value)
        return q

##########################################################################

    def fwd(self, dt):
        ''' Calculate the continuous forward rate at the forward date. '''
        t = inputTime(dt, self)
        dt = 0.000001
        df1 = self.df(t)
        df2 = self.df(t+dt)
        fwd = np.log(df1/df2)/dt
        return fwd

##########################################################################

    def fwdRate(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according to the specified
        day count convention. '''

        if date1 < self._curveDate:
            raise ValueError("Date1 before curve value date.")

        if date2 < date1:
            raise ValueError("Date2 must not be before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1, date2)
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / yearFrac
        return fwd

##########################################################################

    def plot(self, title):
        ''' Display yield curve. '''

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel('Time to Maturity (years)')
        plt.ylabel('Zero Rate (%)')

        tmax = np.max(self._yearsToMaturity)
        t = np.linspace(0.0, int(tmax+0.5), 100)

        zeroRate = self.zeroRate(t)
        zeroRate = scale(zeroRate, 100.0)
        plt.plot(t, zeroRate, label="Zero Rate Bootstrap", marker='o')
        plt.legend(loc='lower right')
        plt.ylim((min(zeroRate)-0.3, max(zeroRate)*1.1))
        plt.grid(True)

##########################################################################

    def __repr__(self):
        numPoints = len(self._times)
        s = labelToString("TIMES,DISCOUNT FACTORS")
        for i in range(0, numPoints):
            s += labelToString(self._times[i], self._values[i])

        return s

##########################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

##########################################################################
