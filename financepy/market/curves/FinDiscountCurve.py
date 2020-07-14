##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...finutils.FinHelperFunctions import inputTime, tableToString
from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from .FinInterpolate import interpolate, FinInterpMethods
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinCalendar import FinCalendarTypes, FinBusDayAdjustTypes
from ...products.libor.FinLiborSwap import FinLiborSwap


###############################################################################
# TODO: Allow it to take in a vector of dates
###############################################################################


def timesFromDates(dt, valuationDate):

    if isinstance(dt, FinDate):
        numDates = 1
        dateList = [dt]
        times = [None]
        times[0] = (dt - valuationDate) / gDaysInYear
    elif isinstance(dt, list) and isinstance(dt[0], FinDate):
        numDates = len(dt)
        dateList = dt
        times = []
        for i in range(0, numDates):
            t = (dateList[i] - valuationDate) / gDaysInYear
            times.append(t)
    else:
        raise FinError("Discount factor must take dates.")

    times = np.array(times)
    return times

###############################################################################


def pv01Times(t, f):
    ''' Calculate a bond style pv01 by calculating remaining coupon times for a
    bond with t years to maturity and a coupon frequency of f. The order of the
    list is reverse time order - it starts with the last coupon date and ends
    with the first coupon date. '''

    dt = 1.0 / f
    pv01Times = []

    while t >= 0.0:
        pv01Times.append(t)
        t -= dt

    return pv01Times

###############################################################################


class FinDiscountCurve():
    ''' This is a base discount curve which has an internal representation of
    a vector of times and discount factors and an interpolation scheme for
    interpolating between these fixed points. '''

###############################################################################

    def __init__(self,
                 valuationDate: FinDate,
                 discountFactorDates,
                 discountFactors,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Create the discount curve from a vector of times and discount
        factors with an anchor date and specify an interpolation scheme. As we
        are explicity linking dates and discount factors, we do not need to
        specify any compounding convention or day count calculation since
        discount factors are pure prices. We do however need to specify a
        convention for interpolating the discount factors in time.'''

        # Validate curve
        if len(discountFactorDates) < 1:
            raise FinError("Times has zero length")

        if len(discountFactorDates) != len(discountFactors):
            raise FinError("Times and Values are not the same")

        self._times = [0.0]
        self._discountFactors = [1.0]

        numPoints = len(discountFactorDates)

        for i in range(1, numPoints):
            t = (discountFactorDates[i] - valuationDate) / gDaysInYear
            self._times.append(t)
            self._discountFactors.append(discountFactors[i])

        self._valuationDate = valuationDate
        self._times = np.array(self._times)

        if testMonotonicity(self._times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._discountFactors = np.array(self._discountFactors)
        self._interpMethod = interpMethod

###############################################################################

    def zeroRate(self,
                 dt: FinDate,
                 frequencyType=FinFrequencyTypes.CONTINUOUS):
        ''' Calculate the zero rate to maturity date. The compounding frequency
        of the rate defaults to continuous which is useful for supplying a rate
        to theoretical models such as Black-Scholes that require a continuously
        compounded zero rate as input. '''

        times = timesFromDates(dt, self._valuationDate)
        zeroRates = self._zeroRate(times, frequencyType)

        if isinstance(dt, FinDate):
            return zeroRates[0]
        else:
            return np.array(zeroRates)

##########################################################################

    def _zeroRate(self,
                  times: list,
                  frequencyType=FinFrequencyTypes.CONTINUOUS):
        ''' Calculate the zero rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use. The compounding frequency defaults to continuous. '''

        f = FinFrequency(frequencyType)
        zerosList = []

        for t in times:

            if t < 0.0:
                raise FinError("Time is negative!")

            if t < 1e-6:
                t = 1e-6

            df = self._df(t)

            if frequencyType == FinFrequencyTypes.SIMPLE:
                zeroRate = (1.0/df-1.0) / t
            elif frequencyType == FinFrequencyTypes.CONTINUOUS:
                zeroRate = -np.log(df) / t
            else:
                zeroRate = (df**(-1.0/t/f) - 1) * f

            zerosList.append(zeroRate)

        return np.array(zerosList)

###############################################################################

    def parRate(self,
                dt: FinDate,
                frequencyType=FinFrequencyTypes.CONTINUOUS):
        ''' Calculate the par rate to maturity date. This is the rate paid by a
        bond that has a price of par today. For a par swap rate, use the swap
        rate function that takes in the swap details. '''

        times = timesFromDates(dt, self._valuationDate)
        parRates = self._parRate(times, frequencyType)

        if isinstance(dt, FinDate):
            return parRates[0]
        else:
            return np.array(parRates)

##########################################################################

    def _parRate(self,
                 times: list,
                 frequencyType=FinFrequencyTypes.ANNUAL):
        ''' Calculate the zero rate to maturity date but with times as inputs.
        This function is used internally and should be discouraged for external
        use. The compounding frequency defaults to continuous. '''

        if frequencyType == FinFrequencyTypes.SIMPLE:
            raise FinError("Cannot calculate par rate with simple yield freq.")
        elif frequencyType == FinFrequencyTypes.CONTINUOUS:
            raise FinError("Cannot calculate par rate with continuous freq.")

        f = FinFrequency(frequencyType)
        dt = 1.0 / f

        if np.any(times < 0.0):
            raise FinError("Unable to calculate par rate as one t < 0.0")

        parRateList = []

        for t in times:

            flowTimes = pv01Times(t, f)
            pv01 = 0.0
            for tFlow in flowTimes:
                pv01 += self._df(tFlow)

            pv01 = pv01 * dt
            dft = self._df(t)
            parRate = (1.0 - dft) / pv01
            parRateList.append(parRate)

        return np.array(parRateList)

##########################################################################

    def swapRate(self,
                 valuationDate,
                 settlementDate,
                 maturityDate,
                 fixedFrequencyType,
                 fixedDayCountType):
        ''' Calculate the breakeven swap rate for an interest rate swap that
        starts on the settlement date (which may be forward starting) with a
        specified frequenct and day count convention. Have omitted calendar and
        other input choices for the moment so default values being used.'''

        coupon = 1.0
        swap = FinLiborSwap(settlementDate,
                            maturityDate,
                            coupon,
                            fixedFrequencyType,
                            fixedDayCountType)

        swapRate = swap.parCoupon(valuationDate, self)
        return swapRate

##########################################################################

    def df(self, dt):
        ''' Function to calculate a discount factor from a date or a
        vector of dates. '''

        times = timesFromDates(dt, self._valuationDate)
        dfs = self._df(times)

        if isinstance(dt, FinDate):
            return dfs[0]
        else:
            return np.array(dfs)

##########################################################################

    def _df(self, t):
        ''' Hidden function to calculate a discount factor from a time or a
        vector of times. Discourage usage in favour of passing in dates. '''
        z = interpolate(t,
                        self._times,
                        self._discountFactors,
                        self._interpMethod.value)
        return z

##########################################################################

    def survProb(self, dt):
        return self.df(dt)

##########################################################################

    def fwd(self, dt):
        ''' Calculate the continuously compounded forward rate at the forward
        FinDate provided. This is done by perturbing the time by a small amount
        and measuring the change in the log of the discount factor divided by
        the time increment dt.'''

        times = timesFromDates(dt, self._valuationDate)
        fwds = self._fwd(times)

        if isinstance(dt, FinDate):
            return fwds[0]
        else:
            return np.array(fwds)

##########################################################################

    def _fwd(self, t):
        ''' Calculate the continuously compounded forward rate at the forward
        time provided. This is done by perturbing the time by a small amount
        and measuring the change in the log of the discount factor divided by
        the time increment dt.'''
        dt = 0.000001

        # Do a double sided calculation if possible
        if np.all(t > dt):
            df1 = self._df(t-dt)
            df2 = self._df(t+dt)
            fwd = np.log(df1/df2)/(2.0*dt)
        else:
            df1 = self._df(t)
            df2 = self._df(t+dt)
            fwd = np.log(df1/df2)/(dt)

        return fwd

##########################################################################

    def bump(self, bumpSize):
        ''' Calculate the continuous forward rate at the forward date. '''

        times = self._times.copy()
        values = self._values.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            values[i] = values[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._curveDate,
                                     times,
                                     values,
                                     self._interpMethod)

        return discCurve

##########################################################################

    def fwdRate(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according to the specified
        day count convention. '''

        if date1 < self._valuationDate:
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

    def __repr__(self):
        header = "TIMES,DISCOUNT FACTORS"
        valueTable = [self._times, self._values]
        precision = "10.7f"

        return tableToString(header, valueTable, precision)

##########################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

#######################################################################
