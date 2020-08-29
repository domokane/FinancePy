##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from .FinInterpolate import interpolate, FinInterpTypes
from ...products.libor.FinLiborSwap import FinLiborSwap
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinHelperFunctions import timesFromDates
from ...finutils.FinHelperFunctions import pv01Times
from ...finutils.FinHelperFunctions import labelToString

###############################################################################


class FinDiscountCurve():
    ''' This is a base discount curve which has an internal representation of
    a vector of times and discount factors and an interpolation scheme for
    interpolating between these fixed points. '''

###############################################################################

    def __init__(self,
                 valuationDate: FinDate,
                 dfDates: list,
                 dfValues: np.ndarray,
                 interpType: FinInterpTypes = FinInterpTypes.FLAT_FORWARDS):
        ''' Create the discount curve from a vector of times and discount
        factors with an anchor date and specify an interpolation scheme. As we
        are explicity linking dates and discount factors, we do not need to
        specify any compounding convention or day count calculation since
        discount factors are pure prices. We do however need to specify a
        convention for interpolating the discount factors in time.'''

        checkArgumentTypes(self.__init__, locals())

        # Validate curve
        if len(dfDates) < 1:
            raise FinError("Times has zero length")

        if len(dfDates) != len(dfValues):
            raise FinError("Times and Values are not the same")

        self._times = [0.0]
        self._dfValues = [1.0]
        self._dfDates = dfDates

        numPoints = len(dfDates)

        startIndex = 0
        if numPoints > 0:
            if dfDates[0] == valuationDate:
                self._dfValues[0] = dfValues[0]
                startIndex = 1

        for i in range(startIndex, numPoints):
            t = (dfDates[i] - valuationDate) / gDaysInYear
            self._times.append(t)
            self._dfValues.append(dfValues[i])

        self._times = np.array(self._times)

        if testMonotonicity(self._times) is False:
            print(self._times)
            raise FinError("Times are not sorted in increasing order")

        self._valuationDate = valuationDate
        self._dfValues = np.array(self._dfValues)
        self._interpType = interpType
        self._frequencyType = FinFrequencyTypes.CONTINUOUS

###############################################################################

    def _dfToZero(self,
                  maturityDt: (FinDate, list),
                  frequencyType: FinFrequencyTypes,
                  dayCountType: FinDayCountTypes):
        ''' Convert a discount factor to a zero rate with specific compounding
        frequency which may be continuous, simple, or compounded at a specific
        frequency which are all choices of FinFrequencyTypes. Returns a list of
        discount factor. '''

        f = FinFrequency(frequencyType)
        dcCounter = FinDayCount(dayCountType)

        if isinstance(maturityDt, FinDate):
            dateList = [maturityDt]
        else:
            dateList = maturityDt

        zeroRates = []
        for dt in dateList:
            df = self.df(dt)

            t = dcCounter.yearFrac(self._valuationDate, dt)
            t = max(t, gSmall)

            if frequencyType == FinFrequencyTypes.CONTINUOUS:
                r = -np.log(df)/t
            elif frequencyType == FinFrequencyTypes.SIMPLE:
                r = (1.0/df - 1.0)/t
            else:
                r = (np.power(df, -1.0/(t * f))-1.0) * f
            zeroRates.append(r)

        return zeroRates

###############################################################################

    def zeroRate(self,
                 dt: (list, FinDate),
                 frequencyType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360):
        ''' Calculation of zero rates with specified frequency. This
        function can return a vector of zero rates given a vector of
        times so must use Numpy functions. Default frequency is a
        continuously compounded rate. '''

        if isinstance(frequencyType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if isinstance(dayCountType, FinDayCountTypes) is False:
            raise FinError("Invalid Day Count type.")

        # Calculate the zero rate taking into account frequency AND day count
        zeroRates = self._dfToZero(dt, frequencyType, dayCountType)

        if isinstance(dt, FinDate):
            return zeroRates[0]
        else:
            return np.array(zeroRates)

###############################################################################

    # def _zeroRate(self,
    #               times: (float, np.ndarray),
    #               frequencyType=FinFrequencyTypes):
    #     ''' Calculate the zero rate to maturity date but with times as inputs.
    #     This function is used internally and should be discouraged for external
    #     use. The compounding frequency defaults to continuous. '''

    #     if isinstance(frequencyType, FinFrequencyTypes) is False:
    #         raise FinError("Invalid Frequency type.")

    #     if isinstance(times, float):
    #         times = np.array([times])

    #     if np.any(times < 0.0):
    #         raise FinError("All times must be positive")

    #     times = np.maximum(times, 1e-6)
    #     dfs = self._df(times)
    #     zeroRates = dfToZero(dfs, times, frequencyType)
    #     return zeroRates

###############################################################################

    def parRate(self,
                dt: FinDate,
                frequencyType=FinFrequencyTypes.ANNUAL):
        ''' Calculate the par rate to maturity date. This is the rate paid by a
        bond that has a price of par today. For a par swap rate, use the swap
        rate function that takes in the swap details. '''

        if isinstance(frequencyType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

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

        if isinstance(frequencyType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

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
                 valuationDate: FinDate,
                 settlementDate: FinDate,
                 maturityDate: FinDate,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes):
        ''' Calculate the breakeven swap rate for an interest rate swap that
        starts on the settlement date (which may be forward starting) with a
        specified frequency and day count convention. Have omitted calendar and
        other input choices for the moment so default values being used.'''

        if isinstance(fixedFrequencyType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        coupon = 1.0
        swap = FinLiborSwap(settlementDate,
                            maturityDate,
                            coupon,
                            fixedFrequencyType,
                            fixedDayCountType)

        swapRate = swap.parCoupon(valuationDate, self)
        return swapRate

##########################################################################

    def df(self,
           dt: (list, FinDate)):
        ''' Function to calculate a discount factor from a date or a
        vector of dates. '''

        times = timesFromDates(dt, self._valuationDate)
        dfs = self._df(times)

        if isinstance(dt, FinDate):
            return dfs[0]
        else:
            return np.array(dfs)

##########################################################################

    def _df(self,
            t: (float, np.ndarray)):
        ''' Hidden function to calculate a discount factor from a time or a
        vector of times. Discourage usage in favour of passing in dates. '''
        z = interpolate(t,
                        self._times,
                        self._dfValues,
                        self._interpType.value)
        return z

##########################################################################

    def survProb(self,
                 dt: FinDate):
        ''' This returns a survival probability to a specified date based on
        the assumption that the continuously compounded rate is a default
        hazard rate in which case the survival probability is directly
        analagous to a discount factor. '''

        return self.df(dt)

##########################################################################

    def fwd(self,
            dt: FinDate):
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

    def _fwd(self,
             times: (np.ndarray, float)):
        ''' Calculate the continuously compounded forward rate at the forward
        time provided. This is done by perturbing the time by a small amount
        and measuring the change in the log of the discount factor divided by
        the time increment dt.'''

        dt = 1e-6
        times = np.maximum(times, dt)

        df1 = self._df(times-dt)
        df2 = self._df(times+dt)
        fwd = np.log(df1/df2)/(2.0*dt)
        return fwd

##########################################################################

    def bump(self,
             bumpSize: float):
        ''' Adjust the continuously compounded forward rates by a perturbation
        upward equal to the bump size and return a curve objet with this bumped
        curve. This is used for interest rate risk. '''

        times = self._times.copy()
        values = self._discountFactors.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            values[i] = values[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._valuationDate,
                                     times,
                                     values,
                                     self._interpType)

        return discCurve

##########################################################################

    def fwdRate(self,
                startDate: (list, FinDate),
                endDate: (list, FinDate),
                dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360):
        ''' Calculate the forward rate between two forward dates according to
        the specified day count convention. This defaults to Actual 360. '''

        if isinstance(startDate, FinDate) and isinstance(endDate, FinDate):
            startDates = []
            startDates.append(startDate)
            endDates = []
            endDates.append(endDate)
        elif isinstance(startDate, list) and isinstance(endDate, list):
            if len(startDate) != len(endDate):
                raise FinError("Start, end date lists must have same length")

            startDates = startDate
            endDates = endDate
        else:
            raise FinError("Start date and end date must be same types.")

        dayCount = FinDayCount(dayCountType)

        numDates = len(startDates)
        fwdRates = []
        for i in range(0, numDates):
            dt1 = startDates[i]
            dt2 = endDates[i]
            yearFrac = dayCount.yearFrac(dt1, dt2)
            df1 = self.df(dt1)
            df2 = self.df(dt2)
            fwdRate = (df1 / df2 - 1.0) / yearFrac
            fwdRates.append(fwdRate)

        if isinstance(startDate, FinDate):
            return fwdRates[0]
        else:
            return fwdRates

###############################################################################

    def __repr__(self):
        s = type(self).__name__ + "\n"
        numPoints = len(self._dfDates)
        s += labelToString("DATES", "DISCOUNT FACTORS")
        for i in range(0, numPoints):
            s += labelToString(self._dfDates[i], self._dfValues[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
