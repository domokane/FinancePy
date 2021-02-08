##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from .FinInterpolator import FinInterpolator, FinInterpTypes, interpolate

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinHelperFunctions import timesFromDates
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
                 interpType: FinInterpTypes = FinInterpTypes.FLAT_FWD_RATES):
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
        self._dfs = [1.0]
        self._dfDates = dfDates

        numPoints = len(dfDates)

        startIndex = 0
        if numPoints > 0:
            if dfDates[0] == valuationDate:
                self._dfs[0] = dfValues[0]
                startIndex = 1

        for i in range(startIndex, numPoints):
            t = (dfDates[i] - valuationDate) / gDaysInYear
            self._times.append(t)
            self._dfs.append(dfValues[i])

        self._times = np.array(self._times)

        if testMonotonicity(self._times) is False:
            print(self._times)
            raise FinError("Times are not sorted in increasing order")

        self._valuationDate = valuationDate
        self._dfs = np.array(self._dfs)
        self._interpType = interpType
        self._freqType = FinFrequencyTypes.CONTINUOUS
        self._dayCountType = None  # Not needed for this curve
        self._interpolator = FinInterpolator(self._interpType)
        self._interpolator.fit(self._times, self._dfs)

###############################################################################

    def _zeroToDf(self,
                  valuationDate: FinDate,
                  rates: (float, np.ndarray),
                  times: (float, np.ndarray),
                  freqType: FinFrequencyTypes,
                  dayCountType: FinDayCountTypes):
        ''' Convert a zero with a specified compounding frequency and day count
        convention to a discount factor for a single maturity date or a list of
        dates. The day count is used to calculate the elapsed year fraction.'''

        if isinstance(times, float):
            times = np.array([times])

        t = np.maximum(times, gSmall)

        f = FinFrequency(freqType)

        if freqType == FinFrequencyTypes.CONTINUOUS:
            df = np.exp(-rates*t)
        elif freqType == FinFrequencyTypes.SIMPLE:
            df = 1.0 / (1.0 + rates * t)
        elif freqType == FinFrequencyTypes.ANNUAL or \
            freqType == FinFrequencyTypes.SEMI_ANNUAL or \
                freqType == FinFrequencyTypes.QUARTERLY or \
                freqType == FinFrequencyTypes.MONTHLY:
            df = 1.0 / np.power(1.0 + rates/f, f * t)
        else:
            raise FinError("Unknown Frequency type")

        return df

###############################################################################

    def _dfToZero(self,
                  dfs: (float, np.ndarray),
                  maturityDts: (FinDate, list),
                  freqType: FinFrequencyTypes,
                  dayCountType: FinDayCountTypes):
        ''' Given a dates this first generates the discount factors. It then
        converts the discount factors to a zero rate with a chosen compounding
        frequency which may be continuous, simple, or compounded at a specific
        frequency which are all choices of FinFrequencyTypes. Returns a list of
        discount factor. '''

        f = FinFrequency(freqType)

        if isinstance(maturityDts, FinDate):
            dateList = [maturityDts]
        else:
            dateList = maturityDts

        if isinstance(dfs, float):
            dfList = [dfs]
        else:
            dfList = dfs

        if len(dateList) != len(dfList):
            raise FinError("Date list and df list do not have same length")

        numDates = len(dateList)
        zeroRates = []

        times = timesFromDates(dateList, self._valuationDate, dayCountType)

        for i in range(0, numDates):

            df = dfList[i]

            t = max(times[i], gSmall)

            if freqType == FinFrequencyTypes.CONTINUOUS:
                r = -np.log(df)/t
            elif freqType == FinFrequencyTypes.SIMPLE:
                r = (1.0/df - 1.0)/t
            else:
                r = (np.power(df, -1.0/(t * f))-1.0) * f

            zeroRates.append(r)

        return np.array(zeroRates)

###############################################################################

    def zeroRate(self,
                 dts: (list, FinDate),
                 freqType: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360):
        ''' Calculation of zero rates with specified frequency. This
        function can return a vector of zero rates given a vector of
        dates so must use Numpy functions. Default frequency is a
        continuously compounded rate. '''

        if isinstance(freqType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if isinstance(dayCountType, FinDayCountTypes) is False:
            raise FinError("Invalid Day Count type.")

        dfs = self.df(dts)
        zeroRates = self._dfToZero(dfs, dts, freqType, dayCountType)

        if isinstance(dts, FinDate):
            return zeroRates[0]
        else:
            return np.array(zeroRates)

        return zeroRates

###############################################################################

    def ccRate(self,
               dts: (list, FinDate), 
               dayCountType: FinDayCountTypes = FinDayCountTypes.SIMPLE):
        ''' Calculation of zero rates with continuous compounding. This
        function can return a vector of cc rates given a vector of
        dates so must use Numpy functions. '''

        ccRates = self.zeroRate(dts, FinFrequencyTypes.CONTINUOUS, dayCountType)
        return ccRates

###############################################################################

    def swapRate(self,
                 effectiveDate: FinDate,
                 maturityDate: (list, FinDate),
                 freqType=FinFrequencyTypes.ANNUAL,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360):
        ''' Calculate the swap rate to maturity date. This is the rate paid by
        a swap that has a price of par today. This is the same as a Libor swap
        rate except that we do not do any business day adjustments. '''

        # Note that this function does not call the FinIborSwap class to
        # calculate the swap rate since that will create a circular dependency.
        # I therefore recreate the actual calculation of the swap rate here.

        if effectiveDate < self._valuationDate:
            raise FinError("Swap starts before the curve valuation date.")

        if isinstance(freqType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if isinstance(freqType, FinFrequencyTypes) is False:
            raise FinError("Invalid Frequency type.")

        if freqType == FinFrequencyTypes.SIMPLE:
            raise FinError("Cannot calculate par rate with simple yield freq.")
        elif freqType == FinFrequencyTypes.CONTINUOUS:
            raise FinError("Cannot calculate par rate with continuous freq.")

        if isinstance(maturityDate, FinDate):
            maturityDates = [maturityDate]
        else:
            maturityDates = maturityDate

        parRates = []

        for maturityDt in maturityDates:

            if maturityDt <= effectiveDate:
                raise FinError("Maturity date is before the swap start date.")

            schedule = FinSchedule(effectiveDate,
                                   maturityDt,
                                   freqType)

            flowDates = schedule._generate()
            flowDates[0] = effectiveDate
            
            dayCounter = FinDayCount(dayCountType)
            prevDt = flowDates[0]
            pv01 = 0.0
            df = 1.0

            for nextDt in flowDates[1:]:                
                df = self.df(nextDt)
                alpha = dayCounter.yearFrac(prevDt, nextDt)[0]
                pv01 += alpha * df
                prevDt = nextDt

            if abs(pv01) < gSmall:
                parRate = 0.0
            else:
                dfStart = self.df(effectiveDate)
                parRate = (dfStart - df) / pv01

            parRates.append(parRate)
 
        parRates = np.array(parRates)

        if isinstance(maturityDate, FinDate):
            return parRates[0]
        else:
            return parRates

###############################################################################

    def df(self,
           dt: (list, FinDate)):
        ''' Function to calculate a discount factor from a date or a
        vector of dates. '''

        times = timesFromDates(dt, self._valuationDate, self._dayCountType)
        dfs = self._df(times)

        if isinstance(dfs, float):
            return dfs
        else:
            return np.array(dfs)

###############################################################################

    def _df(self,
            t: (float, np.ndarray)):
        ''' Hidden function to calculate a discount factor from a time or a
        vector of times. Discourage usage in favour of passing in dates. '''

        if self._interpType is FinInterpTypes.FLAT_FWD_RATES or\
            self._interpType is FinInterpTypes.LINEAR_ZERO_RATES or\
                self._interpType is FinInterpTypes.LINEAR_FWD_RATES:
                    
                    df = interpolate(t,
                                     self._times,
                                     self._dfs,
                                     self._interpType.value)        

        else:

             df = self._interpolator.interpolate(t)

        return df

###############################################################################

    def survProb(self,
                 dt: FinDate):
        ''' This returns a survival probability to a specified date based on
        the assumption that the continuously compounded rate is a default
        hazard rate in which case the survival probability is directly
        analagous to a discount factor. '''

        q = self.df(dt)
        return q

###############################################################################

    def fwd(self,
            dts: FinDate):
        ''' Calculate the continuously compounded forward rate at the forward
        FinDate provided. This is done by perturbing the time by one day only
        and measuring the change in the log of the discount factor divided by
        the time increment dt. I am assuming continuous compounding over the
        one date. '''

        if isinstance(dts, FinDate):
            dtsPlusOneDays = [dts.addDays(1)]
        else:
            dtsPlusOneDays = []
            for dt in dts:
                dtsPlusOneDay = dt.addDays(1)
                dtsPlusOneDays.append(dtsPlusOneDay)

        df1 = self.df(dts)
        df2 = self.df(dtsPlusOneDays)
        dt = 1.0 / gDaysInYear
        fwd = np.log(df1/df2)/(1.0*dt)

        if isinstance(dts, FinDate):
            return fwd[0]
        else:
            return np.array(fwd)

        return fwd

###############################################################################

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

###############################################################################

    def bump(self,
             bumpSize: float):
        ''' Adjust the continuously compounded forward rates by a perturbation
        upward equal to the bump size and return a curve objet with this bumped
        curve. This is used for interest rate risk. '''

        times = self._times.copy()
        values = self._dfs.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            values[i] = values[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._valuationDate,
                                     times,
                                     values,
                                     self._interpType)

        return discCurve

###############################################################################

    def fwdRate(self,
                startDate: (list, FinDate),
                dateOrTenor: (FinDate, str),
                dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360):
        ''' Calculate the forward rate between two forward dates according to
        the specified day count convention. This defaults to Actual 360. The
        first date is specified and the second is given as a date or as a tenor
        which is added to the first date. '''

        if isinstance(startDate, FinDate):
            startDates = []
            startDates.append(startDate)
        elif isinstance(startDate, list):
            startDates = startDate
        else:
            raise FinError("Start date and end date must be same types.")

        dayCount = FinDayCount(dayCountType)

        numDates = len(startDates)
        fwdRates = []
        for i in range(0, numDates):
            dt1 = startDates[i]

            if isinstance(dateOrTenor, str):
                dt2 = dt1.addTenor(dateOrTenor)
            elif isinstance(dateOrTenor, FinDate):
                dt2 = dateOrTenor
            elif isinstance(dateOrTenor, list):
                dt2 = dateOrTenor[i]

            yearFrac = dayCount.yearFrac(dt1, dt2)[0]
            df1 = self.df(dt1)
            df2 = self.df(dt2)
            fwdRate = (df1 / df2 - 1.0) / yearFrac
            fwdRates.append(fwdRate)

        if isinstance(startDate, FinDate):
            return fwdRates[0]
        else:
            return np.array(fwdRates)

###############################################################################

    def __repr__(self):

        s = labelToString("OBJECT TYPE", type(self).__name__)
        numPoints = len(self._dfDates)
        s += labelToString("DATES", "DISCOUNT FACTORS")
        for i in range(0, numPoints):
            s += labelToString("%12s" % self._dfDates[i],
                               "%12.8f" % self._dfs[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
