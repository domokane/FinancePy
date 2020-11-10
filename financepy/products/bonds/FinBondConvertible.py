##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO - MUST ADD ACCRUED INTEREST TO MODEL!!!!

from math import exp, sqrt
from numba import njit
import numpy as np
from typing import List

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes

from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...market.curves.FinInterpolator import FinInterpTypes, _uinterpolate

###############################################################################


@njit(fastmath=True, cache=True)
def _valueConvertible(tmat,
                      faceAmount,
                      couponTimes,
                      couponFlows,
                      callTimes,
                      callPrices,
                      putTimes,
                      putPrices,
                      convRatio,
                      startConvertTime,
                      # Market inputs
                      stockPrice,
                      dfTimes,
                      dfValues,
                      dividendTimes,
                      dividendYields,
                      stockVolatility,
                      creditSpread,
                      recRate,
                      # Tree details
                      numStepsPerYear):

    interp = FinInterpTypes.FLAT_FWD_RATES.value

    if len(couponTimes) > 0:
        if couponTimes[-1] > tmat:
            raise FinError("Coupon after maturity")

    if len(callTimes) > 0:
        if callTimes[-1] > tmat:
            raise FinError("Call times after maturity")

    if len(putTimes) > 0:
        if putTimes[-1] > tmat:
            raise FinError("Put times after maturity")

    if len(dfTimes) > 0:
        if dfTimes[-1] > tmat:
            raise FinError("Discount times after maturity")

    if len(dividendTimes) > 0:
        if dividendTimes[-1] > tmat:
            raise FinError("Dividend times after maturity")

    if creditSpread < 0.0:
        raise FinError("Credit spread negative.")

    if recRate < 0.0 or recRate > 1.0:
        raise FinError("Recovery rate should be between 0 and 1.")

    if stockVolatility < 0.0:
        raise FinError("Stock volatility cannot be negative.")

    if numStepsPerYear < 1:
        raise FinError("Num Steps per year must more than 1.")

    if len(dividendTimes) > 0.0:
        if dividendTimes[-1] > tmat:
            raise FinError("Last dividend is after bond maturity.")

    if recRate > 0.999 or recRate < 0.0:
        raise FinError("Recovery rate must be between 0 and 0.999.")

    numTimes = int(numStepsPerYear * tmat) + 1  # add one for today time 0
    numTimes = numStepsPerYear  # XXXXXXXX!!!!!!!!!!!!!!!!!!!!!

    if numTimes < 5:
        raise FinError("Numsteps must be greater than 5.")

    numLevels = numTimes

    # this is the size of the step
    dt = tmat / (numTimes-1)

    treeTimes = np.linspace(0.0, tmat, numTimes)
    treeDfs = np.zeros(numTimes)
    for i in range(0, numTimes):
        df = _uinterpolate(treeTimes[i], dfTimes, dfValues, interp)
        treeDfs[i] = df

    h = creditSpread/(1.0 - recRate)
    survProb = exp(-h*dt)

    # map coupons onto tree but preserve their present value using risky dfs
    treeFlows = np.zeros(numTimes)
    numCoupons = len(couponTimes)
    for i in range(0, numCoupons):
        flowTime = couponTimes[i]
        n = int(round(flowTime/dt, 0))
        treeTime = treeTimes[n]
        df_flow = _uinterpolate(flowTime, dfTimes, dfValues, interp)
        df_flow *= exp(-h*flowTime)
        df_tree = _uinterpolate(treeTime, dfTimes, dfValues, interp)
        df_tree *= exp(-h*treeTime)
        treeFlows[n] += couponFlows[i] * 1.0 * df_flow / df_tree

    # map call onto tree - must have no calls at high value
    treeCallValue = np.ones(numTimes) * faceAmount * 1000.0
    numCalls = len(callTimes)
    for i in range(0, numCalls):
        callTime = callTimes[i]
        n = int(round(callTime/dt, 0))
        treeCallValue[n] = callPrices[i]

    # map puts onto tree
    treePutValue = np.zeros(numTimes)
    numPuts = len(putTimes)
    for i in range(0, numPuts):
        putTime = putTimes[i]
        n = int(round(putTime/dt, 0))
        treePutValue[n] = putPrices[i]

    # map discrete dividend yields onto tree dates when they are made
    treeDividendYield = np.zeros(numTimes)
    numDividends = len(dividendTimes)
    for i in range(0, numDividends):
        dividendTime = dividendTimes[i]
        n = int(round(dividendTime/dt, 0))
        treeDividendYield[n] = dividendYields[i]

    # Set up the tree of stock prices using a 2D matrix (half the matrix is
    # unused but this may be a cost worth bearing for simpler code. Review.
    treeStockValue = np.zeros(shape=(numTimes, numLevels))
    e = stockVolatility**2 - h
    if e < 0.0:
        raise FinError("Volatility squared minus the hazard rate is negative.")

    u = exp(sqrt(e*dt))
    d = 1.0 / u
    u2 = u*u
    treeStockValue[0, 0] = stockPrice
    for iTime in range(1, numTimes):
        s = treeStockValue[iTime-1, 0] * d
        treeStockValue[iTime, 0] = s

        for iNode in range(1, iTime + 1):
            s = s * u2
            treeStockValue[iTime, iNode] = s

        # we now reduce all stocks by the same yield amount at the same date
        y = treeDividendYield[iTime]
        for iNode in range(0, iTime + 1):
            treeStockValue[iTime, iNode] *= (1.0-y)

    # set up the tree of conversion values. Before allowed to convert the
    # conversion value must be set equal to zero

    treeConvertValue = np.zeros(shape=(numTimes, numLevels))
    for iTime in range(0, numTimes):
        if treeTimes[iTime] >= startConvertTime:
            for iNode in range(0, iTime + 1):
                s = treeStockValue[iTime, iNode]
                treeConvertValue[iTime, iNode] = s * convRatio * 1.0

#    printTree(treeConvertValue)

    treeConvBondValue = np.zeros(shape=(numTimes, numLevels))

    # store probability of up move as a function of time on the tree
    treeProbsUp = np.zeros(numTimes)
    treeProbsDn = np.zeros(numTimes)
    q = 0.0  # we have discrete dividends paid as dividend yields only
    for iTime in range(1, numTimes):
        a = treeDfs[iTime-1]/treeDfs[iTime] * exp(-q*dt)
        treeProbsUp[iTime] = (a - d*survProb) / (u - d)
        treeProbsDn[iTime] = (u*survProb - a) / (u - d)
#        r = log(a)/dt
#        n_min = r*r / stockVolatility / stockVolatility

    if np.any(treeProbsUp > 1.0):
        raise FinError("pUp > 1.0. Increase time steps.")

    ###########################################################################
    # work backwards by first setting values at bond maturity date
    ###########################################################################

    flow = treeFlows[numTimes-1]
    bulletPV = (1.0 + flow) * faceAmount
    for iNode in range(0, numLevels):
        convValue = treeConvertValue[numTimes-1, iNode]
        treeConvBondValue[numTimes-1, iNode] = max(bulletPV, convValue)

    #  begin backward steps from expiry
    for iTime in range(numTimes-2, -1, -1):

        pUp = treeProbsUp[iTime+1]
        pDn = treeProbsDn[iTime+1]
        pDef = 1.0 - survProb
        df = treeDfs[iTime+1]/treeDfs[iTime]
        call = treeCallValue[iTime]
        put = treePutValue[iTime]
        flow = treeFlows[iTime]

        for iNode in range(0, iTime+1):
            futValueUp = treeConvBondValue[iTime+1, iNode+1]
            futValueDn = treeConvBondValue[iTime+1, iNode]
            hold = pUp * futValueUp + pDn * futValueDn  # pUp already embeds Q
            holdPV = df * hold + pDef * df * recRate * faceAmount + flow * faceAmount
            conv = treeConvertValue[iTime, iNode]
            value = min(max(holdPV, conv, put), call)
            treeConvBondValue[iTime, iNode] = value

        bulletPV = df * bulletPV * survProb
        bulletPV += pDef * df * recRate * faceAmount
        bulletPV += flow * faceAmount

    price = treeConvBondValue[0, 0]
    delta = (treeConvBondValue[1, 1] - treeConvBondValue[1, 0]) / \
        (treeStockValue[1, 1] - treeStockValue[1, 0])
    deltaUp = (treeConvBondValue[2, 3] - treeConvBondValue[2, 2]) / \
        (treeStockValue[2, 3] - treeStockValue[2, 2])
    deltaDn = (treeConvBondValue[2, 2] - treeConvBondValue[2, 1]) / \
        (treeStockValue[2, 2] - treeStockValue[2, 1])
    gamma = (deltaUp - deltaDn) / (treeStockValue[1, 1] - treeStockValue[1, 0])
    theta = (treeConvBondValue[2, 2] - treeConvBondValue[0, 0]) / (2.0 * dt)
    results = np.array([price, bulletPV, delta, gamma, theta])
    return results

###############################################################################


class FinBondConvertible(object):
    ''' Class for convertible bonds. These bonds embed rights to call and put
    the bond in return for equity. Until then they are bullet bonds which
    means they have regular coupon payments of a known size that are paid on
    known dates plus a payment of par at maturity. As the options are price
    based, the decision to convert to equity depends on the stock price,
    the credit quality of the issuer and the level of interest rates.'''

    def __init__(self,
                 maturityDate: FinDate,  # bond maturity date
                 coupon: float,  # annual coupon
                 freqType: FinFrequencyTypes,  # coupon frequency type
                 startConvertDate: FinDate,  # conversion starts on this date
                 conversionRatio: float,  # num shares per face of notional
                 callDates: List[FinDate],  # list of call dates
                 callPrices: List[float],  # list of call prices
                 putDates: List[FinDate],  # list of put dates
                 putPrices: List[float],  # list of put prices
                 accrualType: FinDayCountTypes,  # day count type for accrued
                 faceAmount: float = 100.0):  # face amount
        ''' Create FinBondConvertible object by providing the bond Maturity
        date, coupon, frequency type, accrual convention type and then all of
        the details regarding the conversion option including the list of the
        call and put dates and the corresponding list of call and put prices. 
        '''

        checkArgumentTypes(self.__init__, locals())

        if startConvertDate > maturityDate:
            raise FinError("Start convert date is after bond maturity.")

        self._maturityDate = maturityDate
        self._coupon = coupon
        self._accrualType = accrualType
        self._frequency = FinFrequency(freqType)
        self._freqType = freqType

        self._callDates = callDates
        self._callPrices = callPrices

        if len(self._callDates) != len(self._callPrices):
            raise FinError("Call dates and prices not same length.")

        self._putDates = putDates
        self._putPrices = putPrices

        if len(self._putDates) != len(self._putPrices):
            raise FinError("Put dates and prices not same length.")

        if len(putDates) > 0:
            if putDates[-1] > maturityDate:
                raise FinError("Last put is after bond maturity.")

        if len(callDates) > 0:
            if callDates[-1] > maturityDate:
                raise FinError("Last call is after bond maturity.")

        self._startConvertDate = startConvertDate

        if conversionRatio < 0.0:
            raise FinError("Conversion ratio is negative.")

        self._conversionRatio = conversionRatio
        self._faceAmount = faceAmount

        self._settlementDate = FinDate(1, 1, 1900)
        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. '''

        self._accruedInterest = None
        self._accruedDays = 0.0
        self._alpha = 0.0
       
###############################################################################

    def _calculateFlowDates(self,
                            settlementDate: FinDate):
        ''' Determine the convertible bond cashflow payment dates. '''

        # No need to generate flows if settlement date has not changed
        if settlementDate == self._settlementDate:
            return

        self._settlementDate = settlementDate
        calendarType = FinCalendarTypes.NONE
        busDayRuleType = FinBusDayAdjustTypes.NONE
        dateGenRuleType = FinDateGenRuleTypes.BACKWARD

        self._flowDates = FinSchedule(settlementDate,
                                      self._maturityDate,
                                      self._freqType,
                                      calendarType,
                                      busDayRuleType,
                                      dateGenRuleType)._generate()

        self._pcd = self._flowDates[0]
        self._ncd = self._flowDates[1]
        self.calcAccruedInterest(settlementDate)

###############################################################################

    def value(self,
              settlementDate: FinDate,
              stockPrice: float,
              stockVolatility: float,
              dividendDates: List[FinDate],
              dividendYields: List[float],
              discountCurve: FinDiscountCurve,
              creditSpread: float,
              recoveryRate: float = 0.40,
              numStepsPerYear: int = 100):
        '''
        A binomial tree valuation model for a convertible bond that captures
        the embedded equity option due to the existence of a conversion option
        which can be invoked after a specific date.

        The model allows the user to enter a schedule of dividend payment
        dates but the size of the payments must be in yield terms i.e. a known
        percentage of currently unknown future stock price is paid. Not a
        fixed amount. A fixed yield. Following this payment the stock is
        assumed to drop by the size of the dividend payment.

        The model also captures the stock dependent credit risk of the cash
        flows in which the bond price can default at any time with a hazard
        rate implied by the credit spread and an associated recovery rate.
        This is the model proposed by Hull (OFODS 6th edition,.page 522).

        The model captures both the issuer's call schedule which is assumed
        to apply on a list of dates provided by the user, along with a call
        price. It also captures the embedded owner's put schedule of prices.
        '''

        if stockPrice <= 0.0:
            stockPrice = 1e-10  # Avoid overflows in delta calc

        if stockVolatility <= 0.0:
            stockVolatility = 1e-10  # Avoid overflows in delta calc

        self._calculateFlowDates(settlementDate)
        tmat = (self._maturityDate - settlementDate) / gDaysInYear
        if tmat <= 0.0:
            raise FinError("Maturity must not be on or before the value date.")

        # We include time zero in the coupon times and flows
        couponTimes = [0.0]
        couponFlows = [0.0]
        cpn = self._coupon/self._frequency
        for dt in self._flowDates[1:]:
            flowTime = (dt - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

        couponTimes = np.array(couponTimes)
        couponFlows = np.array(couponFlows)

        if np.any(couponTimes < 0.0):
            raise FinError("No coupon times can be before the value date.")

        if np.any(couponTimes > tmat):
            raise FinError("No coupon times can be after the maturity date.")

        callTimes = []
        for dt in self._callDates:
            callTime = (dt - settlementDate) / gDaysInYear
            callTimes.append(callTime)
        callTimes = np.array(callTimes)
        callPrices = np.array(self._callPrices)

        if np.any(callTimes < 0.0):
            raise FinError("No call times can be before the value date.")

        if np.any(callTimes > tmat):
            raise FinError("No call times can be after the maturity date.")

        putTimes = []
        for dt in self._putDates:
            putTime = (dt - settlementDate) / gDaysInYear
            putTimes.append(putTime)
        putTimes = np.array(putTimes)
        putPrices = np.array(self._putPrices)

        if np.any(putTimes > tmat):
            raise FinError("No put times can be after the maturity date.")

        if np.any(putTimes <= 0.0):
            raise FinError("No put times can be on or before value date.")

        if len(dividendYields) != len(dividendDates):
            raise FinError("Number of dividend yields and dates not same.")

        dividendTimes = []
        for dt in dividendDates:
            dividendTime = (dt - settlementDate) / gDaysInYear
            dividendTimes.append(dividendTime)
        dividendTimes = np.array(dividendTimes)
        dividendYields = np.array(dividendYields)

        # If it's before today it starts today
        tconv = (self._startConvertDate - settlementDate) / gDaysInYear
        tconv = max(tconv, 0.0)

        discountFactors = []
        for t in couponTimes:
            df = discountCurve._df(t)
            discountFactors.append(df)

        discountTimes = np.array(couponTimes)
        discountFactors = np.array(discountFactors)

        if testMonotonicity(couponTimes) is False:
            raise FinError("Coupon times not monotonic")

        if testMonotonicity(callTimes) is False:
            raise FinError("Coupon times not monotonic")

        if testMonotonicity(putTimes) is False:
            raise FinError("Coupon times not monotonic")

        if testMonotonicity(discountTimes) is False:
            raise FinError("Coupon times not monotonic")

        if testMonotonicity(dividendTimes) is False:
            raise FinError("Coupon times not monotonic")

        v1 = _valueConvertible(tmat,
                               self._faceAmount,
                               couponTimes,
                               couponFlows,
                               callTimes,
                               callPrices,
                               putTimes,
                               putPrices,
                               self._conversionRatio,
                               tconv,
                               # Market inputs
                               stockPrice,
                               discountTimes,
                               discountFactors,
                               dividendTimes,
                               dividendYields,
                               stockVolatility,
                               creditSpread,
                               recoveryRate,
                               # Tree details
                               numStepsPerYear)

        v2 = _valueConvertible(tmat,
                               self._faceAmount,
                               couponTimes,
                               couponFlows,
                               callTimes,
                               callPrices,
                               putTimes,
                               putPrices,
                               self._conversionRatio,
                               tconv,
                               # Market inputs
                               stockPrice,
                               discountTimes,
                               discountFactors,
                               dividendTimes,
                               dividendYields,
                               stockVolatility,
                               creditSpread,
                               recoveryRate,
                               # Tree details
                               numStepsPerYear + 1)

        cbprice = (v1[0] + v2[0])/2.0
        bond = (v1[1] + v2[1])/2.0
        delta = (v1[2] + v2[2])/2.0
        gamma = (v1[3] + v2[3])/2.0
        theta = (v1[4] + v2[4])/2.0

        results = {"cbprice": cbprice, "bond": bond, "delta": delta,
                   "gamma": gamma, "theta": theta}

        return results

###############################################################################

    def accruedDays(self, 
                    settlementDate: FinDate):
        ''' Calculate number days from previous coupon date to settlement.'''
        self._calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        return settlementDate - self._pcd

###############################################################################

    def calcAccruedInterest(self,
                            settlementDate: FinDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        if settlementDate != self._settlementDate:
            self._calculateFlowDates(settlementDate)

        if len(self._flowDates) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = FinDayCount(self._accrualType)

        (accFactor, num, _) = dc.yearFrac(self._pcd,
                                          settlementDate,
                                          self._ncd, 
                                          self._frequency)

        self._alpha = 1.0 - accFactor * self._frequency

        self._accrued = accFactor * self._faceAmount * self._coupon
        self._accruedDays = num
        return self._accruedInterest

###############################################################################

    def currentYield(self,
                     cleanPrice: float):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._faceAmount / cleanPrice
        return y

###############################################################################

    def __repr__(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount)
        s += labelToString("CONVERSION RATIO", self._conversionRatio)
        s += labelToString("START CONVERT DATE", self._startConvertDate)
        s += labelToString("CALL", "DATES")

        for i in range(0, len(self._callDates)):
            s += labelToString(self._callDates[i],
                               self._callPrices[i])

        s += labelToString("PUT", "DATES")

        for i in range(0, len(self._putDates)):
            s += labelToString(self._putDates[i],
                               self._putPrices[i])

        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# TEST PV OF CASHFLOW MAPPING
#    if 1==0:
#        pv = 0.0
#        for i in range(0, numCoupons):
#            t = couponTimes[i]
#            df = uinterpolate(t, discountTimes, discountFactors, interp)
#            pv += df * couponAmounts[i]
#            print(i, t, couponAmounts[i], df, pv)
#        pv += df
#
#        print("ACTUAL PV",pv)
#
#        pv = 0.0
#        for i in range(0, numTimes):
#            t = treeTimes[i]
#            df = uinterpolate(t, discountTimes, discountFactors, interp)
#            pv += df * treeFlows[i]
#            print(i, t, treeFlows[i], df, pv)
#        pv += df
#
#        print("ACTUAL PV",pv)
###############################################################################
###############################################################################
###############################################################################


###############################################################################


def printTree(array):
    n1, n2 = array.shape
    for i in range(0, n1):
        for j in range(0, n2):
            x = array[j, n1-1-i]
            if x != 0.0:
                print("%10.2f" % array[j, n1-i-1], end="")
            else:
                print("%10s" % '-', end="")
        print("")

###############################################################################
