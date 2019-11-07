# -*- coding: utf-8 -*-

# TODO

from math import exp, sqrt
from numba import njit, jit
import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinInterpolate import FinInterpMethods, interpolate
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinGlobalVariables import gDaysInYear
from .FinBond import FinBondAccruedTypes

###############################################################################


''' A binomial tree valuation model for a convertible bond that captures the
embedded equity option due to the existence of a conversion option which can
be invoked after a specific date.

The model allows the user to enter a schedule of dividend payment dates but
the size of the payments must be in yield terms i.e. a known percentage of
currently unknown future stock price is paid. Not a fixed amount. A fixed
yield. Following this payment the stock is assumed to drop by the size of
the dividend payment.

The model also captures the stock dependent credit risk of the cash flows in
which the bond price can default at any time with a hazard rate implied by
the credit spread and an associated recovery rate.

The model captures both the issuer's call schedule which is assumed to apply
on a list of dates provided by the user, along with a call price. It also
captures the embedded owner's put schedule of prices.
'''


@njit(fastmath=True, cache=True)
def valueConvertible(tmat,
                     face,
                     couponTimes,
                     couponAmounts,
                     callTimes,
                     callPrices,
                     putTimes,
                     putPrices,
                     convRatio,
                     startConvertTime,
                     # Market inputs
                     stockPrice,
                     discountTimes,
                     discountFactors,
                     dividendTimes,
                     dividendYields,
                     stockVolatility,
                     creditSpread,
                     recRate,
                     # Tree details
                     numStepsPerYear):

    interp = FinInterpMethods.FLAT_FORWARDS.value

    if len(couponTimes) > 0:
        if couponTimes[-1] > tmat:
            raise ValueError("Coupon after maturity")

    if len(callTimes) > 0:
        if callTimes[-1] > tmat:
            raise ValueError("Coupon after maturity")

    if len(putTimes) > 0:
        if putTimes[-1] > tmat:
            raise ValueError("Coupon after maturity")

    if len(discountTimes) > 0:
        if discountTimes[-1] > tmat:
            raise ValueError("Coupon after maturity")

    if len(dividendTimes) > 0:
        if dividendTimes[-1] > tmat:
            raise ValueError("Coupon after maturity")

    if creditSpread < 0.0:
        raise ValueError("Credit spread negative.")

    if recRate < 0.0 or recRate > 1.0:
        raise ValueError("Recovery rate should be between 0 and 1.")

    if stockVolatility < 0.0:
        raise ValueError("Stock volatility cannot be negative.")

    if numStepsPerYear < 1:
        raise ValueError("Num Steps per year must more than 1.")

    if putTimes[-1] > tmat:
        raise ValueError("Last put is after bond maturity.")

    if callTimes[-1] > tmat:
        raise ValueError("Last call is after bond maturity.")

    if dividendTimes[-1] > tmat:
        raise ValueError("Last dividend is after bond maturity.")

    numTimes = int(numStepsPerYear * tmat) + 1  # add one for today time 0
    numLevels = numTimes

    # this is the size of the step
    dt = tmat / (numTimes-1)
    treeTimes = np.linspace(0.0, tmat, numTimes)
    treeDfs = np.zeros(numTimes)
    for i in range(0, numTimes):
        df = interpolate(treeTimes[i], discountTimes, discountFactors, interp)
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
        df_flow = interpolate(flowTime, discountTimes, discountFactors, interp)
        df_flow *= exp(-h*flowTime)
        df_tree = interpolate(treeTime, discountTimes, discountFactors, interp)
        df_tree *= exp(-h*treeTime)
        treeFlows[n] += couponAmounts[i] * 100.0 * df_flow / df_tree

    # map call onto tree - must have no calls at high value
    treeCallValue = np.ones(numTimes) * 9999.0
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
                treeConvertValue[iTime, iNode] = s * convRatio * 100.0 / face

    treeConvBondValue = np.zeros(shape=(numTimes, numLevels))

    # store probability of up move as a function of time on the tree
    treeProbsUp = np.zeros(numTimes)
    treeProbsDn = np.zeros(numTimes)
    q = 0.0  # we have discrete dividends paid as dividend yields only
    for iTime in range(1, numTimes):
        a = treeDfs[iTime-1]/treeDfs[iTime] * exp(-q*dt)
        treeProbsUp[iTime] = (a - d*survProb) / (u - d)
        treeProbsDn[iTime] = (u*survProb - a) / (u - d)

    ###########################################################################
    # work backwards by first setting values at bond maturity date
    ###########################################################################

    flow = treeFlows[numTimes-1]
    bulletPV = 100.0 + flow

    for iNode in range(0, numLevels):
        convValue = treeConvertValue[numTimes-1, iNode]
        treeConvBondValue[numTimes-1, iNode] = max(100.0 + flow, convValue)

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
            hold = df * (pUp * futValueUp + pDn * futValueDn + pDef * recRate)
            conv = treeConvBondValue[iTime][iNode]
            value = min(max(hold, conv, put), call) + flow
            treeConvBondValue[iTime][iNode] = value

        bulletPV = bulletPV * df * survProb + (1.0 - survProb) * recRate
        bulletPV += flow

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


class FinConvertibleBond(object):
    ''' Class for convertible bonds. These bonds embed rights to call and put
    the bond in return for equity. Until then they are bullet bonds which
    means they have regular coupon payments of a known size that are paid on
    known dates plus a payment of par at maturity. As the options are price
    based, the decision to convert to equity depends on the stock price,
    the credit quality of the issuer and the level of interest rates.'''

    def __init__(self,
                 maturityDate, # bond maturity date
                 coupon, # annual coupon
                 frequencyType, # coupon frequency type
                 startConvertDate, # date after which conversion is possible
                 conversionRatio,  # number of shares per face of notional
                 callDates, # list of call dates
                 callPrices, # list of call prices
                 putDates, # list of put dates
                 putPrices, # list of put prices
                 accrualType, # day count type for accrued interest
                 face=100.0 # face amount
                 ):
        ''' Create FinBond object by providing Maturity Date, Frequency,
        coupon and the accrual convention type. '''

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Invalid Frequency:" + str(frequencyType))
            return

        if accrualType not in FinBondAccruedTypes:
            raise FinError("Unknown Bond Accrued Convention type " +
                           str(accrualType))

        self._maturityDate = maturityDate
        self._coupon = coupon
        self._accrualType = accrualType
        self._frequency = FinFrequency(frequencyType)
        self._callDates = callDates
        self._callPrices = callPrices
        self._putDates = putDates
        self._putPrices = putPrices
        self._startConvertDate = startConvertDate
        self._conversionRatio = conversionRatio
        self._face = face
        self._settlementDate = None
        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. '''

##########################################################################

    def calculateFlowDates(self, settlementDate):
        ''' Determine the bond cashflow payment dates. '''

        self._settlementDate = settlementDate
        self._flowDates = []

        nextDate = self._maturityDate
        months = int(12 / self._frequency)

        while nextDate._excelDate > settlementDate._excelDate:
            self._flowDates.append(nextDate)
            nextDate = nextDate.addMonths(-months)

        self._flowDates.reverse()
        return

##########################################################################

    def value(self,
              settlementDate,
              stockPrice,
              stockVolatility,
              dividendDates,
              dividendYields,
              discountCurve,
              creditSpread,
              recoveryRate,
              numStepsPerYear):
        ''' Valuation of a convertible bond using binomial tree model with a
        credit spread component as proposed by Hull (6th edition,.page 522) '''

        self.calculateFlowDates(settlementDate)
        tmat = (self._maturityDate - settlementDate) / gDaysInYear

        # We include time zero in the coupon times and flows
        couponTimes = [0.0]
        couponFlows = [0.0]
        cpn = self._coupon/self._frequency
        for dt in self._flowDates:
            flowTime = (dt - settlementDate) / gDaysInYear
            couponTimes.append(flowTime)
            couponFlows.append(cpn)

        couponTimes = np.array(couponTimes)
        couponFlows = np.array(couponFlows)

        callTimes = []
        for dt in self._callDates:
            callTime = (dt - settlementDate) / gDaysInYear
            callTimes.append(callTime)
        callTimes = np.array(callTimes)
        callPrices = np.array(self._callPrices)

        putTimes = []
        for dt in self._putDates:
            putTime = (dt - settlementDate) / gDaysInYear
            putTimes.append(putTime)
        putTimes = np.array(putTimes)
        putPrices = np.array(self._putPrices)

        if len(dividendYields) != len(dividendDates):
            raise ValueError("Number of dividend yields and dates not same.")

        dividendTimes = []
        for dt in dividendDates:
            dividendTime = (dt - settlementDate) / gDaysInYear
            dividendTimes.append(dividendTime)
        dividendTimes = np.array(dividendTimes)
        dividendYields = np.array(dividendYields)

        tconv = (self._startConvertDate - settlementDate) / gDaysInYear

        discountFactors = []
        for t in couponTimes:
            df = discountCurve.df(t)
            discountFactors.append(df)

        discountTimes = np.array(couponTimes)
        discountFactors = np.array(discountFactors)

        if testMonotonicity(couponTimes) is False:
            return ValueError("Coupon times not monotonic")

        if testMonotonicity(callTimes) is False:
            return ValueError("Coupon times not monotonic")

        if testMonotonicity(putTimes) is False:
            return ValueError("Coupon times not monotonic")

        if testMonotonicity(discountTimes) is False:
            return ValueError("Coupon times not monotonic")

        if testMonotonicity(dividendTimes) is False:
            return ValueError("Coupon times not monotonic")

        v = valueConvertible(tmat,
                             self._face,
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

        return v

##########################################################################

    def accruedDays(self, settlementDate):
        ''' Calculate number days from previous coupon date to settlement.'''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        return settlementDate - self.pcd(settlementDate)

##########################################################################

    def pcd(self, settlementDate):
        ''' Determine the previous coupon date before the settlement date. '''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12 / self._frequency)
        ncd = self._flowDates[0]
        pcd = FinDate.addMonths(ncd, -months)

        return pcd

##########################################################################

    def accruedInterest(self, settlementDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12 / self._frequency)
        ncd = self._flowDates[0]
        pcd = FinDate.addMonths(ncd, -months)

        f = self._frequency
        dt1 = pcd
        dt2 = settlementDate

        if dt1 > dt2:
            raise FinError("First coupon must precede settlement date.")

        d1 = dt1._d
        d2 = dt2._d
        m1 = dt1._m
        m2 = dt2._m
        y1 = dt1._y
        y2 = dt2._y

        if self._accrualType == FinBondAccruedTypes.THIRTY_360:
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            daysInPeriod = 360.0 / f
            accFactor = dayDiff / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.THIRTY_360_BOND:
            d1 = min(d1, 30)
            if d1 == 31 or d1 == 30:
                d2 = min(d2, 30)
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            daysInPeriod = 360.0 / f
            accFactor = dayDiff / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_360:
            accFactor = (dt2 - dt1) / 360.0
            daysInPeriod = 360.0 / f
            accFactor = (dt2 - dt1) / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_365:
            daysInPeriod = 365.0 / f
            accFactor = (dt2 - dt1) / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_ACT:
            daysInPeriod = ncd - pcd
            accFactor = (settlementDate - pcd) / daysInPeriod
        else:
            raise FinError(str(self._type) + " is not one of BondAccrualTypes")

        flow = self._coupon / f
        accrued = accFactor * flow * self._face
        return accrued

##########################################################################

    def currentYield(self, cleanPrice):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._face / cleanPrice
        return y

##########################################################################

    def print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print("Maturity Date:", self._maturityDate)
        print("Coupon:", self._coupon)
        print("Frequency:", self._frequencyType)
        print("Accrual:", self._accrualType)

##########################################################################
