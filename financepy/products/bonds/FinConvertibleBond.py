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
which the bond price tends to zero as the stock price also tends to zero.
NOT TRUE YET

The model captures both the issuer's call schedule which is assumed to apply
on a list of dates provided by the user, alonmg with a call price. It also
captures the embedded owner's put schedule of prices. '''


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
                     noConversionTime,
                     # Market inputs
                     stockPrice,
                     discountTimes,
                     discountFactors,
                     dividendTimes,
                     dividendYields,
                     stockVolatility,
                     creditSpread,
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

    numTimes = numStepsPerYear+1 # int(numStepsPerYear * tmat) + 1 # add one for today time 0
    numLevels = numTimes

    # this is the size of the step
    dt = tmat / (numTimes-1)
    treeTimes = np.linspace(0.0, tmat, numTimes)
    treeDfs = np.zeros(numTimes)
    for i in range(0, numTimes):
        df = interpolate(treeTimes[i], discountTimes, discountFactors, interp)
        treeDfs[i] = df

    # map coupons onto tree but preserve their present value
    treeFlows = np.zeros(numTimes)
    numCoupons = len(couponTimes)
    for i in range(0, numCoupons):
        flowTime = couponTimes[i]
        n = int(round(flowTime/dt, 0))
        treeTime = treeTimes[n]
        df_flow = interpolate(flowTime, discountTimes, discountFactors, interp)
        df_tree = interpolate(treeTime, discountTimes, discountFactors, interp)
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

    # map dividends onto tree
    treeDividendYield = np.zeros(numTimes)
    numDividends = len(dividendTimes)
    for i in range(0, numDividends):
        dividendTime = dividendTimes[i]
        n = int(round(dividendTime/dt, 0))
        treeDividendYield[n] = dividendYields[i]

    # Set up the tree of stock prices using a 2D matrix (half the matrix is
    # unused but this may be a cost worth bearing for simpler code. Review.
    treeStockValue = np.zeros(shape=(numTimes, numLevels))
    treeStockValue[0, 0] = stockPrice
    u = exp(stockVolatility * sqrt(dt))
    d = 1.0 / u
    u2 = u*u
    sLow = stockPrice
    for iTime in range(1, numTimes):
        sLow *= d
        s = sLow
        y = treeDividendYield[iTime]
        for iNode in range(0, iTime + 1):
            treeStockValue[iTime, iNode] = s
            s = s * u2 * (1.0-y)

    # set up the tree of conversion values. Before allowed to convert the
    # conversion value must be set equal to zero
    treeConvertValue = np.zeros(shape=(numTimes, numLevels))
    for iTime in range(0, numTimes):
        if treeTimes[iTime] >= noConversionTime:
            for iNode in range(0, iTime + 1):
                s = treeStockValue[iTime, iNode]
                treeConvertValue[iTime, iNode] = s * convRatio * 100.0 / face

    treeConvBondValue = np.zeros(shape=(numTimes, numLevels))

    # store probability of up move as a function of time on the tree
    treeProbs = np.zeros(numTimes)
    q = 0.0  # we have discrete dividends paid as dividend yields only
    for iTime in range(1, numTimes):
        a = treeDfs[iTime-1]/treeDfs[iTime] * exp(-q*dt)
        treeProbs[iTime] = (a - d) / (u - d)

    ###########################################################################
    # work backwards by first setting values at expiry date

    flow = treeFlows[numTimes-1]
    bulletPV = 100.0 + flow

    for iNode in range(0, numLevels):
        convValue = treeConvertValue[numTimes-1, iNode]
        treeConvBondValue[numTimes-1, iNode] = max(100.0, convValue) + flow

#    print("BOND1")
#    print(treeConvBondValue)
    #  begin backward steps from expiry
    for iTime in range(numTimes-2, -1, -1):

        pUp = treeProbs[iTime+1]
        pDn = 1.0 - pUp
        df = treeDfs[iTime+1]/treeDfs[iTime]
        call = treeCallValue[iTime]
        put = treePutValue[iTime]
        flow = treeFlows[iTime]
        df_credit = exp(-dt*creditSpread)
#        print(iTime,pUp,df)

        for iNode in range(0, iTime+1):
            futureValueUp = treeConvBondValue[iTime+1, iNode+1]
            futureValueDn = treeConvBondValue[iTime+1, iNode]
            futureValue = pUp * futureValueUp + pDn * futureValueDn
            hold = df * futureValue * df_credit
#            print("X",iTime,df,iNode,futureValueUp,futureValueDn,futureValue,hold)
            conv = treeConvBondValue[iTime][iNode]
            value = min(max(hold, conv, put), call) + flow
            treeConvBondValue[iTime][iNode] = value

        bulletPV *= df * df_credit
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
                 maturityDate,
                 coupon,
                 frequencyType,
                 startConvertDate,
                 conversionRatio,  # number of shares per 100 of face
                 callDates,
                 callPrices,
                 putDates,
                 putPrices,
                 accrualType,
                 face=100.0):
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
              numStepsPerYear):

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
                             # Tree details
                             numStepsPerYear)

        return v

##########################################################################

    def print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print("Maturity Date:", self._maturityDate)
        print("Coupon:", self._coupon)
        print("Frequency:", self._frequencyType)
        print("Accrual:", self._accrualType)

##########################################################################
