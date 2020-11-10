##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...products.bonds.FinBond import FinBond

from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from enum import Enum
import numpy as np
from typing import List


###############################################################################
# TODO: Make it possible to specify start and end of American Callable/Puttable
###############################################################################


class FinBondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class FinBondOptionTypes(Enum):
    EUROPEAN_CALL = 1
    EUROPEAN_PUT = 2
    AMERICAN_CALL = 3
    AMERICAN_PUT = 4


###############################################################################


class FinBondEmbeddedOption(object):
    ''' Class for fixed coupon bonds with embedded call or put optionality. '''

    def __init__(self,
                 issueDate: FinDate,
                 maturityDate: FinDate,  # FinDate
                 coupon: float,  # Annualised coupon - 0.03 = 3.00%
                 freqType: FinFrequencyTypes,
                 accrualType: FinDayCountTypes,
                 callDates: List[FinDate],
                 callPrices: List[float],
                 putDates: List[FinDate],
                 putPrices: List[float],
                 faceAmount: float = 100.0):
        ''' Create a FinBondEmbeddedOption object with a maturity date, coupon
        and all of the bond inputs. '''

        checkArgumentTypes(self.__init__, locals())

        self._issueDate = issueDate
        self._maturityDate = maturityDate
        self._coupon = coupon
        self._freqType = freqType
        self._accrualType = accrualType

        self._bond = FinBond(issueDate,
                             maturityDate,
                             coupon,
                             freqType,
                             accrualType,
                             faceAmount)

        # Validate call and put schedules
        for dt in callDates:
            if dt > self._maturityDate:
                raise FinError("Call date after bond maturity date")

        if len(callDates) > 0:
            dtprev = callDates[0]
            for dt in callDates[1:]:
                if dt <= dtprev:
                    raise FinError("Call dates not increasing")
                else:
                    dtprev = dt

        for dt in putDates:
            if dt > self._maturityDate:
                raise FinError("Put date after bond maturity date")

        if len(putDates) > 0:
            dtprev = putDates[0]
            for dt in putDates[1:]:
                if dt <= dtprev:
                    raise FinError("Put dates not increasing")
                else:
                    dtprev = dt

        for px in callPrices:
            if px < 0.0:
                raise FinError("Call price must be positive.")

        for px in putPrices:
            if px < 0.0:
                raise FinError("Put price must be positive.")

        if len(callDates) != len(callPrices):
            raise FinError("Number of call dates and call prices not the same")

        if len(putDates) != len(putPrices):
            raise FinError("Number of put dates and put prices not the same")

        self._callDates = callDates
        self._callPrices = callPrices
        self._putDates = putDates
        self._putPrices = putPrices
        self._faceAmount = faceAmount
        self._bond._calculateFlowDates()

###############################################################################

    def value(self,
              settlementDate: FinDate,
              discountCurve: FinDiscountCurve,
              model):
        ''' Value the bond that settles on the specified date that can have
        both embedded call and put options. This is done using the specified
        model and a discount curve. '''

        # Generate bond coupon flow schedule
        cpn = self._bond._coupon/self._bond._frequency
        cpnTimes = []
        cpnAmounts = []

        for flowDate in self._bond._flowDates[1:]:
            if flowDate > settlementDate:
                cpnTime = (flowDate - settlementDate) / gDaysInYear
                cpnTimes.append(cpnTime)
                cpnAmounts.append(cpn)

        cpnTimes = np.array(cpnTimes)
        cpnAmounts = np.array(cpnAmounts)

        # Generate bond call times and prices
        callTimes = []
        for dt in self._callDates:
            if dt > settlementDate:
                callTime = (dt - settlementDate) / gDaysInYear
                callTimes.append(callTime)
        callTimes = np.array(callTimes)
        callPrices = np.array(self._callPrices)

        # Generate bond put times and prices
        putTimes = []
        for dt in self._putDates:
            if dt > settlementDate:
                putTime = (dt - settlementDate) / gDaysInYear
                putTimes.append(putTime)
        putTimes = np.array(putTimes)
        putPrices = np.array(self._putPrices)

        maturityDate = self._bond._maturityDate
        tmat = (maturityDate - settlementDate) / gDaysInYear
        dfTimes = discountCurve._times
        dfValues = discountCurve._dfs

        faceAmount = self._bond._faceAmount

        if isinstance(model, FinModelRatesHW):

            ''' We need to build the tree out to the bond maturity date. To be
            more precise we only need to go out the the last option date but
            we can do that refinement at a later date. '''

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 callTimes, callPrices,
                                                 putTimes, putPrices,
                                                 faceAmount)
            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 callTimes, callPrices,
                                                 putTimes, putPrices,
                                                 faceAmount)
            model._numTimeSteps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

        elif isinstance(model, FinModelRatesBK):

            ''' Because we not have a closed form bond price we need to build
            the tree out to the bond maturity which is after option expiry. '''

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 callTimes, callPrices,
                                                 putTimes, putPrices,
                                                 faceAmount)
            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 callTimes, callPrices,
                                                 putTimes, putPrices,
                                                 faceAmount)
            model._numTimeSteps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
        else:
            raise FinError("Unknown model type")

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issueDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount)

        s += labelToString("NUM CALL DATES", len(self._callDates))
        for i in range(0, len(self._callDates)):
            s += "%12s %12.6f\n" % (self._callDates[i], self._callPrices[i])

        s += labelToString("NUM PUT DATES", len(self._putDates))
        for i in range(0, len(self._putDates)):
            s += "%12s %12.6f\n" % (self._putDates[i], self._putPrices[i])

        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
