# ##############################################################################
# # Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# ##############################################################################

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinOptionTypes import FinOptionExerciseTypes
from ...finutils.FinOptionTypes import FinLiborSwaptionTypes
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

from ...products.libor.FinLiborSwap import FinLiborSwap

from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...models.FinModelRatesBDT import FinModelRatesBDT

###############################################################################


class FinLiborBermudanSwaption(object):
    ''' This is the class for the Bermudan-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. This swaption
    can be exercised on any of the fixed coupon payment dates after the first
    exercise date. '''

    def __init__(self,
                 settlementDate: FinDate,
                 exerciseDate: FinDate,
                 maturityDate: FinDate,
                 swaptionType: FinLiborSwaptionTypes,
                 exerciseType: FinOptionExerciseTypes,
                 fixedCoupon: float,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 notional=ONE_MILLION,
                 floatFrequencyType=FinFrequencyTypes.QUARTERLY,
                 floatDayCountType=FinDayCountTypes.THIRTY_360,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all of the fixed
        # leg coupon dates until the exercise date inclusive. '''

        checkArgumentTypes(self.__init__, locals())

        if settlementDate > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > maturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        self._settlementDate = settlementDate
        self._exerciseDate = exerciseDate
        self._maturityDate = maturityDate
        self._swaptionType = swaptionType
        self._exerciseType = exerciseType

        self._fixedCoupon = fixedCoupon
        self._fixedFrequencyType = fixedFrequencyType
        self._fixedDayCountType = fixedDayCountType
        self._notional = notional
        self._floatFrequencyType = floatFrequencyType
        self._floatDayCountType = floatDayCountType

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None

###############################################################################

    def value(self,
              valuationDate,
              discountCurve,
              model):
        ''' Value the Bermudan swaption using the specified model and a
        discount curve. '''

        floatSpread = 0.0
        payFixedFlag = True

        # The underlying is a swap in which we pay the fixed amount
        swap = FinLiborSwap(self._exerciseDate,
                            self._maturityDate,
                            self._fixedCoupon,
                            self._fixedFrequencyType,
                            self._fixedDayCountType,
                            self._notional,
                            floatSpread,
                            self._floatFrequencyType,
                            self._floatDayCountType,
                            payFixedFlag,
                            self._calendarType,
                            self._busDayAdjustType,
                            self._dateGenRuleType)

        #  I need to do this to generate the fixed leg flows - design issue
        swap.pv01(valuationDate, discountCurve)

        # Set up for Tree
        numFlows = len(swap._adjustedFixedDates)
        cpnTimes = []
        cpnAmounts = []

        texp = (self._exerciseDate - valuationDate) / gDaysInYear
        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        for iFlow in range(1, numFlows):
            flowDate = swap._adjustedFixedDates[iFlow]
            cpnTime = (flowDate - valuationDate) / gDaysInYear
            cpnFlow = swap._fixedFlows[iFlow-1] / self._notional
            cpnTimes.append(cpnTime)
            cpnAmounts.append(cpnFlow)

        cpnTimes = np.array(cpnTimes)
        cpnAmounts = np.array(cpnAmounts)

        callTimes = []
        for cpnTime in cpnTimes:
            if cpnTime >= texp:
                callTimes.append(cpnTime)

#        numCalls = len(callTimes)
#        callTimes = np.array(callTimes)
#        callPrices = np.ones(numCalls)

        dfTimes = discountCurve._times
        dfValues = discountCurve._dfValues

        face = 1.0
        strikePrice = 1.0

        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        if isinstance(model, FinModelRatesHW):

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)

            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)
            model._numTimeSteps -= 1

            if self._swaptionType == FinLiborSwaptionTypes.RECEIVER:
                v = self._notional * (v1['rec'] + v2['rec'])/2.0
                return v
            elif self._swaptionType == FinLiborSwaptionTypes.PAYER:
                v = self._notional * (v1['pay'] + v2['pay'])/2.0
                return v

        elif type(model) == FinModelRatesBK:

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)

            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)
            model._numTimeSteps -= 1

            if self._swaptionType == FinLiborSwaptionTypes.RECEIVER:
                v = self._notional * (v1['rec'] + v2['rec'])/2.0
                return v
            elif self._swaptionType == FinLiborSwaptionTypes.PAYER:
                v = self._notional * (v1['pay'] + v2['pay'])/2.0
                return v

        elif type(model) == FinModelRatesBDT:

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)

            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.bermudanSwaption(texp,
                                        tmat,
                                        strikePrice,
                                        face,
                                        cpnTimes,
                                        cpnAmounts,
                                        self._exerciseType)
            model._numTimeSteps -= 1

            if self._swaptionType == FinLiborSwaptionTypes.RECEIVER:
                v = self._notional * (v1['rec'] + v2['rec'])/2.0
                return v
            elif self._swaptionType == FinLiborSwaptionTypes.PAYER:
                v = self._notional * (v1['pay'] + v2['pay'])/2.0
                return v
        else:
            raise FinError("Unknown model and option combination")

###############################################################################

    def __repr__(self):

        s = labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("SWAPTION TYPE", self._swaptionType)
        s += labelToString("EXERCISE TYPE", self._exerciseType)
        s += labelToString("FIXED COUPON", self._fixedCoupon)
        s += labelToString("FIXED FREQUENCY", self._fixedFrequencyType)
        s += labelToString("FIXED DAYCOUNT TYPE", self._fixedDayCountType)
        s += labelToString("FLOAT FREQUENCY", self._floatFrequencyType)
        s += labelToString("FLOAT DAYCOUNT TYPE", self._floatDayCountType)
        s += labelToString("NOTIONAL", self._notional)

        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
