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
from ...finutils.FinOptionTypes import FinLiborSwapTypes
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

from ...products.libor.FinLiborSwap import FinLiborSwap

from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...models.FinModelRatesBDT import FinModelRatesBDT
from ...models.FinModelBlack import FinModelBlack
from ...models.FinModelBlackShifted import FinModelBlackShifted
from ...models.FinModelSABR import FinModelSABR
from ...models.FinModelSABRShifted import FinModelSABRShifted

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
                 swapMaturityDate: FinDate,
                 swapType: FinLiborSwapTypes,
                 exerciseType: FinOptionExerciseTypes,
                 fixedCoupon: float,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 notional=ONE_MILLION,
                 floatFrequencyType=FinFrequencyTypes.QUARTERLY,
                 floatDayCountType=FinDayCountTypes.THIRTY_E_360,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all of the fixed
        # leg coupon dates until the exercise date inclusive. '''

        checkArgumentTypes(self.__init__, locals())

        if settlementDate > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > swapMaturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        if exerciseType == FinOptionExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self._settlementDate = settlementDate
        self._exerciseDate = exerciseDate
        self._swapMaturityDate = swapMaturityDate
        self._swapType = swapType
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
        discount curve. The choices of model are the Hull-White model, the 
        Black-Karasinski model and the Black-Derman-Toy model. '''

        floatSpread = 0.0

        # The underlying is a swap in which we pay the fixed amount
        swap = FinLiborSwap(self._exerciseDate,
                            self._swapMaturityDate,
                            self._swapType,
                            self._fixedCoupon,
                            self._fixedFrequencyType,
                            self._fixedDayCountType,
                            self._notional,
                            floatSpread,
                            self._floatFrequencyType,
                            self._floatDayCountType,
                            self._calendarType,
                            self._busDayAdjustType,
                            self._dateGenRuleType)

        #  I need to do this to generate the fixed leg flows
        swap.pv01(valuationDate, discountCurve)

        texp = (self._exerciseDate - valuationDate) / gDaysInYear
        tmat = (self._swapMaturityDate - valuationDate) / gDaysInYear

        cpnTimes = [texp]
        cpnFlows = [0.0]

        numFlows = len(swap._adjustedFixedDates)

        # The first flow is always the PCD

        for iFlow in range(1, numFlows):

            pcd = swap._adjustedFixedDates[iFlow-1]
            ncd = swap._adjustedFixedDates[iFlow]

            if ncd > valuationDate:

                if len(cpnTimes) == 0:
                    cpnTime = (pcd - valuationDate) / gDaysInYear
                    cpnFlow = swap._fixedFlows[iFlow-1] / self._notional
                    cpnTimes.append(cpnTime)
                    cpnFlows.append(cpnFlow)
                    
                cpnTime = (ncd - valuationDate) / gDaysInYear
                cpnFlow = swap._fixedFlows[iFlow-1] / self._notional
                cpnTimes.append(cpnTime)
                cpnFlows.append(cpnFlow)

        cpnTimes = np.array(cpnTimes)
        cpnFlows = np.array(cpnFlows)

        callTimes = []
        for cpnTime in cpnTimes:
            if cpnTime >= texp:
                callTimes.append(cpnTime)

        dfTimes = discountCurve._times
        dfValues = discountCurve._dfValues

        faceAmount = 1.0
        strikePrice = 1.0

        #######################################################################
        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        #######################################################################

        if isinstance(model, FinModelBlack) or \
           isinstance(model, FinModelBlackShifted) or \
             isinstance(model, FinModelSABR) or \
               isinstance(model, FinModelSABRShifted):
                  raise FinError("Model is not valid for Bermudan Swaptions")

        model.buildTree(tmat, dfTimes, dfValues)
        v = model.bermudanSwaption(texp,
                                   tmat,
                                   strikePrice,
                                   faceAmount,
                                   cpnTimes,
                                   cpnFlows,
                                   self._exerciseType)

        if self._swapType == FinLiborSwapTypes.RECEIVER:
            v = self._notional * v['rec']
            return v
        elif self._swapType == FinLiborSwapTypes.PAYER:
            v = self._notional * v['pay']
            return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("MATURITY DATE", self._swapMaturityDate)
        s += labelToString("SWAP TYPE", self._swapType)
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
