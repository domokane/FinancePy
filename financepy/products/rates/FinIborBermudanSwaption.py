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
from ...finutils.FinGlobalTypes import FinExerciseTypes
from ...finutils.FinGlobalTypes import FinSwapTypes
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

from ...products.rates.FinIborSwap import FinIborSwap

from ...models.FinModelRatesBDT import FinModelRatesBDT
from ...models.FinModelRatesBK import FinModelRatesBK
from ...models.FinModelRatesHW import FinModelRatesHW

###############################################################################


class FinIborBermudanSwaption(object):
    ''' This is the class for the Bermudan-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. This swaption
    can be exercised on any of the fixed coupon payment dates after the first
    exercise date. '''

    def __init__(self,
                 settlementDate: FinDate,
                 exerciseDate: FinDate,
                 maturityDate: FinDate,
                 fixedLegType: FinSwapTypes,
                 exerciseType: FinExerciseTypes,
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

        if exerciseDate > maturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        if exerciseType == FinExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self._settlementDate = settlementDate
        self._exerciseDate = exerciseDate
        self._maturityDate = maturityDate
        self._fixedLegType = fixedLegType
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
        self._cpnTimes = None
        self._cpnFlows = None
        
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
        self._underlyingSwap = FinIborSwap(self._exerciseDate,
                                           self._maturityDate,
                                           self._fixedLegType,
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
        self._pv01 = self._underlyingSwap.pv01(valuationDate, discountCurve)

        texp = (self._exerciseDate - valuationDate) / gDaysInYear
        tmat = (self._maturityDate - valuationDate) / gDaysInYear

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpnTimes = [texp]
        cpnFlows = [0.0]

        # The first flow is the expiry date
        numFlows = len(self._underlyingSwap._fixedLeg._paymentDates)

        swap = self._underlyingSwap

        for iFlow in range(0, numFlows):

            flowDate = self._underlyingSwap._fixedLeg._paymentDates[iFlow]

            if flowDate > self._exerciseDate:
                cpnTime = (flowDate - valuationDate) / gDaysInYear
                cpnFlow = swap._fixedLeg._payments[iFlow-1] / self._notional
                cpnTimes.append(cpnTime)
                cpnFlows.append(cpnFlow)

        cpnTimes = np.array(cpnTimes)
        cpnFlows = np.array(cpnFlows)

        self._cpnTimes = cpnTimes
        self._cpnFlows = cpnFlows

        # Allow exercise on coupon dates but control this later for europeans
        self._callTimes = cpnTimes

        dfTimes = discountCurve._times
        dfValues = discountCurve._dfs

        faceAmount = 1.0
        strikePrice = 1.0 # Floating leg is assumed to price at par

        #######################################################################
        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        #######################################################################

        if isinstance(model, FinModelRatesBDT) or isinstance(model, FinModelRatesBK) or isinstance(model, FinModelRatesHW):

            model.buildTree(tmat, dfTimes, dfValues)

            v = model.bermudanSwaption(texp,
                                       tmat,
                                       strikePrice,
                                       faceAmount,
                                       cpnTimes,
                                       cpnFlows,
                                       self._exerciseType)
        else:

            raise FinError("Invalid model choice for Bermudan Swaption")

        if self._fixedLegType == FinSwapTypes.RECEIVE:
            v = self._notional * v['rec']
        elif self._fixedLegType == FinSwapTypes.PAY:
            v = self._notional * v['pay']

        return v

###############################################################################

    def printSwaptionValue(self):

        print("SWAP PV01:", self._pv01)
        
        n = len(self._cpnTimes)
        
        for i in range(0,n):
            print("CPN TIME: ", self._cpnTimes[i], "FLOW", self._cpnFlows[i])

        n = len(self._callTimes)

        for i in range(0,n):
            print("CALL TIME: ", self._callTimes[i])

###############################################################################
        
    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("SWAP FIXED LEG TYPE", self._fixedLegType)
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
