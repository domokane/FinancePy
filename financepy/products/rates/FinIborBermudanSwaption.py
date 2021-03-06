# ##############################################################################
# # Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# ##############################################################################

import numpy as np

from ...utils.date import Date
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_variables import gDaysInYear
from ...utils.fin_math import ONE_MILLION
from ...utils.FinGlobalTypes import FinExerciseTypes
from ...utils.FinGlobalTypes import FinSwapTypes
from ...utils.FinError import FinError
from ...utils.helper_functions import labelToString, check_argument_types

from ...products.rates.IborSwap import FinIborSwap

from ...models.rates_bdt_tree import FinModelRatesBDT
from ...models.rates_bk_tree import FinModelRatesBK
from ...models.rates_hull_white_tree import FinModelRatesHW

###############################################################################


class FinIborBermudanSwaption(object):
    """ This is the class for the Bermudan-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. This swaption
    can be exercised on any of the fixed coupon payment dates after the first
    exercise date. """

    def __init__(self,
                 settlement_date: Date,
                 exerciseDate: Date,
                 maturity_date: Date,
                 fixed_legType: FinSwapTypes,
                 exerciseType: FinExerciseTypes,
                 fixedCoupon: float,
                 fixedFrequencyType: FrequencyTypes,
                 fixedDayCountType: DayCountTypes,
                 notional=ONE_MILLION,
                 floatFrequencyType=FrequencyTypes.QUARTERLY,
                 floatDayCountType=DayCountTypes.THIRTY_E_360,
                 calendar_type=CalendarTypes.WEEKEND,
                 bus_day_adjust_type=BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type=DateGenRuleTypes.BACKWARD):
        """ Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all of the fixed
        # leg coupon dates until the exercise date inclusive. """

        check_argument_types(self.__init__, locals())

        if settlement_date > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > maturity_date:
            raise FinError("Exercise date must be before swap maturity date")

        if exerciseType == FinExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self._settlement_date = settlement_date
        self._exerciseDate = exerciseDate
        self._maturity_date = maturity_date
        self._fixed_legType = fixed_legType
        self._exerciseType = exerciseType

        self._fixedCoupon = fixedCoupon
        self._fixedFrequencyType = fixedFrequencyType
        self._fixedDayCountType = fixedDayCountType
        self._notional = notional
        self._floatFrequencyType = floatFrequencyType
        self._floatDayCountType = floatDayCountType

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None
        self._cpnTimes = None
        self._cpnFlows = None
        
###############################################################################

    def value(self,
              valuation_date,
              discount_curve,
              model):
        """ Value the Bermudan swaption using the specified model and a
        discount curve. The choices of model are the Hull-White model, the 
        Black-Karasinski model and the Black-Derman-Toy model. """

        floatSpread = 0.0

        # The underlying is a swap in which we pay the fixed amount
        self._underlyingSwap = FinIborSwap(self._exerciseDate,
                                           self._maturity_date,
                                           self._fixed_legType,
                                           self._fixedCoupon,
                                           self._fixedFrequencyType,
                                           self._fixedDayCountType,
                                           self._notional,
                                           floatSpread,
                                           self._floatFrequencyType,
                                           self._floatDayCountType,
                                           self._calendar_type,
                                           self._bus_day_adjust_type,
                                           self._date_gen_rule_type)

        #  I need to do this to generate the fixed leg flows
        self._pv01 = self._underlyingSwap.pv01(valuation_date, discount_curve)

        texp = (self._exerciseDate - valuation_date) / gDaysInYear
        tmat = (self._maturity_date - valuation_date) / gDaysInYear

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpnTimes = [texp]
        cpnFlows = [0.0]

        # The first flow is the expiry date
        num_flows = len(self._underlyingSwap._fixed_leg._payment_dates)

        swap = self._underlyingSwap

        for iFlow in range(0, num_flows):

            flowDate = self._underlyingSwap._fixed_leg._payment_dates[iFlow]

            if flowDate > self._exerciseDate:
                cpnTime = (flowDate - valuation_date) / gDaysInYear
                cpnFlow = swap._fixed_leg._payments[iFlow-1] / self._notional
                cpnTimes.append(cpnTime)
                cpnFlows.append(cpnFlow)

        cpnTimes = np.array(cpnTimes)
        cpnFlows = np.array(cpnFlows)

        self._cpnTimes = cpnTimes
        self._cpnFlows = cpnFlows

        # Allow exercise on coupon dates but control this later for europeans
        self._call_times = cpnTimes

        dfTimes = discount_curve._times
        df_values = discount_curve._dfs

        face_amount = 1.0
        strikePrice = 1.0 # Floating leg is assumed to price at par

        #######################################################################
        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        #######################################################################

        if isinstance(model, FinModelRatesBDT) or isinstance(model, FinModelRatesBK) or isinstance(model, FinModelRatesHW):

            model.buildTree(tmat, dfTimes, df_values)

            v = model.bermudanSwaption(texp,
                                       tmat,
                                       strikePrice,
                                       face_amount,
                                       cpnTimes,
                                       cpnFlows,
                                       self._exerciseType)
        else:

            raise FinError("Invalid model choice for Bermudan Swaption")

        if self._fixed_legType == FinSwapTypes.RECEIVE:
            v = self._notional * v['rec']
        elif self._fixed_legType == FinSwapTypes.PAY:
            v = self._notional * v['pay']

        return v

###############################################################################

    def printSwaptionValue(self):

        print("SWAP PV01:", self._pv01)
        
        n = len(self._cpnTimes)
        
        for i in range(0,n):
            print("CPN TIME: ", self._cpnTimes[i], "FLOW", self._cpnFlows[i])

        n = len(self._call_times)

        for i in range(0,n):
            print("CALL TIME: ", self._call_times[i])

###############################################################################
        
    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("SWAP FIXED LEG TYPE", self._fixed_legType)
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
