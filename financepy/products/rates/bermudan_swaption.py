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
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...utils.global_types import FinExerciseTypes
from ...utils.global_types import SwapTypes
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

from ...products.rates.ibor_swap import IborSwap

from ...models.bdt_tree import BDTTree
from ...models.bk_tree import BKTree
from ...models.hw_tree import HWTree

###############################################################################


class IborBermudanSwaption:
    """ This is the class for the Bermudan-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. This swaption
    can be exercised on any of the fixed coupon payment dates after the first
    exercise date. """

    def __init__(self,
                 settlement_date: Date,
                 exercise_date: Date,
                 maturity_date: Date,
                 fixed_leg_type: SwapTypes,
                 exercise_type: FinExerciseTypes,
                 fixed_coupon: float,
                 fixed_frequency_type: FrequencyTypes,
                 fixed_day_count_type: DayCountTypes,
                 notional=ONE_MILLION,
                 float_frequency_type=FrequencyTypes.QUARTERLY,
                 float_day_count_type=DayCountTypes.THIRTY_E_360,
                 calendar_type=CalendarTypes.WEEKEND,
                 bus_day_adjust_type=BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type=DateGenRuleTypes.BACKWARD):
        """ Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all of the fixed
        # leg coupon dates until the exercise date inclusive. """

        check_argument_types(self.__init__, locals())

        if settlement_date > exercise_date:
            raise FinError("Settlement date must be before expiry date")

        if exercise_date > maturity_date:
            raise FinError("Exercise date must be before swap maturity date")

        if exercise_type == FinExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self._settlement_date = settlement_date
        self._exercise_date = exercise_date
        self._maturity_date = maturity_date
        self._fixed_leg_type = fixed_leg_type
        self._exercise_type = exercise_type

        self._fixed_coupon = fixed_coupon
        self._fixed_frequency_type = fixed_frequency_type
        self._fixed_day_count_type = fixed_day_count_type
        self._notional = notional
        self._float_frequency_type = float_frequency_type
        self._float_day_count_type = float_day_count_type

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None
        self._cpn_times = None
        self._cpn_flows = None

###############################################################################

    def value(self,
              valuation_date,
              discount_curve,
              model):
        """ Value the Bermudan swaption using the specified model and a
        discount curve. The choices of model are the Hull-White model, the 
        Black-Karasinski model and the Black-Derman-Toy model. """

        float_spread = 0.0

        # The underlying is a swap in which we pay the fixed amount
        self._underlyingSwap = IborSwap(self._exercise_date,
                                        self._maturity_date,
                                        self._fixed_leg_type,
                                        self._fixed_coupon,
                                        self._fixed_frequency_type,
                                        self._fixed_day_count_type,
                                        self._notional,
                                        float_spread,
                                        self._float_frequency_type,
                                        self._float_day_count_type,
                                        self._calendar_type,
                                        self._bus_day_adjust_type,
                                        self._date_gen_rule_type)

        #  I need to do this to generate the fixed leg flows
        self._pv01 = self._underlyingSwap.pv01(valuation_date, discount_curve)

        texp = (self._exercise_date - valuation_date) / gDaysInYear
        tmat = (self._maturity_date - valuation_date) / gDaysInYear

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpn_times = [texp]
        cpn_flows = [0.0]

        # The first flow is the expiry date
        num_flows = len(self._underlyingSwap._fixed_leg._payment_dates)

        swap = self._underlyingSwap

        for iFlow in range(0, num_flows):

            flow_date = self._underlyingSwap._fixed_leg._payment_dates[iFlow]

            if flow_date > self._exercise_date:
                cpn_time = (flow_date - valuation_date) / gDaysInYear
                cpn_flow = swap._fixed_leg._payments[iFlow-1] / self._notional
                cpn_times.append(cpn_time)
                cpn_flows.append(cpn_flow)

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        self._cpn_times = cpn_times
        self._cpn_flows = cpn_flows

        # Allow exercise on coupon dates but control this later for europeans
        self._call_times = cpn_times

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        face_amount = 1.0
        strike_price = 1.0  # Floating leg is assumed to price at par

        #######################################################################
        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        #######################################################################

        if isinstance(model, BDTTree) or isinstance(model, BKTree) or isinstance(model, HWTree):

            model.build_tree(tmat, df_times, df_values)

            v = model.bermudan_swaption(texp,
                                        tmat,
                                        strike_price,
                                        face_amount,
                                        cpn_times,
                                        cpn_flows,
                                        self._exercise_type)
        else:

            raise FinError("Invalid model choice for Bermudan Swaption")

        if self._fixed_leg_type == SwapTypes.RECEIVE:
            v = self._notional * v['rec']
        elif self._fixed_leg_type == SwapTypes.PAY:
            v = self._notional * v['pay']

        return v

###############################################################################

    def print_swaption_value(self):

        print("SWAP PV01:", self._pv01)

        n = len(self._cpn_times)

        for i in range(0, n):
            print("CPN TIME: ", self._cpn_times[i], "FLOW", self._cpn_flows[i])

        n = len(self._call_times)

        for i in range(0, n):
            print("CALL TIME: ", self._call_times[i])

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXERCISE DATE", self._exercise_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("SWAP FIXED LEG TYPE", self._fixed_leg_type)
        s += label_to_string("EXERCISE TYPE", self._exercise_type)
        s += label_to_string("FIXED COUPON", self._fixed_coupon)
        s += label_to_string("FIXED FREQUENCY", self._fixed_frequency_type)
        s += label_to_string("FIXED DAYCOUNT TYPE", self._fixed_day_count_type)
        s += label_to_string("FLOAT FREQUENCY", self._float_frequency_type)
        s += label_to_string("FLOAT DAYCOUNT TYPE", self._float_day_count_type)
        s += label_to_string("NOTIONAL", self._notional)
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
