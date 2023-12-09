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
                 settle_dt: Date,
                 exercise_dt: Date,
                 maturity_dt: Date,
                 fixed_leg_type: SwapTypes,
                 exercise_type: FinExerciseTypes,
                 fixed_coupon: float,
                 fixed_freq_type: FrequencyTypes,
                 fixed_dc_type: DayCountTypes,
                 notional=ONE_MILLION,
                 float_freq_type=FrequencyTypes.QUARTERLY,
                 float_dc_type=DayCountTypes.THIRTY_E_360,
                 cal_type=CalendarTypes.WEEKEND,
                 bd_type=BusDayAdjustTypes.FOLLOWING,
                 dg_type=DateGenRuleTypes.BACKWARD):
        """ Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all the fixed
        # leg coupon dates until the exercise date inclusive. """

        check_argument_types(self.__init__, locals())

        if settle_dt > exercise_dt:
            raise FinError("Settlement date must be before expiry date")

        if exercise_dt > maturity_dt:
            raise FinError("Exercise date must be before swap maturity date")

        if exercise_type == FinExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self._settle_dt = settle_dt
        self._exercise_dt = exercise_dt
        self._maturity_dt = maturity_dt
        self._fixed_leg_type = fixed_leg_type
        self._exercise_type = exercise_type

        self._fixed_coupon = fixed_coupon
        self._fixed_freq_type = fixed_freq_type
        self._fixed_dc_type = fixed_dc_type
        self._notional = notional
        self._float_freq_type = float_freq_type
        self._float_dc_type = float_dc_type

        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type

        self._pv01 = None
        self._fwd_swap_rate = None
        self._forward_df = None
        self._underlying_swap = None
        self._cpn_times = None
        self._cpn_flows = None

###############################################################################

    def value(self,
              value_dt,
              discount_curve,
              model):
        """ Value the Bermudan swaption using the specified model and a
        discount curve. The choices of model are the Hull-White model, the
        Black-Karasinski model and the Black-Derman-Toy model. """

        float_spread = 0.0

        # The underlying is a swap in which we pay the fixed amount
        self._underlying_swap = IborSwap(self._exercise_dt,
                                         self._maturity_dt,
                                         self._fixed_leg_type,
                                         self._fixed_coupon,
                                         self._fixed_freq_type,
                                         self._fixed_dc_type,
                                         self._notional,
                                         float_spread,
                                         self._float_freq_type,
                                         self._float_dc_type,
                                         self._cal_type,
                                         self._bd_type,
                                         self._dg_type)

        #  I need to do this to generate the fixed leg flows
        self._pv01 = self._underlying_swap.pv01(value_dt, discount_curve)

        t_exp = (self._exercise_dt - value_dt) / gDaysInYear
        t_mat = (self._maturity_dt - value_dt) / gDaysInYear

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpn_times = [t_exp]
        cpn_flows = [0.0]

        # The first flow is the expiry date
        num_flows = len(self._underlying_swap._fixed_leg._payment_dts)

        swap = self._underlying_swap

        for iFlow in range(0, num_flows):

            flow_dt = self._underlying_swap._fixed_leg._payment_dts[iFlow]

            if flow_dt > self._exercise_dt:
                cpn_time = (flow_dt - value_dt) / gDaysInYear
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

            model.build_tree(t_mat, df_times, df_values)

            v = model.bermudan_swaption(t_exp,
                                        t_mat,
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
        s += label_to_string("EXERCISE DATE", self._exercise_dt)
        s += label_to_string("MATURITY DATE", self._maturity_dt)
        s += label_to_string("SWAP FIXED LEG TYPE", self._fixed_leg_type)
        s += label_to_string("EXERCISE TYPE", self._exercise_type)
        s += label_to_string("FIXED COUPON", self._fixed_coupon)
        s += label_to_string("FIXED FREQUENCY", self._fixed_freq_type)
        s += label_to_string("FIXED DAYCOUNT TYPE", self._fixed_dc_type)
        s += label_to_string("FLOAT FREQUENCY", self._float_freq_type)
        s += label_to_string("FLOAT DAYCOUNT TYPE", self._float_dc_type)
        s += label_to_string("NOTIONAL", self._notional)
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
