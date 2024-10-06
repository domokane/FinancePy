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
from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION
from ...utils.global_types import FinExerciseTypes
from ...utils.global_types import SwapTypes
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types

from .ibor_swap import IborSwap

from ...models.bdt_tree import BDTTree
from ...models.bk_tree import BKTree
from ...models.hw_tree import HWTree


###############################################################################


class IborBermudanSwaption:
    """This is the class for the Bermudan-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. This swaption
    can be exercised on any of the fixed coupon payment dates after the first
    exercise date."""

    def __init__(
        self,
        settle_dt: Date,
        exercise_dt: Date,
        maturity_dt: Date,
        fixed_leg_type: SwapTypes,
        exercise_type: FinExerciseTypes,
        fixed_cpn: float,
        fixed_freq_type: FrequencyTypes,
        fixed_dc_type: DayCountTypes,
        notional=ONE_MILLION,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.THIRTY_E_360,
        cal_type=CalendarTypes.WEEKEND,
        bd_type=BusDayAdjustTypes.FOLLOWING,
        dg_type=DateGenRuleTypes.BACKWARD,
    ):
        """Create a Bermudan swaption contract. This is an option to enter
        into a payer or receiver swap at a fixed coupon on all the fixed
        # leg coupon dates until the exercise date inclusive."""

        check_argument_types(self.__init__, locals())

        if settle_dt > exercise_dt:
            raise FinError("Settlement date must be before expiry date")

        if exercise_dt > maturity_dt:
            raise FinError("Exercise date must be before swap maturity date")

        if exercise_type == FinExerciseTypes.AMERICAN:
            raise FinError("American optionality not supported.")

        self.settle_dt = settle_dt
        self.exercise_dt = exercise_dt
        self.maturity_dt = maturity_dt
        self.fixed_leg_type = fixed_leg_type
        self.exercise_type = exercise_type

        self.fixed_cpn = fixed_cpn
        self.fixed_freq_type = fixed_freq_type
        self.fixed_dc_type = fixed_dc_type
        self.notional = notional
        self.float_freq_type = float_freq_type
        self.float_dc_type = float_dc_type

        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        self.pv01 = None
        self.fwd_swap_rate = None
        self.forward_df = None
        self.underlying_swap = None
        self.cpn_times = None
        self.cpn_flows = None
        self.call_times = None

    ###########################################################################

    def value(self, value_dt, discount_curve, model):
        """Value the Bermudan swaption using the specified model and a
        discount curve. The choices of model are the Hull-White model, the
        Black-Karasinski model and the Black-Derman-Toy model."""

        float_spread = 0.0

        # The underlying is a swap in which we pay the fixed amount
        self.underlying_swap = IborSwap(
            self.exercise_dt,
            self.maturity_dt,
            self.fixed_leg_type,
            self.fixed_cpn,
            self.fixed_freq_type,
            self.fixed_dc_type,
            self.notional,
            float_spread,
            self.float_freq_type,
            self.float_dc_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

        #  I need to do this to generate the fixed leg flows
        self.pv01 = self.underlying_swap.pv01(value_dt, discount_curve)

        t_exp = (self.exercise_dt - value_dt) / g_days_in_year
        t_mat = (self.maturity_dt - value_dt) / g_days_in_year

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpn_times = [t_exp]
        cpn_flows = [0.0]

        # The first flow is the expiry date
        num_flows = len(self.underlying_swap.fixed_leg.payment_dts)

        swap = self.underlying_swap

        for i_flow in range(0, num_flows):

            flow_dt = self.underlying_swap.fixed_leg.payment_dts[i_flow]

            if flow_dt > self.exercise_dt:
                cpn_time = (flow_dt - value_dt) / g_days_in_year
                cpn_flow = swap.fixed_leg.payments[i_flow - 1] / self.notional
                cpn_times.append(cpn_time)
                cpn_flows.append(cpn_flow)

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        self.cpn_times = cpn_times
        self.cpn_flows = cpn_flows

        # Allow exercise on coupon dates but control this later for europeans
        self.call_times = cpn_times

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        face_amount = 1.0
        strike_price = 1.0  # Floating leg is assumed to price at par

        #######################################################################
        # For both models, the tree needs to extend out to maturity because of
        # the multi-callable nature of the Bermudan Swaption
        #######################################################################

        if (
            isinstance(model, BDTTree)
            or isinstance(model, BKTree)
            or isinstance(model, HWTree)
        ):

            model.build_tree(t_mat, df_times, df_values)

            v = model.bermudan_swaption(
                t_exp,
                t_mat,
                strike_price,
                face_amount,
                cpn_times,
                cpn_flows,
                self.exercise_type,
            )
        else:

            raise FinError("Invalid model choice for Bermudan Swaption")

        if self.fixed_leg_type == SwapTypes.RECEIVE:
            v = self.notional * v["rec"]
        elif self.fixed_leg_type == SwapTypes.PAY:
            v = self.notional * v["pay"]

        return v

    ###########################################################################

    def print_swaption_value(self):

        print("SWAP PV01:", self.pv01)

        n = len(self.cpn_times)

        for i in range(0, n):
            print("CPN TIME: ", self.cpn_times[i], "FLOW", self.cpn_flows[i])

        n = len(self.call_times)

        for i in range(0, n):
            print("CALL TIME: ", self.call_times[i])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXERCISE DATE", self.exercise_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("SWAP FIXED LEG TYPE", self.fixed_leg_type)
        s += label_to_string("EXERCISE TYPE", self.exercise_type)
        s += label_to_string("FIXED COUPON", self.fixed_cpn)
        s += label_to_string("FIXED FREQUENCY", self.fixed_freq_type)
        s += label_to_string("FIXED DAYCOUNT TYPE", self.fixed_dc_type)
        s += label_to_string("FLOAT FREQUENCY", self.float_freq_type)
        s += label_to_string("FLOAT DAYCOUNT TYPE", self.float_dc_type)
        s += label_to_string("NOTIONAL", self.notional)
        return s

    ###########################################################################

    def _print(self):
        print(self)


###############################################################################
