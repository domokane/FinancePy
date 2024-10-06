##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


##############################################################################
# TODO: Allow an index curve to be passed in that is not same as discount curve
# TODO: Extend to allow term structure of volatility
# TODO: Extend to allow two fixed legs in underlying swap
# TODO: Cash settled swaptions
##############################################################################

import numpy as np

from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date

from ...products.rates.ibor_swap import IborSwap

from ...models.black import Black
from ...models.black_shifted import BlackShifted
from ...models.sabr import SABR
from ...models.sabr_shifted import SABRShifted
from ...models.hw_tree import HWTree
from ...models.bk_tree import BKTree
from ...models.bdt_tree import BDTTree

from ...utils.global_types import OptionTypes
from ...utils.global_types import SwapTypes
from ...utils.global_types import FinExerciseTypes

###############################################################################


class IborSwaption:
    """This is the class for the European-style swaption, an option to enter
    into a swap (payer or receiver of the fixed cpn), that starts in the
    future and with a fixed maturity, at a swap rate fixed today."""

    def __init__(
        self,
        settle_dt: Date,
        exercise_dt: Date,
        maturity_dt: Date,
        fixed_leg_type: SwapTypes,
        fixed_cpn: float,
        fixed_freq_type: FrequencyTypes,
        fixed_dc_type: DayCountTypes,
        notional: float = ONE_MILLION,
        float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        float_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed cpn and the details of the fixed and the
        floating leg payment schedules. Bermudan style swaption should be
        priced using the IborBermudanSwaption class."""

        check_argument_types(self.__init__, locals())

        if settle_dt > exercise_dt:
            raise FinError("Settlement date must be before expiry date")

        if exercise_dt > maturity_dt:
            raise FinError("Exercise date must be before swap maturity date")

        self.settle_dt = settle_dt
        self.exercise_dt = exercise_dt
        self.maturity_dt = maturity_dt
        self.fixed_leg_type = fixed_leg_type

        self.notional = notional

        self.fixed_cpn = fixed_cpn
        self.fixed_freq_type = fixed_freq_type
        self.fixed_dc_type = fixed_dc_type
        self.float_freq_type = float_freq_type
        self.float_dc_type = float_dc_type

        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        self.pv01 = None
        self.fwd_swap_rate = None
        self.forward_df = None
        self.underlying_swap = None
        self.swap_type = None

    ###########################################################################

    def value(self, value_dt, discount_curve, model):
        """Valuation of a Ibor European-style swaption using a choice of
        models on a specified valuation date. Models include FinModelBlack,
        FinModelBlackShifted, SABR, SABRShifted, FinModelHW,
        FinModelBK and FinModelBDT. The last two involved a tree-based
        valuation."""

        float_spread = 0.0

        # We create a swap that starts on the exercise date.
        swap = IborSwap(
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

        k = self.fixed_cpn

        # The pv01 is the value of the swap cash flows as of the curve date
        pv01 = swap.pv01(value_dt, discount_curve)

        # We need to calculate the forward swap rate on the swaption exercise
        # date that makes the forward swap worth par including principal
        s = swap.swap_rate(value_dt, discount_curve)

        t_exp = (self.exercise_dt - self.settle_dt) / g_days_in_year
        t_mat = (self.maturity_dt - self.settle_dt) / g_days_in_year

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        #######################################################################
        # For the tree models we need to generate a vector of the cpns
        #######################################################################

        cpn_times = [t_exp]
        cpn_flows = [0.0]

        # The first flow is on the day after the expiry date
        num_flows = len(swap.fixed_leg.payment_dts)

        for i_flow in range(0, num_flows):

            flow_dt = swap.fixed_leg.payment_dts[i_flow]

            # Only flows occurring after option expiry are counted.
            # Flows on the expiry date are not included
            if flow_dt > self.exercise_dt:
                cpn_time = (flow_dt - value_dt) / g_days_in_year
                cpn_flow = swap.fixed_leg.payments[i_flow] / self.notional
                cpn_times.append(cpn_time)
                cpn_flows.append(cpn_flow)

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        if np.any(cpn_times < 0.0):
            raise FinError("No cpn times can be before the value date.")

        strike_price = 1.0
        face_amount = 1.0
        swaption_price = None

        #######################################################################

        if isinstance(model, Black):

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_CALL
                )
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_PUT
                )

        elif isinstance(model, BlackShifted):

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_CALL
                )
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_PUT
                )

        elif isinstance(model, SABR):

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_CALL
                )
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_PUT
                )

        elif isinstance(model, SABRShifted):

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_CALL
                )
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_PUT
                )

        elif isinstance(model, HWTree):

            swaption_px = model.european_bond_option_jamshidian(
                t_exp,
                strike_price,
                face_amount,
                cpn_times,
                cpn_flows,
                df_times,
                df_values,
            )

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaption_px["put"]
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaption_px["call"]
            else:
                raise FinError(
                    "Unknown swaption option type" + str(self.swap_type)
                )

            # Cancel the multiplication at the end below
            swaption_price /= pv01

        elif isinstance(model, BKTree):

            model.build_tree(t_mat, df_times, df_values)
            swaption_px = model.bermudan_swaption(
                t_exp,
                t_mat,
                strike_price,
                face_amount,
                cpn_times,
                cpn_flows,
                FinExerciseTypes.EUROPEAN,
            )

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaption_px["pay"]
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaption_px["rec"]

            swaption_price /= pv01

        elif isinstance(model, BDTTree):

            model.build_tree(t_mat, df_times, df_values)
            swaption_px = model.bermudan_swaption(
                t_exp,
                t_mat,
                strike_price,
                face_amount,
                cpn_times,
                cpn_flows,
                FinExerciseTypes.EUROPEAN,
            )

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaption_px["pay"]
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaption_px["rec"]

            swaption_price /= pv01
        else:
            raise FinError("Unknown swaption model " + str(model))

        self.pv01 = pv01
        self.fwd_swap_rate = s
        self.forward_df = discount_curve.df(self.exercise_dt)
        self.underlying_swap = swap

        # The exchange of cash occurs on the settlement date. However the
        # actual value is that on the specified valuation date which could
        # be the swaption settlement date.
        df_settle = discount_curve.df(self.settle_dt)
        swaption_price = swaption_price * pv01 * self.notional / df_settle
        return swaption_price

    ###########################################################################

    def cash_settled_value(
        self, value_dt: Date, discount_curve, swap_rate: float, model
    ):
        """Valuation of a Ibor European-style swaption using a cash settled
        approach which is a market convention that used Black's model and that
        discounts all of the future payments at a flat swap rate. Note that the
        Black volatility for this valuation should in general not equal the
        Black volatility for the standard arbitrage-free valuation."""

        float_spread = 0.0

        swap = IborSwap(
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

        k = self.fixed_cpn
        s = swap_rate

        pv01 = swap.cash_settled_pv01(
            value_dt, swap_rate, self.fixed_freq_type
        )

        t_exp = (self.exercise_dt - self.settle_dt) / g_days_in_year

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        if isinstance(model, Black):

            if self.fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_CALL
                )
            elif self.fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(
                    s, k, t_exp, df, OptionTypes.EUROPEAN_PUT
                )
        else:
            raise FinError(
                "Cash settled swaptions must be priced using"
                + " Black's model."
            )

        self.fwd_swap_rate = swap_rate
        self.forward_df = discount_curve.df(self.exercise_dt)
        self.underlying_swap = swap
        # The annuity needs to be discounted to today using the correct df
        self.pv01 = pv01 * self.forward_df

        # The exchange of cash occurs on the settlement date but we need to
        # value the swaption on the provided valuation date - which could be
        # the settlement date or may be a different date.
        df_value_dt = discount_curve.df(value_dt)
        swaption_price = (
            swaption_price * self.pv01 * self.notional / df_value_dt
        )
        return swaption_price

    ###########################################################################

    def print_swap_fixed_leg(self):

        if self.underlying_swap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self.underlying_swap.print_fixed_leg_pv()

    ###########################################################################

    def print_swap_float_leg(self):

        if self.underlying_swap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self.underlying_swap.print_float_leg_pv()

    ###########################################################################

    def __repr__(self):
        """Function to allow us to print the swaption details."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SETTLEMENT DATE", self.settle_dt)
        s += label_to_string("EXERCISE DATE", self.exercise_dt)
        s += label_to_string("SWAP FIXED LEG TYPE", str(self.fixed_leg_type))
        s += label_to_string("SWAP MATURITY DATE", self.maturity_dt)
        s += label_to_string("SWAP NOTIONAL", self.notional)
        s += label_to_string("FIXED cpn", self.fixed_cpn * 100)
        s += label_to_string("FIXED FREQUENCY", str(self.fixed_freq_type))
        s += label_to_string("FIXED DAY COUNT", str(self.fixed_dc_type))
        s += label_to_string("FLOAT FREQUENCY", str(self.float_freq_type))
        s += label_to_string("FLOAT DAY COUNT", str(self.float_dc_type))

        if self.pv01 is not None:
            s += label_to_string("PV01", self.pv01)
            s += label_to_string("FWD SWAP RATE", self.fwd_swap_rate * 100)
            s += label_to_string("FWD DF TO EXPIRY", self.forward_df, "")

        return s

    ###########################################################################

    def _print(self):
        """Alternative print method."""

        print(self)


###############################################################################
