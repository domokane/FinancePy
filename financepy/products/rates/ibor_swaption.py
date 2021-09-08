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
from ...utils.global_vars import gDaysInYear
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


class IborSwaption():
    """ This is the class for the European-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. """

    def __init__(self,
                 settlement_date: Date,
                 exercise_date: Date,
                 maturity_date: Date,
                 fixed_leg_type: SwapTypes,
                 fixed_coupon: float,
                 fixed_frequency_type: FrequencyTypes,
                 fixed_day_count_type: DayCountTypes,
                 notional: float = ONE_MILLION,
                 float_frequency_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 float_day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed coupon and the details of the fixed and the
        floating leg payment schedules. Bermudan style swaption should be
        priced using the IborBermudanSwaption class. """

        check_argument_types(self.__init__, locals())

        if settlement_date > exercise_date:
            raise FinError("Settlement date must be before expiry date")

        if exercise_date > maturity_date:
            raise FinError("Exercise date must be before swap maturity date")

        self._settlement_date = settlement_date
        self._exercise_date = exercise_date
        self._maturity_date = maturity_date
        self._fixed_leg_type = fixed_leg_type

        self._notional = notional

        self._fixed_coupon = fixed_coupon
        self._fixed_frequency_type = fixed_frequency_type
        self._fixed_day_count_type = fixed_day_count_type
        self._float_frequency_type = float_frequency_type
        self._float_day_count_type = float_day_count_type

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None

###############################################################################

    def value(self,
              valuation_date,
              discount_curve,
              model):
        """ Valuation of a Ibor European-style swaption using a choice of
        models on a specified valuation date. Models include FinModelBlack,
        FinModelBlackShifted, SABR, SABRShifted, FinModelHW,
        FinModelBK and FinModelBDT. The last two involved a tree-based
        valuation. """

        float_spread = 0.0

        # We create a swap that starts on the exercise date.
        swap = IborSwap(self._exercise_date,
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

        k = self._fixed_coupon

        # The pv01 is the value of the swap cash flows as of the curve date
        pv01 = swap.pv01(valuation_date, discount_curve)

        # We need to calculate the forward swap rate on the swaption exercise
        # date that makes the forward swap worth par including principal
        s = swap.swap_rate(valuation_date, discount_curve)

        texp = (self._exercise_date - self._settlement_date) / gDaysInYear
        tmat = (self._maturity_date - self._settlement_date) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpn_times = [texp]
        cpn_flows = [0.0]

        # The first flow is on the day after the expiry date
        num_flows = len(swap._fixed_leg._payment_dates)

        for iFlow in range(0, num_flows):

            flow_date = swap._fixed_leg._payment_dates[iFlow]

            # Only flows occurring after option expiry are counted.
            # Flows on the expiry date are not included
            if flow_date > self._exercise_date:
                cpn_time = (flow_date - valuation_date) / gDaysInYear
                cpn_flow = swap._fixed_leg._payments[iFlow] / self._notional
                cpn_times.append(cpn_time)
                cpn_flows.append(cpn_flow)

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        if np.any(cpn_times < 0.0):
            raise FinError("No coupon times can be before the value date.")

        strike_price = 1.0
        face_amount = 1.0
        swaption_price = None

        #######################################################################

        if isinstance(model, Black):

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_CALL)
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, BlackShifted):

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_CALL)
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABR):

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_CALL)
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABRShifted):

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_CALL)
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, HWTree):

            swaptionPx = model.european_bond_option_jamshidian(texp,
                                                               strike_price,
                                                               face_amount,
                                                               cpn_times,
                                                               cpn_flows,
                                                               df_times,
                                                               df_values)

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaptionPx['put']
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaptionPx['call']
            else:
                raise FinError("Unknown swaption option type" +
                               str(self._swapType))

            # Cancel the multiplication at the end below
            swaption_price /= pv01

        elif isinstance(model, BKTree):

            model.build_tree(tmat, df_times, df_values)
            swaptionPx = model.bermudan_swaption(texp,
                                                 tmat,
                                                 strike_price,
                                                 face_amount,
                                                 cpn_times,
                                                 cpn_flows,
                                                 FinExerciseTypes.EUROPEAN)

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaptionPx['pay']
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaptionPx['rec']

            swaption_price /= pv01

        elif isinstance(model, BDTTree):

            model.build_tree(tmat, df_times, df_values)
            swaptionPx = model.bermudan_swaption(texp,
                                                 tmat,
                                                 strike_price,
                                                 face_amount,
                                                 cpn_times,
                                                 cpn_flows,
                                                 FinExerciseTypes.EUROPEAN)

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = swaptionPx['pay']
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = swaptionPx['rec']

            swaption_price /= pv01
        else:
            raise FinError("Unknown swaption model " + str(model))

        self._pv01 = pv01
        self._fwdSwapRate = s
        self._forwardDf = discount_curve.df(self._exercise_date)
        self._underlyingSwap = swap

        # The exchange of cash occurs on the settlement date. However the
        # actual value is that on the specified valuation date which could
        # be the swaption settlement date.
        dfSettlement = discount_curve.df(self._settlement_date)
        swaption_price = swaption_price * pv01 * self._notional / dfSettlement
        return swaption_price

###############################################################################

    def cash_settled_value(self,
                           valuation_date: Date,
                           discount_curve,
                           swap_rate: float,
                           model):
        """ Valuation of a Ibor European-style swaption using a cash settled
        approach which is a market convention that used Black's model and that
        discounts all of the future payments at a flat swap rate. Note that the
        Black volatility for this valuation should in general not equal the
        Black volatility for the standard arbitrage-free valuation. """

        float_spread = 0.0

        swap = IborSwap(self._exercise_date,
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

        k = self._fixed_coupon
        s = swap_rate

        pv01 = swap.cash_settled_pv01(valuation_date,
                                      swap_rate,
                                      self._fixed_frequency_type)

        texp = (self._exercise_date - self._settlement_date) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        if isinstance(model, Black):

            if self._fixed_leg_type == SwapTypes.PAY:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_CALL)
            elif self._fixed_leg_type == SwapTypes.RECEIVE:
                swaption_price = model.value(s, k, texp, df,
                                             OptionTypes.EUROPEAN_PUT)
        else:
            raise FinError("Cash settled swaptions must be priced using"
                           + " Black's model.")

        self._fwdSwapRate = swap_rate
        self._forwardDf = discount_curve.df(self._exercise_date)
        self._underlyingSwap = swap
        # The annuity needs to be discounted to today using the correct df
        self._pv01 = pv01 * self._forwardDf

        # The exchange of cash occurs on the settlement date but we need to
        # value the swaption on the provided valuation date - which could be
        # the settlement date or may be a different date.
        dfValuation = discount_curve.df(valuation_date)
        swaption_price = swaption_price * self._pv01 * self._notional / dfValuation
        return swaption_price

###############################################################################

    def print_swap_fixed_leg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.print_fixed_leg_pv()

###############################################################################

    def print_swap_float_leg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.print_float_leg_pv()

###############################################################################

    def __repr__(self):
        """ Function to allow us to print the swaption details. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("SETTLEMENT DATE", self._settlement_date)
        s += label_to_string("EXERCISE DATE", self._exercise_date)
        s += label_to_string("SWAP FIXED LEG TYPE", str(self._fixed_leg_type))
        s += label_to_string("SWAP MATURITY DATE", self._maturity_date)
        s += label_to_string("SWAP NOTIONAL", self._notional)
        s += label_to_string("FIXED COUPON", self._fixed_coupon * 100)
        s += label_to_string("FIXED FREQUENCY",
                             str(self._fixed_frequency_type))
        s += label_to_string("FIXED DAY COUNT",
                             str(self._fixed_day_count_type))
        s += label_to_string("FLOAT FREQUENCY",
                             str(self._float_frequency_type))
        s += label_to_string("FLOAT DAY COUNT",
                             str(self._float_day_count_type))

        if self._pv01 is not None:
            s += label_to_string("PV01", self._pv01)
            s += label_to_string("FWD SWAP RATE", self._fwdSwapRate*100)
            s += label_to_string("FWD DF TO EXPIRY", self._forwardDf, "")

        return s

###############################################################################

    def _print(self):
        """ Alternative print method. """

        print(self)

###############################################################################
