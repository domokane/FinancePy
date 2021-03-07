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
from ...utils.helpers import labelToString, check_argument_types
from ...utils.date import Date

from ...products.rates.IborSwap import FinIborSwap

from ...models.black import FinModelBlack
from ...models.black_shifted import FinModelBlackShifted
from ...models.sabr import FinModelSABR
from ...models.sabr_shifted import FinModelSABRShifted
from ...models.rates_hull_white_tree import FinModelRatesHW
from ...models.rates_bk_tree import FinModelRatesBK
from ...models.rates_bdt_tree import FinModelRatesBDT

from ...utils.global_types import FinOptionTypes
from ...utils.global_types import FinSwapTypes
from ...utils.global_types import FinExerciseTypes

###############################################################################


class FinIborSwaption():
    """ This is the class for the European-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. """

    def __init__(self,
                 settlement_date: Date,
                 exerciseDate: Date,
                 maturity_date: Date,
                 fixed_legType: FinSwapTypes,
                 fixedCoupon: float,
                 fixedFrequencyType: FrequencyTypes,
                 fixedDayCountType: DayCountTypes,
                 notional: float = ONE_MILLION,
                 floatFrequencyType: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 floatDayCountType: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed coupon and the details of the fixed and the
        floating leg payment schedules. Bermudan style swaption should be
        priced using the FinIborBermudanSwaption class. """

        check_argument_types(self.__init__, locals())

        if settlement_date > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > maturity_date:
            raise FinError("Exercise date must be before swap maturity date")

        self._settlement_date = settlement_date
        self._exerciseDate = exerciseDate
        self._maturity_date = maturity_date
        self._fixed_legType = fixed_legType

        self._notional = notional

        self._fixedCoupon = fixedCoupon
        self._fixedFrequencyType = fixedFrequencyType
        self._fixedDayCountType = fixedDayCountType
        self._floatFrequencyType = floatFrequencyType
        self._floatDayCountType = floatDayCountType

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
        FinModelBlackShifted, FinModelSABR, FinModelSABRShifted, FinModelHW,
        FinModelBK and FinModelBDT. The last two involved a tree-based
        valuation. """

        floatSpread = 0.0

        # We create a swap that starts on the exercise date.
        swap = FinIborSwap(self._exerciseDate,
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

        k = self._fixedCoupon

        # The pv01 is the value of the swap cash flows as of the curve date
        pv01 = swap.pv01(valuation_date, discount_curve)

        # We need to calculate the forward swap rate on the swaption exercise
        # date that makes the forward swap worth par including principal
        s = swap.swap_rate(valuation_date, discount_curve)

        texp = (self._exerciseDate - self._settlement_date) / gDaysInYear
        tmat = (self._maturity_date - self._settlement_date) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpn_times = [texp]
        cpnFlows = [0.0]

        # The first flow is on the day after the expiry date
        num_flows = len(swap._fixed_leg._payment_dates)

        for iFlow in range(0, num_flows):
            
            flowDate = swap._fixed_leg._payment_dates[iFlow]

            # Only flows occurring after option expiry are counted. 
            # Flows on the expiry date are not included
            if flowDate > self._exerciseDate:
                cpn_time = (flowDate - valuation_date) / gDaysInYear
                cpnFlow = swap._fixed_leg._payments[iFlow] / self._notional
                cpn_times.append(cpn_time)
                cpnFlows.append(cpnFlow)

        cpn_times = np.array(cpn_times)
        cpnFlows = np.array(cpnFlows)

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        if np.any(cpn_times < 0.0):
            raise FinError("No coupon times can be before the value date.")

        strike_price = 1.0
        face_amount = 1.0

        #######################################################################

        if isinstance(model, FinModelBlack):

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBlackShifted):

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABR):

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABRShifted):

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelRatesHW):

            swaptionPx = model.europeanBondOptionJamshidian(texp,
                                                  strike_price,
                                                  face_amount,
                                                  cpn_times,
                                                  cpnFlows,
                                                  df_times,
                                                  df_values)

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['put']
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['call']
            else:
                raise FinError("Unknown swaption option type" +
                               str(self._swapType))

            # Cancel the multiplication at the end below
            swaptionPrice /= pv01

        elif isinstance(model, FinModelRatesBK):

            model.buildTree(tmat, df_times, df_values)
            swaptionPx = model.bermudanSwaption(texp,
                                                tmat,
                                                strike_price,
                                                face_amount,
                                                cpn_times,
                                                cpnFlows,
                                                FinExerciseTypes.EUROPEAN)

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['pay']
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['rec']

            swaptionPrice /= pv01

        elif isinstance(model, FinModelRatesBDT):

            model.buildTree(tmat, df_times, df_values)
            swaptionPx = model.bermudanSwaption(texp,
                                                tmat,
                                                strike_price,
                                                face_amount,
                                                cpn_times,
                                                cpnFlows,
                                                FinExerciseTypes.EUROPEAN)

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['pay']
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['rec']

            swaptionPrice /= pv01
        else:
            raise FinError("Unknown swaption model " + str(model))

        self._pv01 = pv01
        self._fwdSwapRate = s
        self._forwardDf = discount_curve.df(self._exerciseDate)
        self._underlyingSwap = swap

        # The exchange of cash occurs on the settlement date. However the
        # actual value is that on the specified valuation date which could
        # be the swaption settlement date.
        dfSettlement = discount_curve.df(self._settlement_date)
        swaptionPrice = swaptionPrice * pv01 * self._notional / dfSettlement
        return swaptionPrice

###############################################################################

    def cashSettledValue(self,
                         valuation_date: Date,
                         discount_curve,
                         swap_rate: float,
                         model):
        """ Valuation of a Ibor European-style swaption using a cash settled
        approach which is a market convention that used Black's model and that
        discounts all of the future payments at a flat swap rate. Note that the
        Black volatility for this valuation should in general not equal the
        Black volatility for the standard arbitrage-free valuation. """

        floatSpread = 0.0

        swap = FinIborSwap(self._exerciseDate,
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

        k = self._fixedCoupon
        s = swap_rate

        pv01 = swap.cashSettledPV01(valuation_date,
                                    swap_rate,
                                    self._fixedFrequencyType)

        texp = (self._exerciseDate - self._settlement_date) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        if isinstance(model, FinModelBlack):

            if self._fixed_legType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixed_legType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)
        else:
            raise FinError("Cash settled swaptions must be priced using"
                           + " Black's model.")

        self._fwdSwapRate = swap_rate
        self._forwardDf = discount_curve.df(self._exerciseDate)
        self._underlyingSwap = swap
        # The annuity needs to be discounted to today using the correct df
        self._pv01 = pv01 * self._forwardDf

        # The exchange of cash occurs on the settlement date but we need to 
        # value the swaption on the provided valuation date - which could be
        # the settlement date or may be a different date.
        dfValuation = discount_curve.df(valuation_date)
        swaptionPrice = swaptionPrice * self._pv01 * self._notional / dfValuation
        return swaptionPrice

###############################################################################

    def printSwapFixedLeg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.printFixedLegPV()

###############################################################################

    def printSwapFloatLeg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.printFloatLegPV()

###############################################################################

    def __repr__(self):
        """ Function to allow us to print the swaption details. """

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("SETTLEMENT DATE", self._settlement_date)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("SWAP FIXED LEG TYPE", str(self._fixed_legType))
        s += labelToString("SWAP MATURITY DATE", self._maturity_date)
        s += labelToString("SWAP NOTIONAL", self._notional)
        s += labelToString("FIXED COUPON", self._fixedCoupon * 100)
        s += labelToString("FIXED FREQUENCY", str(self._fixedFrequencyType))
        s += labelToString("FIXED DAY COUNT", str(self._fixedDayCountType))
        s += labelToString("FLOAT FREQUENCY", str(self._floatFrequencyType))
        s += labelToString("FLOAT DAY COUNT", str(self._floatDayCountType))

        if self._pv01 is not None:
            s += labelToString("PV01", self._pv01)
            s += labelToString("FWD SWAP RATE", self._fwdSwapRate*100)
            s += labelToString("FWD DF TO EXPIRY", self._forwardDf, "")

        return s

###############################################################################

    def _print(self):
        """ Alternative print method. """

        print(self)

###############################################################################
