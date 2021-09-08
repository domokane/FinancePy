# ##############################################################################
# # Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# ##############################################################################

# from ...utils.FinGlobalVariables import gDaysInYear
# from ...models.HWTree import HWTree
# from ...models.BKTree import BKTree
# from ...utils.FinError import FinError
# from ...utils.FinFrequency import FrequencyTypes
# from ...utils.FinDayCount import FinDayCount
# from ...utils.FinDayCount import DayCountTypes
# from ...products.bonds.Bond import Bond


# from ...utils.Date import Date
# from ...utils.FinCalendar import FinCalendar, CalendarTypes
# from ...utils.FinCalendar import BusDayAdjustTypes, DateGenRuleTypes
# from ...utils.FinSchedule import FinSchedule
# from ...utils.FinMath import ONE_MILLION

# from ...utils.FinHelperFunctions import label_to_string

# from enum import Enum
# import numpy as np

# ###############################################################################


# class FinSwaptionModelTypes(Enum):
#     HULL_WHITE = 1
#     BLACK_KARASINSKI = 2

# ###############################################################################


# class BondOptionTypes(Enum):
#     EUROPEAN_CALL = 1
#     EUROPEAN_PUT = 2
#     AMERICAN_CALL = 3
#     AMERICAN_PUT = 4


# class FinSwaptionTypes(Enum):
#     PAY = 1
#     RECEIVE = 2

# class FinSwaptionExerciseTypes(Enum):
#     EUROPEAN = 1
#     BERMUDAN = 2

# ###############################################################################


# class IborBermudanSwaption
#     """ Class for fixed coupon bonds with embedded call or put optionality. """

#     def __init__(self,
#                  exercise_date,
#                  exercise_type,
#                  maturity_date,
#                  swaptionType,
#                  fixed_coupon,
#                  fixed_frequency_type,
#                  fixed_day_count_type,
#                  notional=ONE_MILLION,
#                  float_frequency_type=FrequencyTypes.QUARTERLY,
#                  float_day_count_type=DayCountTypes.THIRTY_E_360,
#                  calendar_type=CalendarTypes.WEEKEND,
#                  bus_day_adjust_type=BusDayAdjustTypes.FOLLOWING,
#                  date_gen_rule_type=DateGenRuleTypes.BACKWARD):
#         """ Create a Bermudan swaption contract. This is an option to enter
#         into a swap at a fixed coupon on all of the fixed leg coupon dates
#         until the exercise date. """

#         if exercise_date > maturity_date:
#             raise FinError("Exercise date must be before swap maturity date")

#         if exercise_type not in FinSwaptionExerciseTypes:
#             raise FinError("Exercise type must be a FinSwaptionExerciseTypes")

#         if fixed_day_count_type not in DayCountTypes:
#             raise FinError(
#                 "Unknown Fixed Day Count Rule type " +
#                 str(fixed_day_count_type))

#         if fixed_frequency_type not in FrequencyTypes:
#             raise FinError(
#                 "Unknown Fixed Frequency type " +
#                 str(fixed_frequency_type))

#         if float_day_count_type not in DayCountTypes:
#             raise FinError(
#                 "Unknown Float Day Count Rule type " +
#                 str(float_day_count_type))

#         if float_frequency_type not in FrequencyTypes:
#             raise FinError(
#                 "Unknown Float Frequency type " +
#                 str(fixed_frequency_type))

#         if calendar_type not in CalendarTypes:
#             raise FinError("Unknown Calendar type " + str(calendar_type))

#         if bus_day_adjust_type not in BusDayAdjustTypes:
#             raise FinError(
#                 "Unknown Business Day Adjust type " +
#                 str(bus_day_adjust_type))

#         if date_gen_rule_type not in DateGenRuleTypes:
#             raise FinError(
#                 "Unknown Date Gen Rule type " +
#                 str(date_gen_rule_type))

#         self._exercise_date = exercise_date
#         self._maturity_date = maturity_date
#         self._fixed_coupon = fixed_coupon
#         self._fixed_frequency_type = fixed_frequency_type
#         self._fixed_day_count_type = fixed_day_count_type
#         self._notional = notional
#         self._float_frequency_type = float_frequency_type
#         self._float_day_count_type = float_day_count_type

#         self._calendar_type = calendar_type
#         self._bus_day_adjust_type = bus_day_adjust_type
#         self._date_gen_rule_type = date_gen_rule_type

#         self._pv01 = None
#         self._fwdSwapRate = None
#         self._forwardDf = None
#         self._underlyingSwap = None

# ###############################################################################

#     def value(self,
#               valuation_date,
#               discount_curve,
#               model):
#         """ Value the bermuda swaption. This is done using the specified
#         model and a discount curve. """

#         float_spread = 0.0
#         payFixedFlag = True

#         # The underlying is a swap in which we pay the fixed amount
#         swap = IborSwap(self._exercise_date,
#                             self._maturity_date,
#                             self._fixed_coupon,
#                             self._fixed_frequency_type,
#                             self._fixed_day_count_type,
#                             self._notional,
#                             float_spread,
#                             self._float_frequency_type,
#                             self._float_day_count_type,
#                             payFixedFlag,
#                             self._calendar_type,
#                             self._bus_day_adjust_type,
#                             self._date_gen_rule_type)

#         swap.generate_flows()
#         cpn_times = []
#         cpn_amounts = []

#         for iFlow in range(1, len(self._swap._adjustedFixedDates)):
#             flow_date= swap._adjustedFixedDates[iFlow]
#             cpn_time = (flow_date - settlement_date) / gDaysInYear
#             cpn_flow = swap._fixedFlows[iFlow-1] / self._notional
#             cpn_times.append(cpn_time)
#             cpn_amounts.append(cpn_flow)

#         cpn_times = np.array(cpn_times)
#         cpn_amounts = np.array(cpn_amounts)

#         # Generate bond call times and prices
#         call_times = []
#         for dt in self._call_dates:
#             call_time = (dt - settlement_date) / gDaysInYear
#             call_times.append(call_time)
#         call_times = np.array(call_times)
#         call_prices = np.array(self._call_prices)

#         # Generate bond put times and prices
#         if self._swaptionType == IborSwaptionType.PAY:
#             call_price = 100.0
#             putPrice = 1e10
#         else:
#             call_price = 1e10
#             putPrice = 100.0

#         put_times = []
#         for putDate in swap._adjustedFixedDates[1:]:
#             if putDate <= self._exercise_date:
#                 put_time = (putDate - settlement_date) / gDaysInYear
#                 put_times.append(put_time)

#         put_times = np.array(put_times)
#         put_prices = np.array(self._put_prices)

#         maturity_date = self._bond._maturity_date
#         tmat = (maturity_date - settlement_date) / gDaysInYear
#         df_times = discount_curve._times
#         df_values = discount_curve._values

#         face = self._bond._face

#         if type(model) is HWTree:

#             """ We need to build the tree out to the bond maturity date. To be
#             more precise we only need to go out the the last option date but
#             we can do that refinement at a later date. """

#             model.buildTree(tmat, df_times, df_values)
#             v1 = model.callablePuttableBond_Tree(cpn_times, cpn_amounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices, face)
#             model._num_time_steps += 1
#             model.buildTree(tmat, df_times, df_values)
#             v2 = model.callablePuttableBond_Tree(cpn_times, cpn_amounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices, face)
#             model._num_time_steps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

#         elif type(model) == BKTree:

#             """ Because we not have a closed form bond price we need to build
#             the tree out to the bond maturity which is after option expiry. """

#             model.buildTree(tmat, df_times, df_values)
#             v1 = model.callablePuttableBond_Tree(cpn_times, cpn_amounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices,
#                                                  face)
#             model._num_time_steps += 1
#             model.buildTree(tmat, df_times, df_values)
#             v2 = model.callablePuttableBond_Tree(cpn_times, cpn_amounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices,
#                                                  face)
#             model._num_time_steps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
#         else:
#             raise FinError("Unknown model type")

# ###############################################################################

#     def __repr__(self):

#         s = label_to_string("MATURITY DATE", self._maturity_date)
#         s += label_to_string("EXERCISE DATE", self._exercise_date)
#         s += label_to_string("COUPON", self._coupon)
#         s += label_to_string("FREQUENCY", self._freq_type)
#         s += label_to_string("ACCRUAL TYPE", self._accrual_type)
#         s += label_to_string("FACE AMOUNT", self._face)
#         s += label_to_string("CONVERSION RATIO", self._conversion_ratio)
#         s += label_to_string("START CONVERT DATE", self._start_convert_date)

#         for i in range(0, len(self._call_dates)):
#             s += label_to_string("CALL DATE AND PRICE", self._call_dates[i],
#                                self._call_prices[i])

#         for i in range(0, len(self._put_dates)):
#             s += label_to_string("PUT DATE AND PRICE", self._put_dates[i],
#                                self._put_prices[i])

#         return s

# ###############################################################################

#     def print(self):
#         print(self)

# ###############################################################################
