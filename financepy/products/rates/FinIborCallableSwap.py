# ##############################################################################
# # Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# ##############################################################################

# from ...utils.FinGlobalVariables import gDaysInYear
# from ...models.FinModelRatesHW import FinModelRatesHW
# from ...models.FinModelRatesBK import FinModelRatesBK
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

# from ...utils.FinHelperFunctions import labelToString

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


# class FinIborBermudanSwaption(object):
#     """ Class for fixed coupon bonds with embedded call or put optionality. """

#     def __init__(self,
#                  exerciseDate,
#                  exerciseType,
#                  maturity_date,
#                  swaptionType,
#                  fixedCoupon,
#                  fixedFrequencyType,
#                  fixedDayCountType,
#                  notional=ONE_MILLION,
#                  floatFrequencyType=FrequencyTypes.QUARTERLY,
#                  floatDayCountType=DayCountTypes.THIRTY_E_360,
#                  calendar_type=CalendarTypes.WEEKEND,
#                  bus_day_adjust_type=BusDayAdjustTypes.FOLLOWING,
#                  date_gen_rule_type=DateGenRuleTypes.BACKWARD):
#         """ Create a Bermudan swaption contract. This is an option to enter
#         into a swap at a fixed coupon on all of the fixed leg coupon dates
#         until the exercise date. """

#         if exerciseDate > maturity_date:
#             raise FinError("Exercise date must be before swap maturity date")

#         if exerciseType not in FinSwaptionExerciseTypes:
#             raise FinError("Exercise type must be a FinSwaptionExerciseTypes")

#         if fixedDayCountType not in DayCountTypes:
#             raise FinError(
#                 "Unknown Fixed Day Count Rule type " +
#                 str(fixedDayCountType))

#         if fixedFrequencyType not in FrequencyTypes:
#             raise FinError(
#                 "Unknown Fixed Frequency type " +
#                 str(fixedFrequencyType))

#         if floatDayCountType not in DayCountTypes:
#             raise FinError(
#                 "Unknown Float Day Count Rule type " +
#                 str(floatDayCountType))

#         if floatFrequencyType not in FrequencyTypes:
#             raise FinError(
#                 "Unknown Float Frequency type " +
#                 str(fixedFrequencyType))

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

#         self._exerciseDate = exerciseDate
#         self._maturity_date = maturity_date
#         self._fixedCoupon = fixedCoupon
#         self._fixedFrequencyType = fixedFrequencyType
#         self._fixedDayCountType = fixedDayCountType
#         self._notional = notional
#         self._floatFrequencyType = floatFrequencyType
#         self._floatDayCountType = floatDayCountType

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

#         floatSpread = 0.0
#         payFixedFlag = True

#         # The underlying is a swap in which we pay the fixed amount
#         swap = FinIborSwap(self._exerciseDate,
#                             self._maturity_date,
#                             self._fixedCoupon,
#                             self._fixedFrequencyType,
#                             self._fixedDayCountType,
#                             self._notional,
#                             floatSpread,
#                             self._floatFrequencyType,
#                             self._floatDayCountType,
#                             payFixedFlag,
#                             self._calendar_type,
#                             self._bus_day_adjust_type,
#                             self._date_gen_rule_type)

#         swap.generateFlows()
#         cpnTimes = []
#         cpnAmounts = []

#         for iFlow in range(1, len(self._swap._adjustedFixedDates)):
#             flowDate= swap._adjustedFixedDates[iFlow]
#             cpnTime = (flowDate - settlement_date) / gDaysInYear
#             cpnFlow = swap._fixedFlows[iFlow-1] / self._notional
#             cpnTimes.append(cpnTime)
#             cpnAmounts.append(cpnFlow)

#         cpnTimes = np.array(cpnTimes)
#         cpnAmounts = np.array(cpnAmounts)

#         # Generate bond call times and prices
#         call_times = []
#         for dt in self._call_dates:
#             call_time = (dt - settlement_date) / gDaysInYear
#             call_times.append(call_time)
#         call_times = np.array(call_times)
#         call_prices = np.array(self._call_prices)

#         # Generate bond put times and prices
#         if self._swaptionType == FinIborSwaptionType.PAY:
#             call_price = 100.0
#             putPrice = 1e10
#         else:
#             call_price = 1e10
#             putPrice = 100.0

#         put_times = []
#         for putDate in swap._adjustedFixedDates[1:]:
#             if putDate <= self._exerciseDate: 
#                 put_time = (putDate - settlement_date) / gDaysInYear
#                 put_times.append(put_time)

#         put_times = np.array(put_times)
#         put_prices = np.array(self._put_prices)

#         maturity_date = self._bond._maturity_date
#         tmat = (maturity_date - settlement_date) / gDaysInYear
#         dfTimes = discount_curve._times
#         df_values = discount_curve._values

#         face = self._bond._face

#         if type(model) is FinModelRatesHW:

#             """ We need to build the tree out to the bond maturity date. To be
#             more precise we only need to go out the the last option date but
#             we can do that refinement at a later date. """

#             model.buildTree(tmat, dfTimes, df_values)
#             v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices, face)
#             model._numTimeSteps += 1
#             model.buildTree(tmat, dfTimes, df_values)
#             v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices, face)
#             model._numTimeSteps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

#         elif type(model) == FinModelRatesBK:

#             """ Because we not have a closed form bond price we need to build
#             the tree out to the bond maturity which is after option expiry. """

#             model.buildTree(tmat, dfTimes, df_values)
#             v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices,
#                                                  face)
#             model._numTimeSteps += 1
#             model.buildTree(tmat, dfTimes, df_values)
#             v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  call_times, call_prices,
#                                                  put_times, put_prices,
#                                                  face)
#             model._numTimeSteps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
#         else:
#             raise FinError("Unknown model type")

# ###############################################################################

#     def __repr__(self):

#         s = labelToString("MATURITY DATE", self._maturity_date)
#         s += labelToString("EXERCISE DATE", self._exerciseDate)
#         s += labelToString("COUPON", self._coupon)
#         s += labelToString("FREQUENCY", self._freq_type)
#         s += labelToString("ACCRUAL TYPE", self._accrual_type)
#         s += labelToString("FACE AMOUNT", self._face)
#         s += labelToString("CONVERSION RATIO", self._conversionRatio)
#         s += labelToString("START CONVERT DATE", self._startConvertDate)

#         for i in range(0, len(self._call_dates)):
#             s += labelToString("CALL DATE AND PRICE", self._call_dates[i],
#                                self._call_prices[i])

#         for i in range(0, len(self._put_dates)):
#             s += labelToString("PUT DATE AND PRICE", self._put_dates[i],
#                                self._put_prices[i])

#         return s

# ###############################################################################

#     def print(self):
#         print(self)

# ###############################################################################
