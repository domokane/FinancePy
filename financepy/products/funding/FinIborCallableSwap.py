# ##############################################################################
# # Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# ##############################################################################

# from ...finutils.FinGlobalVariables import gDaysInYear
# from ...models.FinModelRatesHW import FinModelRatesHW
# from ...models.FinModelRatesBK import FinModelRatesBK
# from ...finutils.FinError import FinError
# from ...finutils.FinFrequency import FinFrequencyTypes
# from ...finutils.FinDayCount import FinDayCount
# from ...finutils.FinDayCount import FinDayCountTypes
# from ...products.bonds.FinBond import FinBond


# from ...finutils.FinDate import FinDate
# from ...finutils.FinCalendar import FinCalendar, FinCalendarTypes
# from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
# from ...finutils.FinSchedule import FinSchedule
# from ...finutils.FinMath import ONE_MILLION

# from ...finutils.FinHelperFunctions import labelToString

# from enum import Enum
# import numpy as np

# ###############################################################################


# class FinSwaptionModelTypes(Enum):
#     HULL_WHITE = 1
#     BLACK_KARASINSKI = 2

# ###############################################################################


# class FinBondOptionTypes(Enum):
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
#     ''' Class for fixed coupon bonds with embedded call or put optionality. '''

#     def __init__(self,
#                  exerciseDate,
#                  exerciseType,
#                  maturityDate,
#                  swaptionType,
#                  fixedCoupon,
#                  fixedFrequencyType,
#                  fixedDayCountType,
#                  notional=ONE_MILLION,
#                  floatFrequencyType=FinFrequencyTypes.QUARTERLY,
#                  floatDayCountType=FinDayCountTypes.THIRTY_E_360,
#                  calendarType=FinCalendarTypes.WEEKEND,
#                  busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
#                  dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
#         ''' Create a Bermudan swaption contract. This is an option to enter
#         into a swap at a fixed coupon on all of the fixed leg coupon dates
#         until the exercise date. '''

#         if exerciseDate > maturityDate:
#             raise FinError("Exercise date must be before swap maturity date")

#         if exerciseType not in FinSwaptionExerciseTypes:
#             raise FinError("Exercise type must be a FinSwaptionExerciseTypes")

#         if fixedDayCountType not in FinDayCountTypes:
#             raise FinError(
#                 "Unknown Fixed Day Count Rule type " +
#                 str(fixedDayCountType))

#         if fixedFrequencyType not in FinFrequencyTypes:
#             raise FinError(
#                 "Unknown Fixed Frequency type " +
#                 str(fixedFrequencyType))

#         if floatDayCountType not in FinDayCountTypes:
#             raise FinError(
#                 "Unknown Float Day Count Rule type " +
#                 str(floatDayCountType))

#         if floatFrequencyType not in FinFrequencyTypes:
#             raise FinError(
#                 "Unknown Float Frequency type " +
#                 str(fixedFrequencyType))

#         if calendarType not in FinCalendarTypes:
#             raise FinError("Unknown Calendar type " + str(calendarType))

#         if busDayAdjustType not in FinBusDayAdjustTypes:
#             raise FinError(
#                 "Unknown Business Day Adjust type " +
#                 str(busDayAdjustType))

#         if dateGenRuleType not in FinDateGenRuleTypes:
#             raise FinError(
#                 "Unknown Date Gen Rule type " +
#                 str(dateGenRuleType))

#         self._exerciseDate = exerciseDate
#         self._maturityDate = maturityDate
#         self._fixedCoupon = fixedCoupon
#         self._fixedFrequencyType = fixedFrequencyType
#         self._fixedDayCountType = fixedDayCountType
#         self._notional = notional
#         self._floatFrequencyType = floatFrequencyType
#         self._floatDayCountType = floatDayCountType

#         self._calendarType = calendarType
#         self._busDayAdjustType = busDayAdjustType
#         self._dateGenRuleType = dateGenRuleType

#         self._pv01 = None
#         self._fwdSwapRate = None
#         self._forwardDf = None
#         self._underlyingSwap = None

# ###############################################################################

#     def value(self,
#               valuationDate,
#               discountCurve,
#               model):
#         ''' Value the bermuda swaption. This is done using the specified
#         model and a discount curve. '''

#         floatSpread = 0.0
#         payFixedFlag = True

#         # The underlying is a swap in which we pay the fixed amount
#         swap = FinIborSwap(self._exerciseDate,
#                             self._maturityDate,
#                             self._fixedCoupon,
#                             self._fixedFrequencyType,
#                             self._fixedDayCountType,
#                             self._notional,
#                             floatSpread,
#                             self._floatFrequencyType,
#                             self._floatDayCountType,
#                             payFixedFlag,
#                             self._calendarType,
#                             self._busDayAdjustType,
#                             self._dateGenRuleType)

#         swap.generateFlows()
#         cpnTimes = []
#         cpnAmounts = []

#         for iFlow in range(1, len(self._swap._adjustedFixedDates)):
#             flowDate= swap._adjustedFixedDates[iFlow]
#             cpnTime = (flowDate - settlementDate) / gDaysInYear
#             cpnFlow = swap._fixedFlows[iFlow-1] / self._notional
#             cpnTimes.append(cpnTime)
#             cpnAmounts.append(cpnFlow)

#         cpnTimes = np.array(cpnTimes)
#         cpnAmounts = np.array(cpnAmounts)

#         # Generate bond call times and prices
#         callTimes = []
#         for dt in self._callDates:
#             callTime = (dt - settlementDate) / gDaysInYear
#             callTimes.append(callTime)
#         callTimes = np.array(callTimes)
#         callPrices = np.array(self._callPrices)

#         # Generate bond put times and prices
#         if self._swaptionType == FinIborSwaptionType.PAY:
#             callPrice = 100.0
#             putPrice = 1e10
#         else:
#             callPrice = 1e10
#             putPrice = 100.0

#         putTimes = []
#         for putDate in swap._adjustedFixedDates[1:]:
#             if putDate <= self._exerciseDate: 
#                 putTime = (putDate - settlementDate) / gDaysInYear
#                 putTimes.append(putTime)

#         putTimes = np.array(putTimes)
#         putPrices = np.array(self._putPrices)

#         maturityDate = self._bond._maturityDate
#         tmat = (maturityDate - settlementDate) / gDaysInYear
#         dfTimes = discountCurve._times
#         dfValues = discountCurve._values

#         face = self._bond._face

#         if type(model) is FinModelRatesHW:

#             ''' We need to build the tree out to the bond maturity date. To be
#             more precise we only need to go out the the last option date but
#             we can do that refinement at a later date. '''

#             model.buildTree(tmat, dfTimes, dfValues)
#             v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  callTimes, callPrices,
#                                                  putTimes, putPrices, face)
#             model._numTimeSteps += 1
#             model.buildTree(tmat, dfTimes, dfValues)
#             v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  callTimes, callPrices,
#                                                  putTimes, putPrices, face)
#             model._numTimeSteps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

#         elif type(model) == FinModelRatesBK:

#             ''' Because we not have a closed form bond price we need to build
#             the tree out to the bond maturity which is after option expiry. '''

#             model.buildTree(tmat, dfTimes, dfValues)
#             v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  callTimes, callPrices,
#                                                  putTimes, putPrices,
#                                                  face)
#             model._numTimeSteps += 1
#             model.buildTree(tmat, dfTimes, dfValues)
#             v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
#                                                  callTimes, callPrices,
#                                                  putTimes, putPrices,
#                                                  face)
#             model._numTimeSteps -= 1

#             v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
#             v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

#             return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
#         else:
#             raise FinError("Unknown model type")

# ###############################################################################

#     def __repr__(self):

#         s = labelToString("MATURITY DATE", self._maturityDate)
#         s += labelToString("EXERCISE DATE", self._exerciseDate)
#         s += labelToString("COUPON", self._coupon)
#         s += labelToString("FREQUENCY", self._freqType)
#         s += labelToString("ACCRUAL TYPE", self._accrualType)
#         s += labelToString("FACE AMOUNT", self._face)
#         s += labelToString("CONVERSION RATIO", self._conversionRatio)
#         s += labelToString("START CONVERT DATE", self._startConvertDate)

#         for i in range(0, len(self._callDates)):
#             s += labelToString("CALL DATE AND PRICE", self._callDates[i],
#                                self._callPrices[i])

#         for i in range(0, len(self._putDates)):
#             s += labelToString("PUT DATE AND PRICE", self._putDates[i],
#                                self._putPrices[i])

#         return s

# ###############################################################################

#     def print(self):
#         print(self)

# ###############################################################################
