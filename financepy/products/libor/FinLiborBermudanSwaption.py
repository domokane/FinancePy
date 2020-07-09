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
#     PAYER = 1
#     RECEIVER = 2

# class FinSwaptionExerciseTypes(Enum):
#     EUROPEAN = 1
#     BERMUDAN = 2

# ###############################################################################


# class FinLiborBermudanSwaption(object):
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
#                  floatDayCountType=FinDayCountTypes.THIRTY_360,
#                  calendarType=FinCalendarTypes.WEEKEND,
#                  busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
#                  dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
#         ''' Create a Bermudan swaption contract. This is an option to enter
#         into a payer or receiver swap at a fixed coupon on all of the fixed
#         # leg coupon dates until the exercise date inclusive. '''

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
#         swap = FinLiborSwap(self._exerciseDate,
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
#             cpnFlow = swap._fixedFlows[iFlow] / self._notional
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
#         if self._swaptionType == FinLiborSwaptionType.PAYER:
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

#         if isinstance(model, FinModelRatesHW):

#             if self._optionType == FinBondOptionTypes.EUROPEAN_CALL \
#                     and model._useJamshidian is True:

#                 v = model.europeanBondOption_Jamshidian(texp,
#                                                         self._strikePrice,
#                                                         self._face,
#                                                         cpnTimes,
#                                                         cpnAmounts,
#                                                         dfTimes, dfValues)

#                 return v['call']

#             elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT  \
#                     and model._useJamshidian is True:

#                 v = model.europeanBondOption_Jamshidian(texp,
#                                                         self._strikePrice,
#                                                         self._face,
#                                                         cpnTimes,
#                                                         cpnAmounts,
#                                                         dfTimes, dfValues)

#                 return v['put']

#             elif self._optionType == FinBondOptionTypes.EUROPEAN_CALL  \
#                     and model._useJamshidian is False:

#                 model.buildTree(texp, dfTimes, dfValues)
#                 americanExercise = False

#                 v1 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps += 1
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v2 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 v = (v1['call'] + v2['call'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT  \
#                     and model._useJamshidian is False:

#                 americanExercise = False
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v1 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps += 1
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v2 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps -= 1
#                 v = (v1['put'] + v2['put'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

#                 americanExercise = True
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v1 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps += 1
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v2 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps -= 1
#                 v = (v1['call'] + v2['call'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

#                 americanExercise = True
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v1 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps += 1
#                 model.buildTree(texp, dfTimes, dfValues)
#                 v2 = model.americanBondOption_Tree(texp, self._strikePrice,
#                                                    self._face,
#                                                    cpnTimes, cpnAmounts,
#                                                    americanExercise)

#                 model._numTimeSteps -= 1
#                 v = (v1['put'] + v2['put'])/2.0
#                 return v

#         elif type(model) == FinModelRatesBK:

#             maturityDate = self._bond._maturityDate
#             tmat = (maturityDate - valueDate) / gDaysInYear

#             if self._optionType == FinBondOptionTypes.EUROPEAN_CALL:

#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v1 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, False)
#                 model._numTimeSteps += 1
#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v2 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)

#                 v = (v1['call'] + v2['call'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT:

#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v1 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, False)

#                 model._numTimeSteps += 1
#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v2 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)

#                 v = (v1['put'] + v2['put'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v1 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)

#                 model._numTimeSteps += 1
#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v2 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)

#                 v = (v1['call'] + v2['call'])/2.0
#                 return v

#             elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v1 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)
#                 model._numTimeSteps += 1
#                 model.buildTree(tmat, dfTimes, dfValues)
#                 v2 = model.bondOption(texp, self._strikePrice, self._face,
#                                       cpnTimes, cpnAmounts, True)

#                 v = (v1['put'] + v2['put'])/2.0
#                 return v

#         else:
#             raise FinError("Unknown model and option combination")


# ###############################################################################

#     def __repr__(self):

#         s = labelToString("MATURITY DATE", self._maturityDate)
#         s += labelToString("EXERCISE DATE", self._exerciseDate)
#         s += labelToString("COUPON", self._coupon)
#         s += labelToString("FREQUENCY", self._frequencyType)
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
