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

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate

from ...products.rates.FinIborSwap import FinIborSwap

from ...models.FinModelBlack import FinModelBlack
from ...models.FinModelBlackShifted import FinModelBlackShifted
from ...models.FinModelSABR import FinModelSABR
from ...models.FinModelSABRShifted import FinModelSABRShifted
from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...models.FinModelRatesBDT import FinModelRatesBDT

from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinGlobalTypes import FinSwapTypes
from ...finutils.FinGlobalTypes import FinExerciseTypes

###############################################################################


class FinIborSwaption():
    ''' This is the class for the European-style swaption, an option to enter
    into a swap (payer or receiver of the fixed coupon), that starts in the
    future and with a fixed maturity, at a swap rate fixed today. '''

    def __init__(self,
                 settlementDate: FinDate,
                 exerciseDate: FinDate,
                 maturityDate: FinDate,
                 fixedLegType: FinSwapTypes,
                 fixedCoupon: float,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 notional: float = ONE_MILLION,
                 floatFrequencyType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed coupon and the details of the fixed and the
        floating leg payment schedules. Bermudan style swaption should be
        priced using the FinIborBermudanSwaption class. '''

        checkArgumentTypes(self.__init__, locals())

        if settlementDate > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > maturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        self._settlementDate = settlementDate
        self._exerciseDate = exerciseDate
        self._maturityDate = maturityDate
        self._fixedLegType = fixedLegType

        self._notional = notional

        self._fixedCoupon = fixedCoupon
        self._fixedFrequencyType = fixedFrequencyType
        self._fixedDayCountType = fixedDayCountType
        self._floatFrequencyType = floatFrequencyType
        self._floatDayCountType = floatDayCountType

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None

###############################################################################

    def value(self,
              valuationDate,
              discountCurve,
              model):
        ''' Valuation of a Ibor European-style swaption using a choice of
        models on a specified valuation date. Models include FinModelBlack,
        FinModelBlackShifted, FinModelSABR, FinModelSABRShifted, FinModelHW,
        FinModelBK and FinModelBDT. The last two involved a tree-based
        valuation. '''

        floatSpread = 0.0

        # We create a swap that starts on the exercise date.
        swap = FinIborSwap(self._exerciseDate,
                           self._maturityDate,
                           self._fixedLegType,
                           self._fixedCoupon,
                           self._fixedFrequencyType,
                           self._fixedDayCountType,
                           self._notional,
                           floatSpread,
                           self._floatFrequencyType,
                           self._floatDayCountType,
                           self._calendarType,
                           self._busDayAdjustType,
                           self._dateGenRuleType)

        k = self._fixedCoupon

        # The pv01 is the value of the swap cashflows as of the curve date
        pv01 = swap.pv01(valuationDate, discountCurve)

        # We need to calculate the forward swap rate on the swaption exercise
        # date that makes the forward swap worth par including principal
        s = swap.swapRate(valuationDate, discountCurve)

        texp = (self._exerciseDate - self._settlementDate) / gDaysInYear
        tmat = (self._maturityDate - self._settlementDate) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        #######################################################################
        # For the tree models we need to generate a vector of the coupons
        #######################################################################

        cpnTimes = [texp]
        cpnFlows = [0.0]

        # The first flow is on the day after the expiry date
        numFlows = len(swap._fixedLeg._paymentDates)

        for iFlow in range(0, numFlows):
            
            flowDate = swap._fixedLeg._paymentDates[iFlow]

            # Only flows occurring after option expiry are counted. 
            # Flows on the expiry date are not included
            if flowDate > self._exerciseDate:
                cpnTime = (flowDate - valuationDate) / gDaysInYear
                cpnFlow = swap._fixedLeg._payments[iFlow] / self._notional
                cpnTimes.append(cpnTime)
                cpnFlows.append(cpnFlow)

        cpnTimes = np.array(cpnTimes)
        cpnFlows = np.array(cpnFlows)

        dfTimes = discountCurve._times
        dfValues = discountCurve._dfs

        if np.any(cpnTimes < 0.0):
            raise FinError("No coupon times can be before the value date.")

        strikePrice = 1.0
        faceAmount = 1.0

        #######################################################################

        if isinstance(model, FinModelBlack):

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBlackShifted):

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABR):

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABRShifted):

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelRatesHW):

            swaptionPx = model.europeanBondOptionJamshidian(texp,
                                                  strikePrice,
                                                  faceAmount, 
                                                  cpnTimes,
                                                  cpnFlows,
                                                  dfTimes,
                                                  dfValues)

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['put']
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['call']
            else:
                raise FinError("Unknown swaption option type" +
                               str(self._swapType))

            # Cancel the multiplication at the end below
            swaptionPrice /= pv01

        elif isinstance(model, FinModelRatesBK):

            model.buildTree(tmat, dfTimes, dfValues)
            swaptionPx = model.bermudanSwaption(texp,
                                                tmat,
                                                strikePrice,
                                                faceAmount,
                                                cpnTimes,
                                                cpnFlows,
                                                FinExerciseTypes.EUROPEAN)

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['pay']
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['rec']

            swaptionPrice /= pv01

        elif isinstance(model, FinModelRatesBDT):

            model.buildTree(tmat, dfTimes, dfValues)
            swaptionPx = model.bermudanSwaption(texp,
                                                tmat,
                                                strikePrice,
                                                faceAmount,
                                                cpnTimes,
                                                cpnFlows,
                                                FinExerciseTypes.EUROPEAN)

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = swaptionPx['pay']
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = swaptionPx['rec']

            swaptionPrice /= pv01
        else:
            raise FinError("Unknown swaption model " + str(model))

        self._pv01 = pv01
        self._fwdSwapRate = s
        self._forwardDf = discountCurve.df(self._exerciseDate)
        self._underlyingSwap = swap

        # The exchange of cash occurs on the settlement date. However the
        # actual value is that on the specified valuation date which could
        # be the swaption settlement date.
        dfSettlement = discountCurve.df(self._settlementDate)
        swaptionPrice = swaptionPrice * pv01 * self._notional / dfSettlement
        return swaptionPrice

###############################################################################

    def cashSettledValue(self,
                         valuationDate: FinDate,
                         discountCurve,
                         swapRate: float,
                         model):
        ''' Valuation of a Ibor European-style swaption using a cash settled
        approach which is a market convention that used Black's model and that
        discounts all of the future payments at a flat swap rate. Note that the
        Black volatility for this valuation should in general not equal the
        Black volatility for the standard arbitrage-free valuation. '''

        floatSpread = 0.0

        swap = FinIborSwap(self._exerciseDate,
                            self._maturityDate,
                            self._fixedLegType,
                            self._fixedCoupon,
                            self._fixedFrequencyType,
                            self._fixedDayCountType,
                            self._notional,
                            floatSpread,
                            self._floatFrequencyType,
                            self._floatDayCountType,
                            self._calendarType,
                            self._busDayAdjustType,
                            self._dateGenRuleType)

        k = self._fixedCoupon
        s = swapRate

        pv01 = swap.cashSettledPV01(valuationDate,
                                    swapRate,
                                    self._fixedFrequencyType)

        texp = (self._exerciseDate - self._settlementDate) / gDaysInYear

        # Discounting is done via the PV01 annuity so no discounting in Black
        df = 1.0

        if isinstance(model, FinModelBlack):

            if self._fixedLegType == FinSwapTypes.PAY:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_CALL)
            elif self._fixedLegType == FinSwapTypes.RECEIVE:
                swaptionPrice = model.value(s, k, texp, df,
                                            FinOptionTypes.EUROPEAN_PUT)
        else:
            raise FinError("Cash settled swaptions must be priced using"
                           + " Black's model.")

        self._fwdSwapRate = swapRate
        self._forwardDf = discountCurve.df(self._exerciseDate)
        self._underlyingSwap = swap
        # The annuity needs to be discounted to today using the correct df
        self._pv01 = pv01 * self._forwardDf

        # The exchange of cash occurs on the settlement date but we need to 
        # value the swaption on the provided valuation date - which could be
        # the settlement date or may be a different date.
        dfValuation = discountCurve.df(valuationDate)
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
        ''' Function to allow us to print the swaption details. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("SETTLEMENT DATE", self._settlementDate)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("SWAP FIXED LEG TYPE", str(self._fixedLegType))
        s += labelToString("SWAP MATURITY DATE", self._maturityDate)
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
        ''' Alternative print method. '''

        print(self)

###############################################################################
