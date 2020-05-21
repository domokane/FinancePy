##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import log, sqrt

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION, N
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString
from ...products.libor.FinLiborSwap import FinLiborSwap
from ...models.FinModelSABR import blackVolFromSABR

from .FinLiborModelTypes import FinLiborModelBlack
from .FinLiborModelTypes import FinLiborModelSABR

##########################################################################

from enum import Enum


class FinLiborSwaptionType(Enum):
    PAYER = 1
    RECEIVER = 2


class FinLiborSwaptionModelTypes(Enum):
    BLACK = 1
    SABR = 2

##########################################################################


class FinLiborSwaption():

    def __init__(self,
                 exerciseDate,
                 swapMaturityDate,
                 swaptionType,
                 swapFixedCoupon,
                 swapFixedFrequencyType,
                 swapFixedDayCountType,
                 swapNotional=ONE_MILLION,
                 swapFloatFrequencyType=FinFrequencyTypes.QUARTERLY,
                 swapFloatDayCountType=FinDayCountTypes.THIRTY_360,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):

        if exerciseDate > swapMaturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        if swapFixedDayCountType not in FinDayCountTypes:
            raise FinError(
                "Unknown Fixed DayCountRule type " +
                str(swapFixedDayCountType))

        if swapFixedFrequencyType not in FinFrequencyTypes:
            raise FinError(
                "Unknown Fixed Frequency type " +
                str(swapFixedFrequencyType))

        if calendarType not in FinCalendarTypes:
            raise FinError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayAdjustTypes:
            raise FinError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise FinError(
                "Unknown Date Gen Rule type " +
                str(dateGenRuleType))

        self._exerciseDate = exerciseDate
        self._maturityDate = swapMaturityDate
        self._swaptionType = swaptionType
        self._swapFixedCoupon = swapFixedCoupon
        self._swapFixedFrequencyType = swapFixedFrequencyType
        self._swapFixedDayCountType = swapFixedDayCountType
        self._swapNotional = swapNotional
        self._swapFloatFrequencyType = swapFloatFrequencyType
        self._swapFloatDayCountType = swapFloatDayCountType

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._pv01 = None
        self._fwdSwapRate = None
        self._forwardDf = None
        self._underlyingSwap = None

##########################################################################

    def value(self,
              valuationDate,
              liborCurve,
              model):

        floatSpread = 0.0
        payFixedFlag = True

        swap = FinLiborSwap(self._exerciseDate,
                            self._maturityDate,
                            self._swapFixedCoupon,
                            self._swapFixedFrequencyType,
                            self._swapFixedDayCountType,
                            self._swapNotional,
                            floatSpread,
                            self._swapFloatFrequencyType,
                            self._swapFloatDayCountType,
                            payFixedFlag,
                            self._calendarType,
                            self._busDayAdjustType,
                            self._dateGenRuleType)

        k = self._swapFixedCoupon
        pv01 = swap.pv01(valuationDate, liborCurve)
        f = swap.parCoupon(valuationDate, liborCurve)
        t = (self._exerciseDate - liborCurve._curveDate) / gDaysInYear
        forwardDf = liborCurve.df(t)

        if type(model) == FinLiborModelBlack:

            v = model._volatility
            d1 = (log(f / k) + v * v * t / 2.0) / v / sqrt(t)
            d2 = d1 - v * sqrt(t)

            if self._swaptionType == FinLiborSwaptionType.PAYER:
                swaptionPrice = f * N(d1) - k * N(d2)
            elif self._swaptionType == FinLiborSwaptionType.RECEIVER:
                swaptionPrice = k * N(-d2) - f * N(-d1)
            else:
                raise FinError("Unknown swaption option type" +
                               str(self._optionType))

        elif type(model) == FinLiborModelSABR:

            alpha = model._alpha
            beta = model._beta
            rho = model._rho
            nu = model._nu

            v = blackVolFromSABR(alpha, beta, rho, nu, f, k, t)
            d1 = (log(f / k) + v * v * t / 2.0) / v / sqrt(t)
            d2 = d1 - v * sqrt(t)

            if self._swaptionType == FinLiborSwaptionType.PAYER:
                swaptionPrice = f * N(d1) - k * N(d2)
            elif self._swaptionType == FinLiborSwaptionType.RECEIVER:
                swaptionPrice = k * N(-d2) - f * N(-d1)
            else:
                raise FinError("Unknown swaption option type" +
                               str(self._optionType))

        else:
            raise FinError("Unknown swaption model " + str(model))

        self._pv01 = pv01
        self._fwdSwapRate = f
        self._forwardDf = forwardDf
        self._underlyingSwap = swap

        swaptionPrice = swaptionPrice * pv01 * self._swapNotional
        return swaptionPrice

###############################################################################

    def printFixedLeg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.printFixedLeg()

###############################################################################

    def printFloatLeg(self):

        if self._underlyingSwap is None:
            raise FinError("Underlying swap has not been set. Do a valuation.")

        self._underlyingSwap.printFloatLeg()

###############################################################################

    def __repr__(self):

        s = labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("SWAPTION TYPE", str(self._swaptionType))
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("SWAP NOTIONAL", self._swapNotional)
        s += labelToString("FIXED COUPON", self._swapFixedCoupon * 100)
        s += labelToString("FIXED FREQUENCY", str(self._swapFixedFrequencyType))
        s += labelToString("FIXED DAY COUNT", str(self._swapFixedDayCountType))
        s += labelToString("FLOAT FREQUENCY", str(self._swapFixedFrequencyType))
        s += labelToString("FLOAT DAY COUNT", str(self._swapFloatDayCountType))

        if self._pv01 is not None:
            s += labelToString("PV01", self._pv01)
            s += labelToString("FWD SWAP RATE", self._fwdSwapRate*100)
            s += labelToString("FWD DF TO EXPIRY", self._forwardDf, "")

        return s

###############################################################################

    def print(self):
        print(self)

###############################################################################
