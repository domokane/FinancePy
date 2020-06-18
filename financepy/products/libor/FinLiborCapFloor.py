##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Implied volatility
# TODO: Term structure of volatility
# TODO: Check that curve anchor date is valuation date ?

import numpy as np
from typing import Union, Optional

from ...finutils.FinDate import FinDate
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinError import FinError
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...models.FinModelBlack import FinModelBlack
from ...models.FinModelBlackShifted import FinModelBlackShifted
from ...models.FinModelBachelier import FinModelBachelier
from ...models.FinModelSABR import FinModelSABR
from ...models.FinModelSABRShifted import FinModelSABRShifted
from ...models.FinModelRatesHW import FinModelRatesHW

from ...finutils.FinOptionTypes import FinOptionTypes

##########################################################################

from enum import Enum


class FinLiborCapFloorTypes(Enum):
    CAP = 1
    FLOOR = 2


class FinLiborCapFloorModelTypes(Enum):
    BLACK = 1
    SHIFTED_BLACK = 2
    SABR = 3

##########################################################################


class FinLiborCapFloor():

    def __init__(self,
                 startDate: FinDate,
                 maturityDateOrTenor: Union[FinDate, str],
                 optionType: FinLiborCapFloorType,
                 strikeRate: float,
                 lastFixing: Optional[float] = None,
                 frequencyType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360_ISDA,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):

        checkArgumentTypes(self.__init__, locals())

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = startDate.addTenor(maturityDateOrTenor)
            calendar = FinCalendar(self._calendarType)
            maturityDate = calendar.adjust(maturityDate,
                                           self._busDayAdjustType)

        if startDate > maturityDate:
            raise FinError("Start date must be before maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._optionType = optionType
        self._strikeRate = strikeRate
        self._lastFixing = lastFixing
        self._frequencyType = frequencyType
        self._dayCountType = dayCountType
        self._notional = notional
        self._dateGenRuleType = dateGenRuleType

        self._capFloorLetValues = []
        self._capFloorLetAlphas = []
        self._capFloorLetFwdRates = []
        self._capFloorLetIntrinsic = []
        self._capFloorLetDiscountFactors = []
        self._capFloorPV = []

        self._valuationDate = None
        self._dayCounter = None

##########################################################################

    def value(self,
              valuationDate,
              liborCurve,
              model):

        self._valuationDate = valuationDate

        self._capFloorDates = FinSchedule(self._startDate,
                                          self._maturityDate,
                                          self._frequencyType,
                                          self._calendarType,
                                          self._busDayAdjustType,
                                          self._dateGenRuleType).generate()

        self._dayCounter = FinDayCount(self._dayCountType)
        numOptions = len(self._capFloorDates)
        strikeRate = self._strikeRate

        if strikeRate < 0.0:
            raise FinError("Strike < 0.0")

        if numOptions <= 1:
            raise FinError("Number of options in capfloor equals 1")

        #######################################################################

        self._capFloorLetValues = [0]
        self._capFloorLetAlphas = [0]
        self._capFloorLetFwdRates = [0]
        self._capFloorLetIntrinsic = [0]
        self._capFloorLetDiscountFactors = [1.00]
        self._capFloorPV = [0.0]

        #######################################################################

        capFloorValue = 0.0
        capFloorLetValue = 0.0
        # Value the first caplet or floorlet with known payoff

        startDate = self._startDate
        endDate = self._capFloorDates[1]

        if self._lastFixing is None:
            fwdRate = liborCurve.fwdRate(startDate, endDate,
                                         self._dayCountType)
        else:
            fwdRate = self._lastFixing

        alpha = self._dayCounter.yearFrac(startDate, endDate)
        df = liborCurve.df(endDate)

        if self._optionType == FinLiborCapFloorTypes.CAP:
            capFloorLetValue = df * alpha * max(fwdRate - strikeRate, 0)
        elif self._optionType == FinLiborCapFloorTypes.FLOOR:
            capFloorLetValue = df * alpha * max(strikeRate - fwdRate, 0)

        capFloorLetValue *= self._notional
        capFloorValue += capFloorLetValue

        self._capFloorLetFwdRates.append(fwdRate)
        self._capFloorLetValues.append(capFloorLetValue)
        self._capFloorLetAlphas.append(alpha)
        self._capFloorLetIntrinsic.append(capFloorLetValue)
        self._capFloorLetDiscountFactors.append(df)
        self._capFloorPV.append(capFloorValue)

        for i in range(2, numOptions):

            startDate = self._capFloorDates[i - 1]
            endDate = self._capFloorDates[i]
            alpha = self._dayCounter.yearFrac(startDate, endDate)
            df = liborCurve.df(endDate)
            fwdRate = liborCurve.fwdRate(startDate, endDate,
                                         self._dayCountType)

            if self._optionType == FinLiborCapFloorTypes.CAP:
                intrinsicValue = df * alpha * max(fwdRate - strikeRate, 0)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                intrinsicValue = df * alpha * max(strikeRate - fwdRate, 0)

            intrinsicValue *= self._notional

            capFloorLetValue = self.valueCapletFloorlet(valuationDate,
                                                        startDate,
                                                        endDate,
                                                        liborCurve,
                                                        model)

            capFloorLetValue *= self._notional * alpha
            capFloorValue += capFloorLetValue

            self._capFloorLetFwdRates.append(fwdRate)
            self._capFloorLetValues.append(capFloorLetValue)
            self._capFloorLetAlphas.append(alpha)
            self._capFloorLetIntrinsic.append(intrinsicValue)
            self._capFloorLetDiscountFactors.append(df)
            self._capFloorPV.append(capFloorValue)

        return capFloorValue

##########################################################################

    def valueCapletFloorlet(self,
                            valuationDate,
                            startDate,
                            endDate,
                            liborCurve,
                            model):

        # Uncertainty starts after the swaption settles on its start date
        texp = (startDate - self._startDate) / gDaysInYear
        tend = (endDate - self._startDate) / gDaysInYear

        f = liborCurve.fwdRate(startDate, endDate, self._dayCountType)
        k = self._strikeRate
        df = liborCurve.df(tend)

        if k == 0.0:
            k = 1e-10

        if isinstance(model, FinModelBlack):

            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBlackShifted):

            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBachelier):

            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABR):

            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABRShifted):

            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelRatesHW):

            tmat = (endDate - valuationDate) / gDaysInYear
            alpha = self._dayCounter.yearFrac(startDate, endDate)
            strikePrice = 1.0/(1.0 + alpha * self._strikeRate)
            notionalAdj = (1.0 + self._strikeRate * alpha)
            face = 1.0
            dfTimes = liborCurve._times
            dfValues = liborCurve._values

            v = model.optionOnZCB(texp, tmat, strikePrice, face,
                                  dfTimes, dfValues)

            # we divide by alpha to offset the multiplication above
            if self._optionType == FinLiborCapFloorTypes.CAP:
                capFloorLetValue = v['put'] * notionalAdj / alpha
            elif self._optionType == FinLiborCapFloorTypes.FLOOR:
                capFloorLetValue = v['call'] * notionalAdj / alpha

        else:
            raise FinError("Unknown model type " + str(model))

        return capFloorLetValue

###############################################################################

    def printLeg(self):
        ''' Prints the cap floor amounts. '''

        print("START DATE:", self._startDate)
        print("MATURITY DATE:", self._maturityDate)
        print("OPTION TYPE", str(self._optionType))
        print("STRIKE (%):", self._strikeRate * 100)
        print("FREQUENCY:", str(self._frequencyType))
        print("DAY COUNT:", str(self._dayCountType))
        print("VALUATION DATE", self._valuationDate)

        if len(self._capFloorLetValues) == 0:
            print("Caplets not calculated.")
            return

        if self._optionType == FinLiborCapFloorTypes.CAP:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    CAPLET_PV       CUM_PV"
        elif self._optionType == FinLiborCapFloorTypes.FLOOR:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    FLRLET_PV       CUM_PV"

        print(header)

        iFlow = 0

        for paymentDate in self._capFloorDates:
            print("%15s %10.7f  %9.5f %12.2f %12.6f %12.2f %12.2f" %
                  (paymentDate,
                   self._capFloorLetAlphas[iFlow],
                   self._capFloorLetFwdRates[iFlow]*100,
                   self._capFloorLetIntrinsic[iFlow],
                   self._capFloorLetDiscountFactors[iFlow],
                   self._capFloorLetValues[iFlow],
                   self._capFloorPV[iFlow]))

            iFlow += 1

###############################################################################

    def __repr__(self):
        s = labelToString("START DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("STRIKE COUPON", self._strikeRate * 100)
        s += labelToString("OPTION TYPE", str(self._optionType))
        s += labelToString("FREQUENCY", str(self._frequencyType))
        s += labelToString("DAY COUNT", str(self._dayCountType), "")
        return s

###############################################################################

    def print(self):
        print(self)

###############################################################################
