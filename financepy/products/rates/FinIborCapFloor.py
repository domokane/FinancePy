##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Implied volatility
# TODO: Term structure of volatility
# TODO: Check that curve anchor date is valuation date ?

from typing import Optional

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
from ...finutils.FinGlobalTypes import FinCapFloorTypes, FinOptionTypes

##########################################################################

from enum import Enum

class FinIborCapFloorModelTypes(Enum):
    BLACK = 1
    SHIFTED_BLACK = 2
    SABR = 3

##########################################################################


class FinIborCapFloor():
    ''' Class for Caps and Floors. These are contracts which observe a Ibor
    reset L on a future start date and then make a payoff at the end of the
    Ibor period which is Max[L-K,0] for a cap and Max[K-L,0] for a floor.
    This is then day count adjusted for the Ibor period and then scaled by
    the contract notional to produce a valuation. A number of models can be
    selected from.'''

    def __init__(self,
                 startDate: FinDate,
                 maturityDateOrTenor: (FinDate, str),
                 optionType: FinCapFloorTypes,
                 strikeRate: float,
                 lastFixing: Optional[float] = None,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360_ISDA,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Initialise FinIborCapFloor object. '''

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
        self._freqType = freqType
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

###############################################################################

    def _generateDates(self):

        schedule = FinSchedule(self._startDate,
                               self._maturityDate,
                               self._freqType,
                               self._calendarType,
                               self._busDayAdjustType,
                               self._dateGenRuleType)

        self._capFloorLetDates = schedule._adjustedDates

##########################################################################

    def value(self, valuationDate, liborCurve, model):
        ''' Value the cap or floor using the chosen model which specifies
        the volatility of the Ibor rate to the cap start date. '''

        self._valuationDate = valuationDate
        self._generateDates()

        self._dayCounter = FinDayCount(self._dayCountType)
        numOptions = len(self._capFloorLetDates)
        strikeRate = self._strikeRate

        if strikeRate < 0.0:
            raise FinError("Strike < 0.0")

        if numOptions <= 1:
            raise FinError("Number of options in capfloor equals 1")

        self._capFloorLetValues = [0]
        self._capFloorLetAlphas = [0]
        self._capFloorLetFwdRates = [0]
        self._capFloorLetIntrinsic = [0]
        self._capFloorLetDiscountFactors = [1.00]
        self._capFloorPV = [0.0]

        capFloorValue = 0.0
        capFloorLetValue = 0.0
        # Value the first caplet or floorlet with known payoff

        startDate = self._startDate
        endDate = self._capFloorLetDates[1]

        if self._lastFixing is None:
            fwdRate = liborCurve.fwdRate(startDate, endDate,
                                         self._dayCountType)
        else:
            fwdRate = self._lastFixing

        alpha = self._dayCounter.yearFrac(startDate, endDate)[0]
        df = liborCurve.df(endDate)

        if self._optionType == FinCapFloorTypes.CAP:
            capFloorLetValue = df * alpha * max(fwdRate - strikeRate, 0.0)
        elif self._optionType == FinCapFloorTypes.FLOOR:
            capFloorLetValue = df * alpha * max(strikeRate - fwdRate, 0.0)

        capFloorLetValue *= self._notional
        capFloorValue += capFloorLetValue

        self._capFloorLetFwdRates.append(fwdRate)
        self._capFloorLetValues.append(capFloorLetValue)
        self._capFloorLetAlphas.append(alpha)
        self._capFloorLetIntrinsic.append(capFloorLetValue)
        self._capFloorLetDiscountFactors.append(df)
        self._capFloorPV.append(capFloorValue)

        for i in range(2, numOptions):

            startDate = self._capFloorLetDates[i - 1]
            endDate = self._capFloorLetDates[i]
            alpha = self._dayCounter.yearFrac(startDate, endDate)[0]

            df = liborCurve.df(endDate)
            fwdRate = liborCurve.fwdRate(startDate, endDate,
                                         self._dayCountType)

            if self._optionType == FinCapFloorTypes.CAP:
                intrinsicValue = df * alpha * max(fwdRate - strikeRate, 0.0)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                intrinsicValue = df * alpha * max(strikeRate - fwdRate, 0.0)

            intrinsicValue *= self._notional

            capFloorLetValue = self.valueCapletFloorLet(valuationDate,
                                                        startDate,
                                                        endDate,
                                                        liborCurve,
                                                        model)

            capFloorValue += capFloorLetValue

            self._capFloorLetFwdRates.append(fwdRate)
            self._capFloorLetValues.append(capFloorLetValue)
            self._capFloorLetAlphas.append(alpha)
            self._capFloorLetIntrinsic.append(intrinsicValue)
            self._capFloorLetDiscountFactors.append(df)
            self._capFloorPV.append(capFloorValue)

        return capFloorValue

###############################################################################

    def valueCapletFloorLet(self,
                            valuationDate,
                            capletStartDate,
                            capletEndDate,
                            liborCurve,
                            model):
        ''' Value the caplet or floorlet using a specific model. '''

        texp = (capletStartDate - self._startDate) / gDaysInYear

        alpha = self._dayCounter.yearFrac(capletStartDate, capletEndDate)[0]

        f = liborCurve.fwdRate(capletStartDate, capletEndDate,
                               self._dayCountType)

        k = self._strikeRate
        df = liborCurve.df(capletEndDate)

        if k == 0.0:
            k = 1e-10

        if isinstance(model, FinModelBlack):

            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBlackShifted):

            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBachelier):

            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABR):

            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABRShifted):

            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelRatesHW):

            tmat = (capletEndDate - valuationDate) / gDaysInYear
            alpha = self._dayCounter.yearFrac(capletStartDate,
                                              capletEndDate)[0]
            strikePrice = 1.0/(1.0 + alpha * self._strikeRate)
            notionalAdj = (1.0 + self._strikeRate * alpha)
            faceAmount = 1.0
            dfTimes = liborCurve._times
            dfValues = liborCurve._dfs

            v = model.optionOnZCB(texp, tmat, strikePrice, faceAmount,
                                  dfTimes, dfValues)

            # we divide by alpha to offset the multiplication above
            if self._optionType == FinCapFloorTypes.CAP:
                capFloorLetValue = v['put'] * notionalAdj / alpha
            elif self._optionType == FinCapFloorTypes.FLOOR:
                capFloorLetValue = v['call'] * notionalAdj / alpha

        else:
            raise FinError("Unknown model type " + str(model))

        capFloorLetValue *= (self._notional * alpha)

        return capFloorLetValue

###############################################################################

    def printLeg(self):
        ''' Prints the cap floor payment amounts. '''

        print("START DATE:", self._startDate)
        print("MATURITY DATE:", self._maturityDate)
        print("OPTION TYPE", str(self._optionType))
        print("STRIKE (%):", self._strikeRate * 100)
        print("FREQUENCY:", str(self._freqType))
        print("DAY COUNT:", str(self._dayCountType))
        print("VALUATION DATE", self._valuationDate)

        if len(self._capFloorLetValues) == 0:
            print("Caplets not calculated.")
            return

        if self._optionType == FinCapFloorTypes.CAP:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    CAPLET_PV       CUM_PV"
        elif self._optionType == FinCapFloorTypes.FLOOR:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    FLRLET_PV       CUM_PV"

        print(header)

        iFlow = 0

        for paymentDate in self._capFloorLetDates[iFlow:]:
            if iFlow == 0:
                print("%15s %10s %9s %12s %12.6f %12s %12s" %
                      (paymentDate, "-", "-", "-",
                       self._capFloorLetDiscountFactors[iFlow], "-", "-"))
            else:
                print("%15s %10.7f %9.5f %12.2f %12.6f %12.2f %12.2f" %
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
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("STRIKE COUPON", self._strikeRate * 100)
        s += labelToString("OPTION TYPE", str(self._optionType))
        s += labelToString("FREQUENCY", str(self._freqType))
        s += labelToString("DAY COUNT", str(self._dayCountType), "")
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
