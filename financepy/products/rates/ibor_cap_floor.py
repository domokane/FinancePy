##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Implied volatility
# TODO: Term structure of volatility
# TODO: Check that curve anchor date is valuation date ?

from typing import Optional

from ...utils.date import Date
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...utils.error import FinError
from ...utils.schedule import Schedule
from ...utils.helpers import label_to_string, check_argument_types
from ...models.black import FinModelBlack
from ...models.black_shifted import FinModelBlackShifted
from ...models.bachelier import FinModelBachelier
from ...models.sabr import FinModelSABR
from ...models.sabr_shifted import FinModelSABRShifted
from ...models.rates_hull_white_tree import FinModelRatesHW
from ...utils.global_types import FinCapFloorTypes, FinOptionTypes

##########################################################################

from enum import Enum

class FinIborCapFloorModelTypes(Enum):
    BLACK = 1
    SHIFTED_BLACK = 2
    SABR = 3

##########################################################################


class IborCapFloor():
    """ Class for Caps and Floors. These are contracts which observe a Ibor
    reset L on a future start date and then make a payoff at the end of the
    Ibor period which is Max[L-K,0] for a cap and Max[K-L,0] for a floor.
    This is then day count adjusted for the Ibor period and then scaled by
    the contract notional to produce a valuation. A number of models can be
    selected from."""

    def __init__(self,
                 start_date: Date,
                 maturity_date_or_tenor: (Date, str),
                 option_type: FinCapFloorTypes,
                 strikeRate: float,
                 lastFixing: Optional[float] = None,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360_ISDA,
                 notional: float = ONE_MILLION,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Initialise FinIborCapFloor object. """

        check_argument_types(self.__init__, locals())

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            maturity_date = start_date.addTenor(maturity_date_or_tenor)
            calendar = Calendar(self._calendar_type)
            maturity_date = calendar.adjust(maturity_date,
                                           self._bus_day_adjust_type)

        if start_date > maturity_date:
            raise FinError("Start date must be before maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._option_type = option_type
        self._strikeRate = strikeRate
        self._lastFixing = lastFixing
        self._freq_type = freq_type
        self._day_count_type = day_count_type
        self._notional = notional
        self._date_gen_rule_type = date_gen_rule_type

        self._capFloorLetValues = []
        self._capFloorLetAlphas = []
        self._capFloorLetFwdRates = []
        self._capFloorLetIntrinsic = []
        self._capFloorLetDiscountFactors = []
        self._capFloorPV = []

        self._valuation_date = None
        self._day_counter = None

###############################################################################

    def _generateDates(self):

        schedule = Schedule(self._start_date,
                            self._maturity_date,
                            self._freq_type,
                            self._calendar_type,
                            self._bus_day_adjust_type,
                            self._date_gen_rule_type)

        self._capFloorLetDates = schedule._adjusted_dates

##########################################################################

    def value(self, valuation_date, libor_curve, model):
        """ Value the cap or floor using the chosen model which specifies
        the volatility of the Ibor rate to the cap start date. """

        self._valuation_date = valuation_date
        self._generateDates()

        self._day_counter = DayCount(self._day_count_type)
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

        start_date = self._start_date
        end_date = self._capFloorLetDates[1]

        if self._lastFixing is None:
            fwd_rate = libor_curve.fwd_rate(start_date, end_date,
                                         self._day_count_type)
        else:
            fwd_rate = self._lastFixing

        alpha = self._day_counter.year_frac(start_date, end_date)[0]
        df = libor_curve.df(end_date)

        if self._option_type == FinCapFloorTypes.CAP:
            capFloorLetValue = df * alpha * max(fwd_rate - strikeRate, 0.0)
        elif self._option_type == FinCapFloorTypes.FLOOR:
            capFloorLetValue = df * alpha * max(strikeRate - fwd_rate, 0.0)

        capFloorLetValue *= self._notional
        capFloorValue += capFloorLetValue

        self._capFloorLetFwdRates.append(fwd_rate)
        self._capFloorLetValues.append(capFloorLetValue)
        self._capFloorLetAlphas.append(alpha)
        self._capFloorLetIntrinsic.append(capFloorLetValue)
        self._capFloorLetDiscountFactors.append(df)
        self._capFloorPV.append(capFloorValue)

        for i in range(2, numOptions):

            start_date = self._capFloorLetDates[i - 1]
            end_date = self._capFloorLetDates[i]
            alpha = self._day_counter.year_frac(start_date, end_date)[0]

            df = libor_curve.df(end_date)
            fwd_rate = libor_curve.fwd_rate(start_date, end_date,
                                         self._day_count_type)

            if self._option_type == FinCapFloorTypes.CAP:
                intrinsicValue = df * alpha * max(fwd_rate - strikeRate, 0.0)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                intrinsicValue = df * alpha * max(strikeRate - fwd_rate, 0.0)

            intrinsicValue *= self._notional

            capFloorLetValue = self.valueCapletFloorLet(valuation_date,
                                                        start_date,
                                                        end_date,
                                                        libor_curve,
                                                        model)

            capFloorValue += capFloorLetValue

            self._capFloorLetFwdRates.append(fwd_rate)
            self._capFloorLetValues.append(capFloorLetValue)
            self._capFloorLetAlphas.append(alpha)
            self._capFloorLetIntrinsic.append(intrinsicValue)
            self._capFloorLetDiscountFactors.append(df)
            self._capFloorPV.append(capFloorValue)

        return capFloorValue

###############################################################################

    def valueCapletFloorLet(self,
                            valuation_date,
                            capletStartDate,
                            capletEndDate,
                            libor_curve,
                            model):
        """ Value the caplet or floorlet using a specific model. """

        texp = (capletStartDate - self._start_date) / gDaysInYear

        alpha = self._day_counter.year_frac(capletStartDate, capletEndDate)[0]

        f = libor_curve.fwd_rate(capletStartDate, capletEndDate,
                               self._day_count_type)

        k = self._strikeRate
        df = libor_curve.df(capletEndDate)

        if k == 0.0:
            k = 1e-10

        if isinstance(model, FinModelBlack):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBlackShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelBachelier):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABR):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelSABRShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               FinOptionTypes.EUROPEAN_PUT)

        elif isinstance(model, FinModelRatesHW):

            tmat = (capletEndDate - valuation_date) / gDaysInYear
            alpha = self._day_counter.year_frac(capletStartDate,
                                              capletEndDate)[0]
            strike_price = 1.0/(1.0 + alpha * self._strikeRate)
            notionalAdj = (1.0 + self._strikeRate * alpha)
            face_amount = 1.0
            df_times = libor_curve._times
            df_values = libor_curve._dfs

            v = model.optionOnZCB(texp, tmat, strike_price, face_amount,
                                  df_times, df_values)

            # we divide by alpha to offset the multiplication above
            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = v['put'] * notionalAdj / alpha
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = v['call'] * notionalAdj / alpha

        else:
            raise FinError("Unknown model type " + str(model))

        capFloorLetValue *= (self._notional * alpha)

        return capFloorLetValue

###############################################################################

    def printLeg(self):
        """ Prints the cap floor payment amounts. """

        print("START DATE:", self._start_date)
        print("MATURITY DATE:", self._maturity_date)
        print("OPTION TYPE", str(self._option_type))
        print("STRIKE (%):", self._strikeRate * 100)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._day_count_type))
        print("VALUATION DATE", self._valuation_date)

        if len(self._capFloorLetValues) == 0:
            print("Caplets not calculated.")
            return

        if self._option_type == FinCapFloorTypes.CAP:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    CAPLET_PV       CUM_PV"
        elif self._option_type == FinCapFloorTypes.FLOOR:
            header = "PAYMENT_DATE     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    FLRLET_PV       CUM_PV"

        print(header)

        iFlow = 0

        for payment_date in self._capFloorLetDates[iFlow:]:
            if iFlow == 0:
                print("%15s %10s %9s %12s %12.6f %12s %12s" %
                      (payment_date, "-", "-", "-",
                       self._capFloorLetDiscountFactors[iFlow], "-", "-"))
            else:
                print("%15s %10.7f %9.5f %12.2f %12.6f %12.2f %12.2f" %
                      (payment_date,
                       self._capFloorLetAlphas[iFlow],
                       self._capFloorLetFwdRates[iFlow]*100,
                       self._capFloorLetIntrinsic[iFlow],
                       self._capFloorLetDiscountFactors[iFlow],
                       self._capFloorLetValues[iFlow],
                       self._capFloorPV[iFlow]))

            iFlow += 1

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._start_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("STRIKE COUPON", self._strikeRate * 100)
        s += label_to_string("OPTION TYPE", str(self._option_type))
        s += label_to_string("FREQUENCY", str(self._freq_type))
        s += label_to_string("DAY COUNT", str(self._day_count_type), "")
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
