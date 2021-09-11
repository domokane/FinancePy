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
from ...models.black import Black
from ...models.black_shifted import BlackShifted
from ...models.bachelier import Bachelier
from ...models.sabr import SABR
from ...models.sabr_shifted import SABRShifted
from ...models.hw_tree import HWTree
from ...utils.global_types import FinCapFloorTypes, OptionTypes

##########################################################################

from enum import Enum


class IborCapFloorModelTypes(Enum):
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
                 strike_rate: float,
                 lastFixing: Optional[float] = None,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360_ISDA,
                 notional: float = ONE_MILLION,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Initialise IborCapFloor object. """

        check_argument_types(self.__init__, locals())

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            maturity_date = start_date.add_tenor(maturity_date_or_tenor)
            calendar = Calendar(self._calendar_type)
            maturity_date = calendar.adjust(maturity_date,
                                            self._bus_day_adjust_type)

        if start_date > maturity_date:
            raise FinError("Start date must be before maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._option_type = option_type
        self._strike_rate = strike_rate
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

    def _generate_dates(self):

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
        self._generate_dates()

        self._day_counter = DayCount(self._day_count_type)
        num_options = len(self._capFloorLetDates)
        strike_rate = self._strike_rate

        if strike_rate < 0.0:
            raise FinError("Strike < 0.0")

        if num_options <= 1:
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
            capFloorLetValue = df * alpha * max(fwd_rate - strike_rate, 0.0)
        elif self._option_type == FinCapFloorTypes.FLOOR:
            capFloorLetValue = df * alpha * max(strike_rate - fwd_rate, 0.0)

        capFloorLetValue *= self._notional
        capFloorValue += capFloorLetValue

        self._capFloorLetFwdRates.append(fwd_rate)
        self._capFloorLetValues.append(capFloorLetValue)
        self._capFloorLetAlphas.append(alpha)
        self._capFloorLetIntrinsic.append(capFloorLetValue)
        self._capFloorLetDiscountFactors.append(df)
        self._capFloorPV.append(capFloorValue)

        for i in range(2, num_options):

            start_date = self._capFloorLetDates[i - 1]
            end_date = self._capFloorLetDates[i]
            alpha = self._day_counter.year_frac(start_date, end_date)[0]

            df = libor_curve.df(end_date)
            fwd_rate = libor_curve.fwd_rate(start_date, end_date,
                                            self._day_count_type)

            if self._option_type == FinCapFloorTypes.CAP:
                intrinsic_value = df * alpha * max(fwd_rate - strike_rate, 0.0)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                intrinsic_value = df * alpha * max(strike_rate - fwd_rate, 0.0)

            intrinsic_value *= self._notional

            capFloorLetValue = self.value_caplet_floor_let(valuation_date,
                                                           start_date,
                                                           end_date,
                                                           libor_curve,
                                                           model)

            capFloorValue += capFloorLetValue

            self._capFloorLetFwdRates.append(fwd_rate)
            self._capFloorLetValues.append(capFloorLetValue)
            self._capFloorLetAlphas.append(alpha)
            self._capFloorLetIntrinsic.append(intrinsic_value)
            self._capFloorLetDiscountFactors.append(df)
            self._capFloorPV.append(capFloorValue)

        return capFloorValue

###############################################################################

    def value_caplet_floor_let(self,
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

        k = self._strike_rate
        df = libor_curve.df(capletEndDate)

        if k == 0.0:
            k = 1e-10

        if isinstance(model, Black):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, BlackShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, Bachelier):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABR):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABRShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                capFloorLetValue = model.value(f, k, texp, df,
                                               OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, HWTree):

            tmat = (capletEndDate - valuation_date) / gDaysInYear
            alpha = self._day_counter.year_frac(capletStartDate,
                                                capletEndDate)[0]
            strike_price = 1.0/(1.0 + alpha * self._strike_rate)
            notionalAdj = (1.0 + self._strike_rate * alpha)
            face_amount = 1.0
            df_times = libor_curve._times
            df_values = libor_curve._dfs

            v = model.option_on_zcb(texp, tmat, strike_price, face_amount,
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

    def print_leg(self):
        """ Prints the cap floor payment amounts. """

        print("START DATE:", self._start_date)
        print("MATURITY DATE:", self._maturity_date)
        print("OPTION TYPE", str(self._option_type))
        print("STRIKE (%):", self._strike_rate * 100)
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
        s += label_to_string("STRIKE COUPON", self._strike_rate * 100)
        s += label_to_string("OPTION TYPE", str(self._option_type))
        s += label_to_string("FREQUENCY", str(self._freq_type))
        s += label_to_string("DAY COUNT", str(self._day_count_type), "")
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
